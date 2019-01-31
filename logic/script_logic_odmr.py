from collections import OrderedDict
import numpy as np
import time
import datetime
import pickle
import matplotlib.pyplot as plt
from core.util import units
from logic.pulsed.pulse_objects import PulseBlock, PulseBlockEnsemble, PulseSequence

from logic.odmr_logic import ODMRLogic

try:
    optimizerlogic
except NameError:
    manager.startModule('logic', 'optimizerlogic')
try:
    odmrlogic_lowfield
except NameError:
    manager.startModule('logic', 'odmrlogic_lowfield')
try:
    odmrlogic_highfield
except NameError:
    manager.startModule('logic', 'odmrlogic_highfield')



############################################ Static hardware parameters ################################################

# static hardware parameters:
# these will be overwritten by any parameters defined in qm_dict
setup = OrderedDict()
setup['gated'] = False
setup['sampling_freq'] = pulsedmasterlogic.pulse_generator_settings['sample_rate']
setup['bin_width'] = 3.2e-9
setup['wait_time'] = 1.0e-6
setup['laser_delay'] = 510e-9
setup['laser_safety'] = 200e-9

if setup['gated']:
    setup['sync_channel'] = 'd_ch2'
    setup['gate_channel'] = 'd_ch3'
else:
    setup['sync_channel'] = 'd_ch2'
    setup['gate_channel'] = 'd_ch3'

setup['laser_channel'] = 'd_ch1'

setup['laser_length'] = 3e-6
setup['trigger_length'] = 20e-9

setup['delay_length'] = 450e-9

setup['channel_amp'] = 1.0
setup['microwave_channel'] = 'a_ch1'
setup['optimize_channel'] = '/Dev1/PFI9'

setup['readout_end'] = 0.5e-6

setup['max_tau'] = 1e-3
setup['max_tau_start'] = 1e-3
setup['max_rabi_period'] = 1e-3
setup['min_microwave_frequency'] = 1
setup['max_microwave_amplitude'] = 2

setup['measurement_time'] = 600
setup['optimize_time'] = 300
setup['freq_optimize_time'] = None
setup['analysis_interval'] = 3

setup['use_ext_microwave'] = False
setup['ext_microwave_amplitude'] = -30
setup['ext_microwave_frequency'] = 1

###########################################





def perform_odmr_scan(odmr_dict, save_folder, save_after_each_scan=True, save_tag = 'CWODMRscan'):
    """
    A method to scan through GHz frequency ranges.

    num_scans --> number of ODMR scans.

    odmr_dict['num_scans'] = 1
    odmr_dict['runtime'] = 20
    odmr_dict['measurement_setup'] = 'lowfield'

    odmr_dict['awg_freq_start'] = 100e6
    odmr_dict['awg_freq_stop'] = 300e6
    odmr_dict['awg_freq_step'] = 1e6
    odmr_dict['ext_mw_freq_start'] = 2.7e9
    odmr_dict['ext_mw_freq_list'] = list(ext_mw_freq_start + (freq_stop-freq_start)*np.arange(num_scans))
    odmr_dict['ext_mw_power'] = -10

    """
    start_time = time.time()

    if odmr_dict['measurement_setup'] == 'lowfield':
        odmr_logic_local = odmrlogic_lowfield
    elif odmr_dict['measurement_setup'] == 'highfield':
        odmr_logic_local = odmrlogic_highfield
    elif odmr_dict['measurement_setup'] == 'synthHD':
        odmr_logic_local = odmrlogic_synthHD
    else:
        print('measurement setup not recognised! Must be lowfield or highfield')

    cw_odmr_list = []
    freq_data = []
    count_data = []
    position = dict()
    position['x'] = []
    position['y'] = []
    position['z'] = []
    position['time_exp'] = []
    position['time_real'] = []

    for i in range(0, odmr_dict['num_scans']):

        # optimise and log position
        x, y, z = optimize_position_return_position()
        position['x'].append(x)
        position['y'].append(y)
        position['z'].append(z)
        position['time_exp'].append(time.time() - start_time)
        position['time_real'].append(datetime.datetime.now())

        # perform odmr scan
        odmr_dict['ext_mw_freq'] = odmr_dict['ext_mw_freq_list'][i]
        save_tag_scan = save_tag + '_LOfreq' + str(round(odmr_dict['ext_mw_freq_list'][i] * 1e-6)) + 'MHz'
        scan_output = perform_odmr_measurement(odmr_dict,
                                         odmr_dict['fit_function'],
                                         save_after_meas = save_after_each_scan,
                                         save_tag=save_tag_scan)
        cw_odmr_list.append(scan_output)

    # pre-process odmr data
    for i in range(0, odmr_dict['num_scans']):
        for j in range(len(cw_odmr_list[i][0])):
            freq_data.append(cw_odmr_list[i][0][j] + odmr_dict['ext_mw_freq_list'][i])
            # count_data.append(cw_odmr_list[i][1][0][j] / np.mean(cw_odmr_list[i][1][0]))
            count_data.append(cw_odmr_list[i][1][0][j] / np.mean(cw_odmr_list[i][1][0]))

        cw_odmr_data = np.column_stack((freq_data, count_data))

    # pre-process position tracking
    position['x'] = np.array(position['x'])
    position['y'] = np.array(position['y'])
    position['z'] = np.array(position['z'])
    position['time_exp'] = np.array(position['time_exp'])
    position['time_real'] = np.array(position['time_real'])
    position['x_drift'] = position['x']-position['x'][0]
    position['y_drift'] = position['y'] - position['y'][0]
    position['z_drift'] = position['z'] - position['z'][0]

    # save data
    min_freq_str = str(int((odmr_dict['ext_mw_freq_list'][0]+odmr_dict['awg_freq_start']) * 1e-6))
    max_freq_str = str(int((odmr_dict['ext_mw_freq_list'][-1]+odmr_dict['awg_freq_start']) * 1e-6)) + 'MHz'
    with open(save_folder + save_tag+'_'+min_freq_str + 'to' + max_freq_str + '.pkl', 'wb') as f:
        pickle.dump(cw_odmr_data, f, 0)
        pickle.dump(position, f, 0)
    # todo: automatically save in same folder as data - atm this needs to be manually set

    return cw_odmr_data, position

def perform_odmr_measurement(odmr_dict, fit_function='No Fit', save_after_meas=True, save_tag='CWODMR'):
    """ An independant method, which can be called by a task with the proper input values
        to perform an odmr measurement.

        odmr_dict['runtime'] = 20
        odmr_dict['measurement_setup'] = 'lowfield'

        # fit function options are either 'No Fit', or:
        # ['Lorentzian dip', 'Two Lorentzian dips', 'N14', 'N15', 'Two Gaussian dips']
        fitfunction_options = odmrlogic_highfield.get_fit_functions()
        odmr_dict['fit_function'] = fitfunction_options[1]

        odmr_dict['awg_freq_start'] = 100e6
        odmr_dict['awg_freq_stop'] = 300e6
        odmr_dict['awg_freq_step'] = 1e6
        odmr_dict['ext_mw_freq'] = 2.7e9
        odmr_dict['ext_mw_power'] = 0

    @return
    """
    scan_start_time = time.time()

    if odmr_dict['measurement_setup'] == 'lowfield':
        odmr_logic_local = odmrlogic_lowfield
    elif odmr_dict['measurement_setup'] == 'highfield':
        odmr_logic_local = odmrlogic_highfield
    elif odmr_dict['measurement_setup'] == 'synthHD':
        odmr_logic_local = odmrlogic_synthHD
    else:
        print('measurement setup not recognised! Must be lowfield or highfield')

    print('starting odmr measurement')
    timeout = 30
    start_time = time.time()
    while odmr_logic_local.module_state() != 'idle':
        time.sleep(0.5)
        timeout -= (time.time() - start_time)
        if timeout <= 0:
            odmr_logic_local.log.error('perform_odmr_measurement failed. Logic module was still locked '
                           'and 30 sec timeout has been reached.')
            return {}

    # set all relevant parameter:
    freq_start = odmr_dict['awg_freq_start']
    freq_stop = odmr_dict['awg_freq_stop']
    freq_step = odmr_dict['awg_freq_step']
    ext_mw_power = odmr_dict['ext_mw_power']
    runtime = odmr_dict['runtime']
    odmr_logic_local._mw_device.ext_microwave_frequency = odmr_dict['ext_mw_freq']
    odmr_logic_local.set_sweep_parameters(freq_start, freq_stop, freq_step, ext_mw_power)
    odmr_logic_local.set_runtime(runtime)

    # wait for hardware settings to be implemented
    time.sleep(1)

    # start the scan
    odmr_logic_local.sigStartOdmrScan.emit()

    # wait until the scan has started
    while odmr_logic_local.module_state() != 'locked':
        time.sleep(1)
    # print('Scan has started')

    # do periodic refocussing
    control_odmr_measurement(odmr_dict, odmr_logic_local, scan_start_time)

    # wait until the scan has finished
    while odmr_logic_local.module_state() == 'locked':
        time.sleep(1)
        # print('waiting for scan to finish')
    # print('Scan has finished')

    # Perform fit if requested
    if fit_function != 'No Fit':
        odmr_logic_local.do_fit(fit_function)
        fit_params = odmr_logic_local.fc.current_fit_param
    else:
        fit_params = None

    # Save data if requested
    if save_after_meas:
        odmr_logic_local.save_odmr_data(tag=save_tag)
    # print(fit_params)

    return odmr_logic_local.odmr_plot_x, odmr_logic_local.odmr_plot_y, fit_params

def perform_odmr_scan_scannedlo(odmr_dict, save_folder, save_after_each_scan=True, save_tag = 'CWODMRscan'):
    """
    A method to scan through GHz frequency ranges.

    num_scans --> number of ODMR scans.

    odmr_dict['num_scans'] = 1
    odmr_dict['runtime'] = 20
    odmr_dict['measurement_setup'] = 'lowfield'

    odmr_dict['awg_freq_start'] = 100e6
    odmr_dict['awg_freq_stop'] = 300e6
    odmr_dict['awg_freq_step'] = 1e6
    odmr_dict['ext_mw_freq_start'] = 2.7e9
    odmr_dict['ext_mw_freq_list'] = list(ext_mw_freq_start + (freq_stop-freq_start)*np.arange(num_scans))
    odmr_dict['ext_mw_power'] = -10

    """
    start_time = time.time()

    if odmr_dict['measurement_setup'] == 'lowfield':
        odmr_logic_local = odmrlogic_lowfield
    elif odmr_dict['measurement_setup'] == 'highfield':
        odmr_logic_local = odmrlogic_highfield
    elif odmr_dict['measurement_setup'] == 'synthHD':
        odmr_logic_local = odmrlogic_synthHD
    else:
        print('measurement setup not recognised! Must be lowfield or highfield')

    cw_odmr_list = []
    freq_data = []
    count_data = []
    position = dict()
    position['x'] = []
    position['y'] = []
    position['z'] = []
    position['time_exp'] = []
    position['time_real'] = []

    for i in range(0, odmr_dict['num_scans']):

        # optimise and log position
        x, y, z = optimize_position_return_position()
        position['x'].append(x)
        position['y'].append(y)
        position['z'].append(z)
        position['time_exp'].append(time.time() - start_time)
        position['time_real'].append(datetime.datetime.now())

        # perform odmr scan
        odmr_dict['ext_mw_freq'] = odmr_dict['ext_mw_freq_list'][i]
        save_tag_scan = save_tag + '_LOfreq' + str(round(odmr_dict['ext_mw_freq_list'][i] * 1e-6)) + 'MHz'
        scan_output = perform_odmr_measurement_scannedlo(odmr_dict,
                                         odmr_dict['fit_function'],
                                         save_after_meas = save_after_each_scan,
                                         save_tag=save_tag_scan)
        cw_odmr_list.append(scan_output)

    # pre-process odmr data
    for i in range(0, odmr_dict['num_scans']):
        for j in range(len(cw_odmr_list[i][0])):
            freq_data.append(cw_odmr_list[i][0][j])
            # count_data.append(cw_odmr_list[i][1][0][j] / np.mean(cw_odmr_list[i][1][0]))
            count_data.append(cw_odmr_list[i][1][0][j] / np.mean(cw_odmr_list[i][1][0]))

        cw_odmr_data = np.column_stack((freq_data, count_data))

    # pre-process position tracking
    position['x'] = np.array(position['x'])
    position['y'] = np.array(position['y'])
    position['z'] = np.array(position['z'])
    position['time_exp'] = np.array(position['time_exp'])
    position['time_real'] = np.array(position['time_real'])
    position['x_drift'] = position['x']-position['x'][0]
    position['y_drift'] = position['y'] - position['y'][0]
    position['z_drift'] = position['z'] - position['z'][0]

    # save data
    min_freq_str = str(int(scan_output[0][0] * 1e-6))
    max_freq_str = str(int(scan_output[0][-1] * 1e-6)) + 'MHz'
    with open(save_folder + save_tag+'_'+min_freq_str + 'to' + max_freq_str + '.pkl', 'wb') as f:
        pickle.dump(cw_odmr_data, f, 0)
        pickle.dump(position, f, 0)
    # todo: automatically save in same folder as data - atm this needs to be manually set

    return cw_odmr_data, position

def perform_odmr_measurement_scannedlo(odmr_dict, fit_function='No Fit', save_after_meas=True, save_tag='CWODMR'):
    """ An independant method, which can be called by a task with the proper input values
        to perform an odmr measurement.

        odmr_dict['runtime'] = 20
        odmr_dict['measurement_setup'] = 'lowfield'

        # fit function options are either 'No Fit', or:
        # ['Lorentzian dip', 'Two Lorentzian dips', 'N14', 'N15', 'Two Gaussian dips']
        fitfunction_options = odmrlogic_highfield.get_fit_functions()
        odmr_dict['fit_function'] = fitfunction_options[1]

        odmr_dict['awg_freq_start'] = 100e6
        odmr_dict['awg_freq_stop'] = 300e6
        odmr_dict['awg_freq_step'] = 1e6
        odmr_dict['ext_mw_freq'] = 2.7e9
        odmr_dict['ext_mw_power'] = 0

    @return
    """
    scan_start_time = time.time()

    if odmr_dict['measurement_setup'] == 'lowfield':
        odmr_logic_local = odmrlogic_lowfield
    elif odmr_dict['measurement_setup'] == 'highfield':
        odmr_logic_local = odmrlogic_highfield
    elif odmr_dict['measurement_setup'] == 'synthHD':
        odmr_logic_local = odmrlogic_synthHD
    else:
        print('measurement setup not recognised! Must be lowfield or highfield')

    print('starting odmr measurement')
    timeout = 30
    start_time = time.time()
    while odmr_logic_local.module_state() != 'idle':
        time.sleep(0.5)
        timeout -= (time.time() - start_time)
        if timeout <= 0:
            odmr_logic_local.log.error('perform_odmr_measurement failed. Logic module was still locked '
                           'and 30 sec timeout has been reached.')
            return {}

    # set all relevant parameter:
    freq_start = odmr_dict['ext_mw_freq']
    freq_stop = odmr_dict['ext_mw_freq'] + odmr_dict['awg_freq_range']
    freq_step = odmr_dict['awg_freq_step']
    ext_mw_power = odmr_dict['ext_mw_power']
    runtime = odmr_dict['runtime']
    # odmr_logic_local._mw_device.ext_microwave_frequency = odmr_dict['ext_mw_freq']
    odmr_logic_local.set_sweep_parameters(freq_start, freq_stop, freq_step, ext_mw_power)
    odmr_logic_local.set_runtime(runtime)

    print('freqscan start={}, stop={}, step={}'.format(freq_start, freq_stop, freq_step))
    # wait for hardware settings to be implemented
    time.sleep(1)

    # start the scan
    odmr_logic_local.sigStartOdmrScan.emit()

    # wait until the scan has started
    while odmr_logic_local.module_state() != 'locked':
        time.sleep(1)
    # print('Scan has started')

    # do periodic refocussing
    # control_odmr_measurement(odmr_dict, odmr_logic_local, scan_start_time)

    # wait until the scan has finished
    while odmr_logic_local.module_state() == 'locked':
        time.sleep(1)
        # print('waiting for scan to finish')
    # print('Scan has finished')

    # Perform fit if requested
    if fit_function != 'No Fit':
        odmr_logic_local.do_fit(fit_function)
        fit_params = odmr_logic_local.fc.current_fit_param
    else:
        fit_params = None

    # Save data if requested
    if save_after_meas:
        odmr_logic_local.save_odmr_data(tag=save_tag)
    # print(fit_params)

    return odmr_logic_local.odmr_plot_x, odmr_logic_local.odmr_plot_y, fit_params


def control_odmr_measurement(qm_dict, odmr_logic_local, start_time):
    ################# Set the timer and run the measurement #################
    optimize_real_time = start_time

    awg = odmr_logic_local._mw_device.awg()

    while True:
        time.sleep(2)

        # break if runtime has elapsed since starting scan
        if qm_dict['runtime'] is not None:
            if (time.time() - start_time) > qm_dict['runtime']:
                user_terminated = False
                break
        # break if scan is manually terminated
        if odmr_logic_local.module_state() != 'locked':
            user_terminated = True
            break
        ##################### optimize position #######################
        if qm_dict['optimize_time'] is not None:
            if time.time() - optimize_real_time > qm_dict['optimize_time']:
                # get name of current measurement waveform on AWG, before overwriting for optimisation measurement
                awg_asset_name = awg.current_loaded_assets['AWG']
                awg_asset_type = awg.current_loaded_assets_type
                odmr_logic_local.sigStopOdmrScan.emit()
                print('paused measurement time = {}'.format(time.time()))
                # release pulser from 'running'
                additional_time = optimize_position()
                start_time = start_time + additional_time
                if awg_asset_type == 'waveform':
                    additional_time = _reload_awg_waveform(awg_asset_name)
                elif awg_asset_type == 'sequence':
                    additional_time = _reload_awg_sequence(awg_asset_name)
                else:
                    print('error: incorrect awg asset type!')
                    return cause_an_error
                start_time = start_time + additional_time
                odmr_logic_local.sigContinueOdmrScan.emit()
                optimize_real_time = time.time()
                print('continuing measurement time = {}'.format(time.time()))

    time.sleep(0.2)
    return user_terminated