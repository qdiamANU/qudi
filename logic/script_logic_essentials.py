from collections import OrderedDict
import numpy as np
import time
import datetime
import copy
import matplotlib.pyplot as plt
from core.util import units
from logic.pulsed.pulse_objects import PulseBlock, PulseBlockEnsemble, PulseSequence

try:
    pulsedmasterlogic
except NameError:
    manager.startModule('logic', 'pulsedmasterlogic')
try:
    scannerlogic
except NameError:
    manager.startModule('logic', 'scannerlogic')
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
setup['gated'] = True
setup['sampling_freq'] = pulsedmasterlogic.pulse_generator_settings['sample_rate']
setup['bin_width_s'] = 1.0e-9
setup['wait_time'] = 1.0e-6
setup['laser_delay'] = 100.0e-9
setup['laser_safety'] = 150.0e-9

if setup['gated']:
    setup['sync_channel'] = 'd_ch2'
    setup['gate_channel'] = 'd_ch3'
else:
    setup['sync_channel'] = 'd_ch2'
    setup['gate_channel'] = 'd_ch3'

setup['laser_channel'] = 'd_ch1'

setup['laser_length'] = 3e-6
setup['trigger_length'] = 20.0e-9

setup['channel_amp'] = 4.0
setup['microwave_channel'] = 'a_ch1'
setup['rf_channel'] = 'a_ch2'
setup['rf_amplitude'] = 0.5
setup['optimize_channel'] = '/Dev1/PFI9'

setup['readout_end'] = 0.5e-6

setup['max_tau'] = 1e-3
setup['max_tau_start'] = 1e-3
setup['max_rabi_period'] = 1e-3
setup['min_microwave_frequency'] = 1
setup['max_microwave_amplitude'] = 2.0

setup['measurement_time'] = 600
setup['optimize_time'] = 300
setup['freq_optimize_time'] = None
setup['analysis_interval'] = 3

setup['use_ext_microwave'] = False
setup['ext_microwave_amplitude'] = -30
setup['ext_microwave_frequency'] = 1

############################## Standard function for conventional and SSR measurements #################################


def do_experiment(experiment, qm_dict, meas_type, meas_info, generate_new=True, save_tag='',
                                load_tag=''):
    t1 = time.time()

    # want to only control AWG output amplitude via sine wave amplitudes, so set AWG channel amplitudes to maximum
    # pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator().set_analog_level({'a_ch1': 4.0, 'a_ch2': 4.0})
    pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator().set_analog_level({'a_ch1': 4.0})

    # add information necessary for measurement type
    qm_dict = meas_info(experiment, qm_dict)
    # perform sanity checks
    if experiment not in pulsedmasterlogic.generate_methods:
        pulsedmasterlogic.log.error('Unknown Experiment {0}'.format(experiment))
        return -1
    if not perform_sanity_check(qm_dict):
        pulsedmasterlogic.log.error('Dictionary sanity check failed')
        return -1

    # prepare the measurement by generating and loading the sequence/waveform
    t2 = time.time()
    prepare_qm(experiment, qm_dict, generate_new)
    print('do_experiment: prepare_qm time = {} s'.format((time.time() - t2)))
    # perform measurement
    t3 = time.time()
    user_terminated = perform_measurement(qm_dict=qm_dict, meas_type=meas_type, load_tag=load_tag, save_tag=save_tag)
    print('do_experiment: perform_measurement time = {} s'.format((time.time() - t3)))

    print('do_experiment: total time = {} s'.format((time.time() - t1)))
    return user_terminated

def perform_sanity_check(qm_dict):
    if 'tau' in qm_dict and qm_dict['tau']>setup['max_tau']:
        return False
    if 'tau_start' in qm_dict and qm_dict['tau_start']>setup['max_tau_start']:
        return False
    if 'rabi_period' in qm_dict and qm_dict['rabi_period']>setup['max_rabi_period']:
        return False
    if 'microwave_frequency' in qm_dict and qm_dict['microwave_frequency']<setup['min_microwave_frequency']:
        return False
    if 'microwave_amplitude' in qm_dict and qm_dict['microwave_amplitude']>setup['max_microwave_amplitude']:
        return False
    if 'microwave_amplitude' in qm_dict and qm_dict['microwave_amplitude']>setup['max_microwave_amplitude']:
        return False
    if 'rf_duration' in qm_dict and qm_dict['rf_duration']>1:
        return False
    return True

def add_conventional_information(experiment, qm_dict):
    qm_dict['experiment'] = experiment
    qm_dict['gated'] = True
    qm_dict['sequence_mode'] = False
    return qm_dict

###################################  Upload and set parameters functionality #######################################


def prepare_qm(experiment, qm_dict,  generate_new = True):
    ###### Prepare a quantum measurement by generating the sequence and loading it up to the pulser
    if generate_new:
        generate_sample_upload(experiment, qm_dict)
    else:
        print(qm_dict['name'])
        load_into_channel(qm_dict['name'], sequence_mode=qm_dict['sequence_mode'])
        try:
            qm_dict.update(memory_dict[qm_dict['name']])
        except:
            pulsedmasterlogic.log.error('Experiment parameters are not known. Needs to be generate newly.')
            return cause_an_error

    if not qm_dict['sequence_mode']:
        qm_dict['sequence_length'] = \
            pulsedmasterlogic.get_ensemble_info(pulsedmasterlogic.saved_pulse_block_ensembles[qm_dict['name']])[0]
    else:
        # the sequence length obtained directly from pulsedmasterlogic is the length of the entire AWG playback
        # i.e. it is NOT the length of one step of a measurement
        # need to divide by number of steps to obtain appropriate sequence length in current context
        if 'num_of_points' in qm_dict:
            qm_dict['sequence_length'] = \
                pulsedmasterlogic.get_sequence_info(pulsedmasterlogic.saved_pulse_sequences[qm_dict['name']])[0] \
                / qm_dict['num_of_points']
        else:
            qm_dict['sequence_length'] = \
                pulsedmasterlogic.get_sequence_info(pulsedmasterlogic.saved_pulse_sequences[qm_dict['name']])[0]

    # Set the parameters
    set_parameters(qm_dict)

    return qm_dict

def customise_setup(dictionary):
    #final_dict = dictionary
    for key in setup.keys():
        if key not in dictionary.keys():
            dictionary[key] = setup[key]
    # get a subdictionary with the generation parameters and set them  todo: check how this is different to set_parameters()
    subdict = dict([(key, dictionary.get(key)) for key in pulsedmasterlogic.generation_parameters if key in dictionary])
    pulsedmasterlogic.set_generation_parameters(subdict)
    return dictionary

def generate_sample_upload(experiment, qm_dict):
    t1 = time.time()
    qm_dict = customise_setup(qm_dict)
    print('generate_sample_upload: customise_setup time = {} s'.format((time.time() - t1)))

    if not qm_dict['sequence_mode']:
        # make sure a previous ensemble is deleted
        pulsedmasterlogic.delete_block_ensemble(qm_dict['name'])
        try:
            pulsedmasterlogic.generate_predefined_sequence(experiment, qm_dict.copy())
        except:
            pulsedmasterlogic.log.error('Generation failed')
            return cause_an_error
        time.sleep(0.01)
        # sample the ensemble
        while pulsedmasterlogic.status_dict['predefined_generation_busy']: time.sleep(0.01)
        # print('qm_dict[name] = {}, saved pulse ensembles = {}'.format(qm_dict['name'], pulsedmasterlogic.saved_pulse_block_ensembles))
        if qm_dict['name'] not in pulsedmasterlogic.saved_pulse_block_ensembles: cause_an_error
        pulsedmasterlogic.sample_ensemble(qm_dict['name'], True)
    else:
        if 'exchange_parts' not in qm_dict or qm_dict['exchange_parts']=={}:
            t2 = time.time()
            # make sure a previous sequence is deleted
            pulsedmasterlogic.delete_sequence(qm_dict['name'])
            # generate the ensemble or sequence
            try:
                pulsedmasterlogic.generate_predefined_sequence(experiment, qm_dict.copy())  # generate_predefined_sequence(predefined_sequence_name, kwargs_dict)
            except:
                pulsedmasterlogic.log.error('generation of predefined sequence failed')
                return cause_an_error
            print('generate_sample_upload: generate_predefined_sequence time = {} s'.format((time.time() - t2)))

            # time.sleep(0.01)

            # wait until sequence generation is complete
            t5 = time.time()
            while pulsedmasterlogic.status_dict['predefined_generation_busy']: time.sleep(0.01)
            print('generate_sample_upload: predefined_generation_busy time = {} s'.format((time.time() - t5)))

            # Sanity check: make sure that the desired sequence is now stored in the pulsed master logic memory
            if qm_dict['name'] not in pulsedmasterlogic.saved_pulse_sequences:
                print('sequence generation error: sequence not in pulsedmasterlogic.saved_pulse_sequences')
                return cause_an_error

            t3 = time.time()
            # sample the sequence
            pulsedmasterlogic.sample_sequence(qm_dict['name'], True)
            print('generate_sample_upload: sample_sequence time = {} s'.format((time.time() - t3)))
        else:
            # get the sequence information
            sequence = pulsedmasterlogic.saved_pulse_sequences.get(qm_dict['name'])
            # just generate the replacement ensembles
            for name in qm_dict['exchange_parts']:
                #adapt name of ensemble to be exchanged
                tmp_dict = qm_dict['exchange_parts'][name]['qm_dict'].copy()
                tmp_dict['name'] = name
                pulsedmasterlogic.delete_block_ensemble(name)
                try:
                    pulsedmasterlogic.generate_predefined_sequence(qm_dict['exchange_parts'][name]['experiment'],
                                                                   tmp_dict.copy())
                except:
                    pulsedmasterlogic.log.error('Generation failed')
                    return cause_an_error
                while pulsedmasterlogic.status_dict['predefined_generation_busy']: time.sleep(0.2)
                if name not in pulsedmasterlogic.saved_pulse_block_ensembles: cause_an_error
                pulsedmasterlogic.sample_ensemble(name, True)
                while pulsedmasterlogic.status_dict['sampload_busy']: time.sleep(0.2)
            # load the sequence
            myawg.write_sequence(qm_dict['name'], sequence.sampling_information['step_parameters'])
            # time.sleep(0.2)
            myawg.load_sequence(qm_dict['name'])

    # wait until sequence is finished sampling
    t4 = time.time()
    while pulsedmasterlogic.status_dict['sampload_busy']: time.sleep(0.1)
    print('generate_sample_upload: sampload_busy time = {} s'.format((time.time() - t4)))
    return


def load_into_channel(name, sequence_mode):
    # upload ensemble to pulser channel with sanity check
    if sequence_mode:
        if name in pulsedmasterlogic.sampled_sequences:
            pulsedmasterlogic.load_sequence(name)
        else:
            pulsedmasterlogic.log.error('Sequence not found. Cannot load to channel')
    else:
        if name + '_ch1' in pulsedmasterlogic.sampled_waveforms:
            pulsedmasterlogic.load_ensemble(name)
        else:
            pulsedmasterlogic.log.error('Ensemble not found. Cannot load to channel')
    # wait until ensemble is loaded to channel
    while pulsedmasterlogic.status_dict['loading_busy']: time.sleep(0.1)
    return


def set_parameters(qm_dict):
    ''' AJH 24/2/2019: Sets the pulsed logic measurement parameters, I think?
    I.e. this enables 'invoke settings' in the pulsed GUI?
    '''
    if not qm_dict['sequence_mode']:
        qm_dict['params'] = pulsedmasterlogic.saved_pulse_block_ensembles.get(qm_dict['name']).measurement_information
    else:
        qm_dict['params'] = pulsedmasterlogic.saved_pulse_sequences.get(qm_dict['name']).measurement_information
    pulsedmasterlogic.set_measurement_settings(qm_dict['params'])

    return qm_dict



############################### Perform Measurement Methods ############################################################


def perform_measurement(qm_dict, meas_type, load_tag='', save_tag='', analysis_interval = None, measurement_time=None,
                        optimize_time=None, freq_optimize_time=None):
    # FIXME: add the possibility to load previous data saved under load_tag
    # laser_off(pulser_on = False)
    pulsedmasterlogic.do_fit('No Fit')
    # save parameters into memory_dict
    memory_dict[qm_dict['name']] = qm_dict.copy()
    ###################### Adjust running and refocussing time ##################
    if measurement_time is not None:
        qm_dict['measurement_time'] = measurement_time
    if optimize_time is not None:
        qm_dict['optimize_time'] = optimize_time
    if freq_optimize_time is not None:
        qm_dict['freq_optimize_time'] = freq_optimize_time
    if analysis_interval is not None:
        qm_dict['analysis_interval'] = analysis_interval

    ################ Start and perform the measurement #################
    t2 = time.time()
    user_terminated = meas_type(qm_dict)
    print('perform_measurement: meas_type time = {} s'.format((time.time() - t2)))
    ########################## Save data ###############################
    save_parameters(save_tag=save_tag, save_dict=qm_dict)
    # if if desired

    if 'fit_experiment' in qm_dict and qm_dict['fit_experiment']!= 'No fit':
        fit_data, fit_result = pulsedmeasurementlogic.do_fit(qm_dict['fit_experiment'])
    if meas_type == conventional_measurement:
        pulsedmasterlogic.save_measurement_data(save_tag, True)
        time.sleep(1)
    return user_terminated


def conventional_measurement(qm_dict):
    #set up
    set_up_conventional_measurement(qm_dict)
    if qm_dict['use_ext_microwave']==True:
        #pulsedmasterlogic.set_extraction_settings({'method': 'threshold'})
        pulsedmasterlogic.set_ext_microwave_settings(frequency=qm_dict['ext_microwave_frequency'], power=qm_dict['ext_microwave_power'],
                                                     use_ext_microwave=True)
        # pulsedmasterlogic.toggle_ext_microwave(True)  # todo: AJH 26/2/19 - I think this is redundant?
    # perform measurement
    pulsedmasterlogic.toggle_pulsed_measurement(True)
    while not pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
    user_terminated = control_measurement(qm_dict, analysis_method=None)
    pulsedmasterlogic.toggle_pulsed_measurement(False)
    while pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
    if qm_dict['use_ext_microwave']==True:
        pulsedmasterlogic.toggle_ext_microwave(False)
    return user_terminated


def set_up_conventional_measurement(qm_dict):
    # fixme: commented out change_sweep_mode, set_delay_start, change_save_mode as these functions aren't implemented for TimeTagger
    #pulsedmeasurementlogic.fastcounter().change_sweep_mode(False)
    pulsedmasterlogic.set_fast_counter_settings({'bin_width': qm_dict['bin_width_s'],
                                                 'record_length': qm_dict['params']['counting_length']})
    time.sleep(0.2)

    # pulse extraction settings for first (or only) laser pulse
    laser_count_gate_width = qm_dict['laser_length'] + qm_dict['laser_safety']
    laser_time_rising0 = qm_dict['laser_delay']
    laser_time_falling0 = laser_count_gate_width + laser_time_rising0

    # convert timing into units of fastcounter bins and store in settings dictionary
    extraction_settings_dict = dict()
    extraction_settings_dict['num_laser_pulses'] = 1
    extraction_settings_dict['laser_indices_rising0'] = int(laser_time_rising0 / qm_dict['bin_width_s'])
    extraction_settings_dict['laser_indices_falling0'] = int(laser_time_falling0 / qm_dict['bin_width_s'])

    # # send extraction settings to pulsed master logic (and on to pulsed_measurement_logic)
    # pulsedmasterlogic.set_extraction_settings({'method': 'absolute_timing',
    #                                            'extraction_settings': extraction_settings_dict})


    pulsedmasterlogic.set_extraction_settings(
       {'method': 'conv_deriv', 'delay': setup['laser_delay'], 'safety': setup['laser_safety']})
    # pulsedmasterlogic.set_extraction_settings(
    #     {'method': 'gated_conv_deriv', 'delay': setup['laser_delay'], 'safety': setup['laser_safety']})



    pulsedmasterlogic.set_analysis_settings({'method': 'mean_norm', 'signal_start': 0, 'signal_end': 500e-9,
                                             'norm_start': 1.8e-6, 'norm_end': 2.8e-6})
    return


def control_measurement(qm_dict, analysis_method=None):
    ################# Set the timer and run the measurement #################
    start_time = time.time()
    optimize_real_time = start_time
    freq_optimize_real_time = start_time
    real_update_time = start_time

    # get name of current measurement waveform on AWG, before overwriting for optimisation measurements
    awg = pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator()
    awg_asset_name = awg.current_loaded_assets['AWG']
    awg_asset_type = awg.current_loaded_assets_type

    while True:
        time.sleep(0.1)
        # singleshotlogic.log.error('run time = {}'.format(time.time() - start_time))
        if qm_dict['measurement_time'] is not None:
            if (time.time() - start_time) > qm_dict['measurement_time']:
                user_terminated = False
                break
        if not pulsedmasterlogic.status_dict['measurement_running']:
            user_terminated = True
            break
        ##################### optimize position #######################
        if qm_dict['optimize_time'] is not None:
            if time.time() - optimize_real_time > qm_dict['optimize_time']:
                # release pulser from 'running'
                pulsedmeasurementlogic.pause_pulsed_measurement()
                print('paused measurement time = {}'.format(time.time()))
                additional_time = optimize_position()
                start_time += additional_time
                additional_time = _reload_awg(awg_asset_name, awg_asset_type)
                start_time += additional_time
                pulsedmeasurementlogic.continue_pulsed_measurement()
                optimize_real_time = time.time()
                print('continuing measurement time = {}'.format(time.time()))
        ####################### optimize frequency ####################
        if qm_dict['freq_optimize_time'] is not None:
            if time.time() - freq_optimize_real_time > qm_dict['freq_optimize_time']:
                pulsedmeasurementlogic.pause_pulsed_measurement()
                additional_time = optimize_frequency_during_experiment(opt_dict=optimize_freq_dict, qm_dict=qm_dict)
                start_time += additional_time
                pulsedmeasurementlogic.continue_pulsed_measurement()
                freq_optimize_real_time = time.time()

        ####################### analyze data ######################
        if (analysis_method and qm_dict['update_time']) is not None:
            if time.time() - real_update_time > qm_dict['update_time']:
                analysis_method()
                real_update_time = time.time()
    time.sleep(0.2)
    return user_terminated


def perform_measurement_on_condition(qm_dict):
    # set up
    set_up_conventional_measurement(qm_dict)
    pulsedmasterlogic.break_var=False
    #FIXME: ADD the Options for optical and microwave optimization
    if ('fit_method' or 'threshold_parameter' or 'fit_threshold') not in qm_dict:
        pulsedmasterlogic.log.error('Not enough parameters specified for measurement on condition!')
        cause_an_error
    pulsedmasterlogic.toggle_pulsed_measurement(True)
    while not pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
    # before the first fit wait so that there is actually some data
    time.sleep(pulsedmasterlogic.timer_interval*1.5)
    fit_value = 1e18
    # if the fit_value is below 1e-15, there is something funny
    while fit_value > qm_dict['fit_threshold'] or fit_value < 1e-15:
        if not all(points==0 for points in pulsedmasterlogic.signal_data[1]):
            time.sleep(pulsedmasterlogic.timer_interval)
            fit_data, fit_result = pulsedmeasurementlogic.do_fit(qm_dict['fit_method'])
            fit_value = fit_result.result_str_dict[qm_dict['threshold_parameter']]['value']
            # user can break it
            if not pulsedmasterlogic.status_dict['measurement_running']: break
    # stop the measurement
    pulsedmasterlogic.toggle_pulsed_measurement(False)
    while pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
    #save_parameters(save_tag='Condition', save_dict=qm_dict)
    #return pulsedmasterlogic.fit_container
    return fit_data, fit_result

#########################################  Save methods ############################################

def save_parameters(save_tag='', save_dict=None):

    timestamp = datetime.datetime.now()

    if save_dict is None:
        save_dict = OrderedDict()

    if not save_dict['sequence_mode']:
        sequence = pulsedmasterlogic.saved_pulse_block_ensembles.get(save_dict['name'])
    else:
        sequence = pulsedmasterlogic.saved_pulse_sequences.get(save_dict['name'])
    # get measurement information
    meas_info = sequence.measurement_information
    # get sampling information
    samp_info = sequence.sampling_information

    # get fastcounter settings
    fc_dict = pulsedmasterlogic.fast_counter_settings

    final_dict = {'parameters': save_dict, 'measurementinfo': meas_info,
                  'sampling_info': samp_info, 'fastcounter_settings': fc_dict}

    # some empty data file
    signal_dict = OrderedDict()
    signal_dict['None'] = np.array([0])

    savelogic.save_data(signal_dict, timestamp=timestamp,
                        parameters=final_dict, fmt='%.15e',
                        filepath=savelogic.get_path_for_module('PulsedMeasurement'),
                        filelabel=save_tag+ '_parameters',
                        delimiter='\t', plotfig=None)
    return


# def save_qdyne_measurement(save_tag):
#
#     parameters = qmeas
#     analysis_dict = read_and_analyse_lst_file()
#
#     if not qmeas['show_full']:
#         plt.style.use(savelogic.mpl_qd_style)
#         fig, ax = plt.subplots(1, 1)
#         ax.plot(analysis_dict['ft_freqs'], analysis_dict['ft_result'], 'o', color='blue',
#                 linestyle=':', linewidth=0.5, label='FFT signal')
#         ax.plot(analysis_dict['fit_freqs'], analysis_dict['fit_result'], '-', color='red',
#                 linestyle='-', linewidth=1, label='FFT signal')
#         ax.set_xlabel('freq')
#         ax.set_ylabel('signal')
#         fig.tight_layout()
#
#         parameters['fit parameters see here'] = ''
#         parameters.update(analysis_dict['fit_params'])
#         del analysis_dict['fit_params']
#         savelogic.save_data(analysis_dict, timestamp=qmeas['timestamp'], parameters=parameters, fmt='%.15e',
#                             filepath=savelogic.get_path_for_module('T1_qdyne'), filelabel=save_tag, filetype='npz',
#                             delimiter='\t', plotfig=fig)
#     else:
#         savelogic.save_data(analysis_dict, timestamp=qmeas['timestamp'], parameters=parameters, fmt='%.15e',
#                             filepath=savelogic.get_path_for_module('T1_qdyne'), filelabel=save_tag, filetype='npz',
#                             delimiter='\t', plotfig=None)
#     return








######################################## Position optimize and laser functions #########################################

def optimize_position():
    # FIXME: Add the option to pause pulsed measurement during position optimization
    """

    similar to that used in e.g. magnet logic
    """

    time_start_optimize = time.time()

    curr_pos = scannerlogic.get_position()
    optimizerlogic.start_refocus(curr_pos, caller_tag='script_logic')

    # check just the state of the optimizer
    while optimizerlogic.module_state() != 'idle':
        time.sleep(0.1)

    # use the position to move the scanner
    scannerlogic.set_position('script_logic',
                              x=optimizerlogic.optim_pos_x,
                              y=optimizerlogic.optim_pos_y,
                              z=optimizerlogic.optim_pos_z, a=0.0)

    time_stop_optimize = time.time()
    additional_time = (time_stop_optimize - time_start_optimize)
    return additional_time

def optimize_position_return_position():
    # FIXME: Add the option to pause pulsed measurement during position optimization

    time_start_optimize = time.time()
    #pulsedmeasurementlogic.fast_counter_pause()
    nicard.digital_channel_switch(setup['optimize_channel'], mode=True)
    # perform refocus
    scannerlogic.stop_scanning()
    crosshair_pos = scannerlogic.get_position()
    optimizerlogic.start_refocus(initial_pos=crosshair_pos)
    while optimizerlogic.module_state() == 'idle':
        time.sleep(0.2)
    while optimizerlogic.module_state() != 'idle':
        time.sleep(0.2)
    scannerlogic.set_position('optimizer', x=optimizerlogic.optim_pos_x, y=optimizerlogic.optim_pos_y,
                              z=optimizerlogic.optim_pos_z, a=0.0)
    time.sleep(0.5)
    # switch off laser
    nicard.digital_channel_switch(setup['optimize_channel'], mode=False)
    # pulsedmeasurementlogic.fast_counter_continue()
    time_stop_optimize = time.time()
    # laser_off()
    additional_time = (time_stop_optimize - time_start_optimize)
    return optimizerlogic.optim_pos_x, optimizerlogic.optim_pos_y, optimizerlogic.optim_pos_z

def optimize_poi(poi):
    # FIXME: Add the option to pause pulsed measurement during position optimization
    print(poi)
    # laser_on()
    time_start_optimize = time.time()
    #pulsedmeasurementlogic.fast_counter_pause()
    nicard.digital_channel_switch(setup['optimize_channel'], mode=True)
    # perform refocus
    poimanagerlogic.optimise_poi(poi)
    while optimizerlogic.module_state() == 'idle':
        time.sleep(0.2)
    while optimizerlogic.module_state() != 'idle':
        time.sleep(0.2)
    scannerlogic.set_position('optimizer', x=optimizerlogic.optim_pos_x, y=optimizerlogic.optim_pos_y,
                              z=optimizerlogic.optim_pos_z, a=0.0)
    time.sleep(0.5)
    # switch off laser
    nicard.digital_channel_switch(setup['optimize_channel'], mode=False)
    # pulsedmeasurementlogic.fast_counter_continue()
    time_stop_optimize = time.time()
    # laser_off()
    additional_time = (time_stop_optimize - time_start_optimize)
    return additional_time

def _reload_awg(name, asset_type):
    """
    @param str name: name of AWG asset to reload
    @param str asset_type: waveform or sequence
    @return float additional_time: time taken
    """
    time_start_load = time.time()
    if asset_type == 'waveform':
        pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator().load_waveform(name)
    elif asset_type == 'sequence':
        pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator().load_sequence(name)
    else:
        print('_reload_awg: asset type ({}) not recognised'.format(asset_type))
        cause_an_error
    time_stop_load = time.time()
    additional_time = (time_stop_load - time_start_load)
    return additional_time

######################################## Microwave frequency optimize functions #########################################


def optimize_frequency_during_experiment(opt_dict, qm_dict):
    # FIXME: Add the moment only working for conventional measurements - need to implement data stashing for SSR
    time_start_optimize = time.time()

    if 'freq_optimization_method' not in opt_dict:
        pulsedmasterlogic.log.error('Not frequency optimization method specified. Cannot run optimization')
        return -1

    # stop pulsed measurement and stash raw data
    pulsedmasterlogic.toggle_pulsed_measurement(False, qm_dict['name'])
    while pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.2)

    # set the frequency optimization interval to None
    opt_dict['freq_optimize_time'] = None

    # generate sequence, upload it, set the parameters and run optimization experiment
    do_experiment(experiment=opt_dict['freq_optimization_method'], qm_dict=opt_dict, meas_type=opt_dict['meas_type'],
                  meas_info=add_conventional_information,
                  generate_new=opt_dict['generate_new'], save_tag = opt_dict['save_tag'])

    # perform a final fit
    fit_data, fit_result = pulsedmeasurementlogic.do_fit(opt_dict['optimize_fit_method'])

    # update the specified parameters
    for key in opt_dict['parameters2update']:
        qm_dict[opt_dict['parameters2update'][key]] = fit_result.best_values[key]

    # generate, sample and upload the new sequence
    prepare_qm(experiment=qm_dict['experiment'], qm_dict=qm_dict, generate_new=True)
    pulsedmasterlogic.do_fit('No Fit')

    pulsedmasterlogic.toggle_pulsed_measurement(True, qm_dict['name'])
    while not pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.2)

    return time.time()-time_start_optimize


def optimize_frequency(opt_dict):
    # Generate a new dictionary with the measurement parameters
    if 'mw_optimization_method' not in opt_dict:
        pulsedmasterlogic.log.error('Not frequency optimization method specified. Cannot run optimization')
        return -1

    # generate sequence, upload it, set the parameters and run optimization experiment
    do_experiment(experiment=opt_dict['mw_optimization_method'], qm_dict=opt_dict, meas_type=opt_dict['meas_type'],
                  meas_info=opt_dict['meas_info'], generate_new=opt_dict['optimize_generate_new'],
                  save_tag = opt_dict['save_tag'])
    # perform a final fit
    fit_data, fit_result = pulsedmeasurementlogic.do_fit(opt_dict['optimize_fit_method'])

    # FIXME:
    # generate, sample and upload the new sequence
    return fit_result

def optimize_Rabi(opt_dict):
    # Generate a new dictionary with the measurement parameters
    if 'optimization_method' not in opt_dict:
        pulsedmasterlogic.log.error('No fit method specified. Cannot run optimization')
        return -1

    # generate sequence, upload it, set the parameters and run optimization experiment
    do_experiment(experiment=opt_dict['optimization_method'], qm_dict=opt_dict, meas_type=opt_dict['meas_type'],
                  meas_info=opt_dict['meas_info'], generate_new=opt_dict['generate_new'],
                  save_tag = opt_dict['save_tag'])
    # perform a final fit
    fit_data, fit_result = pulsedmeasurementlogic.do_fit(opt_dict['fit_method'])

    # FIXME:
    # generate, sample and upload the new sequence
    return fit_result

########################################## Qdyne Analysis functions #########################################################

# FIXME:
# def read_and_analyse_lst_file(type='qdyne'):
#     if type == 'qdyne':
#         headerlength = 89
#         read_dict = pulseextractionlogic.ungated_qdyne_conventional(qmeas['full_list_name'], headerlength=headerlength,
#                                                                     read_lines=pulseextractionlogic.read_lines,
#                                                                     lines_to_read=None,
#                                                                     time_trace=pulseextractionlogic.time_trace)
#         pulseextractionlogic.read_lines = read_dict['read_lines']
#         pulseextractionlogic.time_trace = read_dict['time_trace']
#
#         analysis_dict = pulseanalysislogic.analyse_qdyne_fft(time_trace=pulseextractionlogic.time_trace,
#                                                              range_around_peak=40)
#         analysis_dict['time_trace'] = read_dict['time_trace']
#
#     elif type == 'qdyne_ssr':
#         time_trace = extract_ssr_qdyne_time_trace(traceanalysislogic.count_trace)
#         time_trace = time_trace - np.mean(time_trace)
#         analysis_dict = pulseanalysislogic.analyse_qdyne_fft(time_trace=time_trace, range_around_peak=40,
#                                                              meas_type=type)
#         analysis_dict['time_trace'] = time_trace
#
#     if len(analysis_dict['ft_freqs_whole']) > 1:
#         if qmeas['show_full']:
#             qdplotlogic.set_data([analysis_dict['ft_freqs_whole']], [analysis_dict['ft_result_whole']])
#         else:
#             lorentz_fit = fitlogic.make_lorentzian_fit(x_axis=analysis_dict['ft_freqs'], data=analysis_dict['ft_result'],
#                                                        estimator=fitlogic.estimate_lorentzian_peak)
#             model, params = fitlogic.make_lorentzian_model()
#             analysis_dict['fit_freqs'] = np.linspace(np.min(analysis_dict['ft_freqs']), np.max(analysis_dict['ft_freqs']),
#                                                      1000)
#             analysis_dict['fit_result'] = model.eval(x=analysis_dict['fit_freqs'], params=lorentz_fit.params)
#             analysis_dict['fit_params'] = lorentz_fit.best_values
#             qdplotlogic.set_data([analysis_dict['ft_freqs'],analysis_dict['fit_freqs']],
#                                  [analysis_dict['ft_result'],analysis_dict['fit_result']])
#
#     qdplotlogic.set_domain()
#     qdplotlogic.set_range()
#
#     return analysis_dict
#
#
# #FIXME
# def read_and_save_time_trace(save_path, save_path_info=None, start_line=None, lines_to_read=None):
#
#     headerlength = 89
#
#
#     if start_line is None:
#         start_line = 1
#     # read file
#     read_dict = pulseextractionlogic.ungated_qdyne_conventional(qmeas['full_list_name'], headerlength=headerlength,
#                                                                 read_lines=start_line-1,
#                                                                 lines_to_read=lines_to_read)
#
#     # save time trace as .npz file
#     np.savez_compressed(save_path, np.array(read_dict['time_trace']))
#
#     read_lines = read_dict['read_lines']
#     read_whole_file = read_dict['read_whole_file']
#     # save info for the time trace
#     if save_path_info != None:
#         with open(save_path_info, 'a') as infofile:
#             infofile.write(save_path + ': ' + str(read_dict['number_of_sweeps']) + ' sweeps \n')
#             infofile.write('read lines of list file: ' + str(read_lines) + '\n')
#             infofile.write('finished to read whole file ')
#             infofile.write(str(read_whole_file))
#             infofile.write('\n \n')
#
#     return read_whole_file, read_lines







################################## Automized measurements for a list of NV centers #####################################

def do_automized_measurements(qm_dict, autoexp):

    # make sure there is an option to break the loop
    pulsedmasterlogic.break_variable = False

    if qm_dict['optimize_at_start']:
        optimize_position()

    # perform all experiments
    for experiment in autoexp:
        if pulsedmasterlogic.break_variable == True: break
        # perform the measurement
        do_experiment(experiment=autoexp[experiment]['type'], qm_dict=autoexp[experiment],
                      meas_type=autoexp[experiment]['meas_type'], meas_info=autoexp[experiment]['meas_info'],
                      generate_new=autoexp[experiment]['generate_new'],
                      save_tag=autoexp[experiment]['name'])

        # fit and update parameters
        if 'fit_experiment' in autoexp[experiment]:
            if autoexp[experiment]['fit_experiment'] != '':
                fit_data, fit_result = pulsedmeasurementlogic.do_fit(autoexp[experiment]['fit_experiment'])
                #pulsedmasterlogic.do_fit(autoexp[experiment]['fit_experiment'])
               # while pulsedmasterlogic.status_dict['fitting_busy']: time.sleep(0.2)
                #time.sleep(1)
                if 'fit_parameter' in autoexp[experiment]:
                    #fit_dict = pulsedmasterlogic.fit_container.current_fit_result.result_str_dict
                    #fit_para = fit_dict[autoexp[experiment]['fit_parameter']]['value']
                    #fit_para = fit_result.best_values[autoexp[experiment]['fit_parameter']]
                    fit_para = fit_result.result_str_dict[autoexp[experiment]['fit_parameter']]['value']
                    if 'update_parameters' in autoexp[experiment]:
                        for key in autoexp[experiment]['update_parameters']:
                            autoexp[key][autoexp[experiment]['update_parameters'][key]] = fit_para
                    if 'update_parameters1' in autoexp[experiment]:
                        fit_para1 = fit_result.result_str_dict[autoexp[experiment]['fit_parameter1']]['value']
                        for key in autoexp[experiment]['update_parameters1']:
                            autoexp[key][autoexp[experiment]['update_parameters1'][key]] = fit_para1
                    if 'update_parameters2' in autoexp[experiment]:
                        if autoexp[experiment]['fit_parameter2'] == 'RF':

                            LO = 10.65e9
                            Bdc = (LO + 2.87e9 + fit_result.result_str_dict['Position 1']['value'])/28e9
                            gammaN14 = 3.0766e6  # [MHz/T]
                            # fit_para2 = 4.945e6 + Bdc*gammaN14 - 25e3
                            fit_para2 = 4.945e6 + Bdc*gammaN14
                            print('Bdc = {} T, rf start freq = {} MHz'.format(Bdc, fit_para2/1e6))
                        else:
                            fit_para2 = fit_result.result_str_dict[autoexp[experiment]['fit_parameter2']]['value']
                        for key in autoexp[experiment]['update_parameters2']:
                            autoexp[key][autoexp[experiment]['update_parameters2'][key]] = fit_para2
                if 'update_frequency' in autoexp[experiment]:

                    freqs_list = list()

                    # if only a single frequency fitted, then naming convention is different:
                    if 'Position' in fit_result.result_str_dict:
                        # append this frequency to a list
                        freqs_list.append(fit_result.result_str_dict['Position']['value']-5e6)
                        for key in autoexp[experiment]['update_frequency']:
                            autoexp[key][autoexp[experiment]['update_frequency'][key]] = freqs_list

                    # otherwise, first find the number of frequencies fitted:
                    else:
                        ii = 0
                        while 'Position {}'.format(ii) in fit_result.result_str_dict:
                            ii += 1
                        n_frequencies = ii

                        # then append these to a list of frequencies
                        for ii in range(n_frequencies):
                            freqs_list.append(fit_result.result_str_dict['Position {}'.format(ii)]['value'])
                        # update subsequent measurements with the frequencies in freq_list
                        for key in autoexp[experiment]['update_frequency']:
                            autoexp[key][autoexp[experiment]['update_frequency'][key]] = str(freqs_list)



        if qm_dict['optimize_between_experiments']:
            optimize_position()


    return


def do_automized_measurements_poi(qm_dict, autoexp):

    # If there is not list of pois specified, take all pois from the current roi
    if not qm_dict['list_pois']:
        qm_dict['list_pois'] = poimanagerlogic.get_all_pois()
        # remove 'crosshair' and 'sample'
        qm_dict['list_pois'].remove('crosshair')
        qm_dict['list_pois'].remove('sample')

    # make sure there is an option to break the loop
    pulsedmasterlogic.break_variable = False
    # check if for the first poi new sequences should be generated
    first_poi = qm_dict['generate_new']
    # loop over all the pois
    for poi in qm_dict['list_pois']:
        if pulsedmasterlogic.break_variable == True: break

        #todo: save waveform name
        # move to current poi and optimize position
        poimanagerlogic.go_to_poi(poi)
        optimize_poi(poi)
        #todo: reload waveform

        # perform all experiments
        for experiment in autoexp:
            if pulsedmasterlogic.break_variable == True: break
            # perform the measurement
            if first_poi:
                do_experiment(experiment=autoexp[experiment]['type'], qm_dict = autoexp[experiment],
                              meas_type=autoexp[experiment]['meas_type'], meas_info=autoexp[experiment]['meas_info'],
                              generate_new=True, save_tag=autoexp[experiment]['name']+poi)
            else:
                do_experiment(experiment=autoexp[experiment]['type'], qm_dict = autoexp[experiment],
                              meas_type=autoexp[experiment]['meas_type'], meas_info=autoexp[experiment]['meas_info'],
                              generate_new=autoexp[experiment]['generate_new'],
                              save_tag=autoexp[experiment]['name']+poi)

            # fit and update parameters
            if 'fit_experiment' in autoexp[experiment]:
                if autoexp[experiment]['fit_experiment'] != '':
                    fit_data, fit_result = pulsedmeasurementlogic.do_fit(autoexp[experiment]['fit_experiment'])
                    #pulsedmasterlogic.do_fit(autoexp[experiment]['fit_experiment'])
                   # while pulsedmasterlogic.status_dict['fitting_busy']: time.sleep(0.2)
                    #time.sleep(1)
                    if 'fit_parameter' in autoexp[experiment]:
                        #fit_dict = pulsedmasterlogic.fit_container.current_fit_result.result_str_dict
                        #fit_para = fit_dict[autoexp[experiment]['fit_parameter']]['value']
                        #fit_para = fit_result.best_values[autoexp[experiment]['fit_parameter']]
                        fit_para = fit_result.result_str_dict[autoexp[experiment]['fit_parameter']]['value']
                        if 'update_parameters' in autoexp[experiment]:
                            for key in autoexp[experiment]['update_parameters']:
                                autoexp[key][autoexp[experiment]['update_parameters'][key]] = fit_para

            if qm_dict['optimize_between_experiments']:
                optimize_poi(poi)
        first_poi = False
    return



################################## Magnet movement #########################################

def get_magnet_pathway(x_range, x_step, y_range, y_step):
    # check if magnet gui and logic are loaded, otherwise load

    position_before = magnet_logic.get_pos()

    # set up pathway and GUI
    pathway, backmap = magnet_logic._create_2d_pathway('x', x_range, x_step, 'y', y_range, y_step, position_before)
    x_start = backmap[0]['x']
    y_start = backmap[0]['y']
    prepared_graph = magnet_logic._prepare_2d_graph(x_start, x_range, x_step, y_start, y_range, y_step)
    magnet_logic._2D_data_matrix, magnet_logic._2D_axis0_data, magnet_logic._2D_axis1_data = prepared_graph
    magnet_logic._2D_add_data_matrix = np.zeros(shape=np.shape(magnet_logic._2D_data_matrix), dtype=object)
    # handle the stupid crosshair
    magnetgui.roi_magnet.setSize([x_range / 10, y_range / 10], True)
    magnetgui.roi_magnet.setPos([position_before['x'], position_before['y']])

    return pathway, backmap, position_before



#################################### arbitrary sequence generation method ##################

def generate_sequence(name, info_dict, rotating_frame = False):
    created_sequences = list()
    element_list = list()
    # generate the indvidual blocks and ensemple and generate sequence
    for block in info_dict:
        #first generate the blocks and ensemples
        pulsedmasterlogic.generate_predefined_sequence(info_dict[block]['method'],info_dict[block]['meas_dict'])
        while info_dict[block]['meas_dict']['name'] not in pulsedmasterlogic.saved_pulse_block_ensembles: time.sleep(0.2)
        # Add the sequence information:
        seq_para = customize_seq_para(info_dict[block]['seq_para'])
        element_list.append([block,seq_para])

    sequence = PulseSequence(name=name, ensemble_list = element_list, rotating_frame = rotating_frame)

    created_sequences.append(sequence)
    return created_sequences



def customize_seq_para(seq_para_dict):
    if 'event_trigger' not in seq_para_dict:
        seq_para_dict['event_trigger'] = 'OFF'
    if 'event_jump_to' not in seq_para_dict:
        seq_para_dict['event_jump_to'] = 0
    if 'wait_for' not in seq_para_dict:
        seq_para_dict['wait_for'] = 'OFF'
    if 'repetitions' not in seq_para_dict:
        seq_para_dict['repetitions'] = 1
    if 'go_to' not in seq_para_dict:
        seq_para_dict['go_to'] = 0
    return seq_para_dict




############################### old method which are currently not used but might be helpful ######################

# def generate_sample_upload(experiment, qmeas_dict):
#     pulsedmasterlogic.generate_predefined_sequence(experiment, qmeas_dict)
#     pulsedmasterlogic.sample_ensemble(qmeas_dict['name'], True)
#     time_tmp = time.time()
#     tries = 0
#     #tries2 = 0
#     while pulsedmasterlogic.status_dict['sampload_busy']:
#         time.sleep(0.1)
#         # unfortunately the software hangs up from time to time. Check if it is already uploaded
#         if time.time() - time_tmp > 10:
#             if qmeas_dict['name'] + '_ch1' in pulsedmasterlogic.sampled_waveforms:
#                 pulsedmasterlogic.status_dict['sampload_busy'] = False
#             time_tmp = time.time()
#             tries = tries + 1
#             # if the upload did not work, try everything again
#             if tries > 2:
#                 tries= 0
#                 pulsedmasterlogic.status_dict['sampload_busy'] = False
#                 generate_sample_upload(experiment, qmeas_dict)
#     return


# #FIXME:
# def set_fastcounter_general(gated, fc_binwidth, length, delay_start=0, save_mode = 0, counts_per_readout=None,
#                                countlength = None, sequences=None, filename=None):
#
#     pulsedmeasurementlogic.fastcounter().change_sweep_mode(gated, countlength, counts_per_readout)
#
#     if gated:
#         # First reduce length to prevent crashes
#         pulsedmasterlogic.set_fast_counter_settings(bin_width=6.4e-9, record_length=100e-9)
#         time.sleep(0.5)
#         # put cycles and preset again. Otherwise for some reason cycles is reset to 1
#
#         singleshotlogic.set_ssr_counter_settings({'countlength': countlength, 'counts_per_readout': counts_per_readout})
#     else:
#         pulsedmasterlogic.set_fast_counter_settings(bin_width=fc_binwidth, record_length=length)
#     time.sleep(0.05)
#
#     pulsedmeasurementlogic.fastcounter().set_delay_start(delay_start)
#     pulsedmeasurementlogic.fastcounter().change_save_mode(save_mode)
#     if filename is not None:
#         pulsedmeasurementlogic.fastcounter().change_filename(filename)
#     return

