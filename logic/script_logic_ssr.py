import numpy as np
from collections import OrderedDict

try:
    singleshotlogic
except NameError:
    manager.startModule('logic', 'singleshotlogic')

awg = sequencegeneratorlogic.pulsegenerator()

####################################### just a SSR readout, no real experiment ###########################

def add_just_ssr_info(experiment, qm_dict):
    qm_dict['experiment'] = experiment
    qm_dict['gated'] = True
    qm_dict['sequence_mode'] = True
    qm_dict['charge_state_selection'] = False
    qm_dict['ssr_counter_rollover'] = False
    # for key in ssr.keys():
    #   if key not in qm_dict.keys():
    for key in ssr:
        if key not in qm_dict:
            qm_dict[key] = ssr[key]
    for key in setup:
        if key not in qm_dict.keys():
            qm_dict[key] = setup[key]
    return qm_dict

def add_ssr_info(experiment, qm_dict):
    qm_dict['experiment'] = experiment
    qm_dict['gated'] = True
    qm_dict['sequence_mode'] = True
    qm_dict['charge_state_selection'] = False
    qm_dict['ssr_counter_rollover'] = True
    for key in ssr:
        if key not in qm_dict:
            qm_dict[key] = ssr[key]
    for key in setup:
        if key not in qm_dict.keys():
            qm_dict[key] = setup[key]
    return qm_dict

def just_ssr_measurement(qm_dict):
    """
    Just-SSR measurements have no sensing or QIP elements - we just have laser initialisation
    and then single shot readout. We can still choose to map the electron state to the nuclear spin,
    using CnNOTe + CeNOTn gates before the readout [ssr_mode=mapping], or just look at nuclear spin flips
    with a CeNOTn gate before readout [ssr_mode = anything else]

    Just-SSR is used for characterising the SSR fidelity, e.g. for monitoring the SSR fidelity to optimise
    the dc magnetic field alignment

    """

    # initalize pulsed GUI
    pulsedmeasurementlogic._initialize_data_arrays()
    if qm_dict['ssr_mode'] == 'mapping':
        pulsedmasterlogic.mapped_state_array = np.zeros(len(qm_dict['params']['controlled_variable']))
    else:
        pulsedmasterlogic.spin_flip_array = np.zeros(len(qm_dict['params']['controlled_variable']))
        pulsedmasterlogic.spin_flip_error = np.zeros(len(qm_dict['params']['controlled_variable']))

    # set parameters
    singleshotlogic.set_sequence_length(qm_dict['sequence_length'])
    qm_dict['num_of_points'] = int(np.floor(qm_dict['measurement_time'] / qm_dict['sequence_length']))
    print('num_of_points = {}'.format(qm_dict['num_of_points']))

    # intialise the SSR GUI and SSR-counter
    set_up_ssr_measurement(qm_dict)

    # do an optical or frequency optimization
    # check_optimization(qm_dict, rr, ii) #LOCALFIX Andrew 10/2/2019: commented out, as frequency optimisation doesn't work atm

    user_terminated = basic_ssr_measurement(qm_dict)
    return


def set_up_ssr_measurement(qm_dict):
    """
    Defines the pulse extraction and pulse analysis settings for both non-normalised and normalised SSR measurements,
    and sets up the fastcounter for SSR operation

    Required parameters:

    PULSE EXTRACTION
    qm_dict['ssr_normalise'] - True/False
    qm_dict['laser_ssr_length'] - length of laser readout pulse during SSR
    qm_dict['laser_delay'] - delay between sending laser TTL trigger and receiving fluorescence counts
    qm_dict['laser_safety'] -
    qm_dict['wait_ssr_length'] - Wait time after SSR laser pulse to allow for NV relaxation from metastable state
    qm_dict['mw_cnot_rabi_period'] - Twice the microwave CNOT gate duration
    qm_dict['bin_width_s'] - bin width (in seconds) of the fast counter

    PULSE ANALYSIS


    COUNTER

    qm_dict['ssr_counter_rollover']


    NOTE: if the normalised SSR readout sequence changes, the timings in this function
    will need to be manually updated
    """

    singleshotlogic.set_parameters(qm_dict)

    #### Set up pulse extraction ###################

    # pulse extraction settings for first (or only) laser pulse
    laser_count_gate_width = qm_dict['laser_ssr_length'] + qm_dict['laser_safety']
    laser_time_rising0 = qm_dict['laser_delay']
    laser_time_falling0 = laser_count_gate_width + laser_time_rising0

    # convert timing into units of fastcounter bins and store in settings dictionary
    extraction_settings_dict = dict()
    extraction_settings_dict['laser_indices_rising0'] = int(laser_time_rising0 / qm_dict['bin_width_s'])
    extraction_settings_dict['laser_indices_falling0'] = int(laser_time_falling0 / qm_dict['bin_width_s'])

    if 'ssr_normalise' in qm_dict and not qm_dict['ssr_normalise']:
        # regular, non-normalised SSR, with alternating CnNOTe and laser pulses
        extraction_settings_dict['num_laser_pulses'] = 1

    elif 'ssr_normalise' in qm_dict and qm_dict['ssr_normalise']:
        # normalised SSR. The CnNOTe gate on alternating nuclear states,
        # and we have 2 laser pulses for each fastcounter 1D timetrace
        extraction_settings_dict['num_laser_pulses'] = 2

        # define timing of second laser pulse relative to first
        laser_time_rising1 = laser_time_falling0 - qm_dict['laser_safety'] + qm_dict['wait_ssr_length'] + \
                             qm_dict['mw_cnot_rabi_period'] / 2 + qm_dict['laser_delay']
        laser_time_falling1 = laser_count_gate_width + laser_time_rising1

        # convert timing into units of fastcounter bins and store in settings dictionary
        extraction_settings_dict['laser_indices_rising1'] = int(laser_time_rising1/qm_dict['bin_width_s'])
        extraction_settings_dict['laser_indices_falling1'] = int(laser_time_falling1/qm_dict['bin_width_s'])

    else:
        print('error: need to define qm_dict[normalise] parameter')
        return cause_an_error

    # send extraction settings to pulsed master logic (and on to pulsed_measurement_logic)
    pulsedmasterlogic.set_extraction_settings({'method': 'absolute_timing',
                                               'extraction_settings': extraction_settings_dict})

    #### Set up pulse analysis ###################

    if 'signal_end' not in qm_dict:
        pulsedmasterlogic.set_analysis_settings({'method': 'sum', 'signal_start': 0.0,
                                             'signal_end': qm_dict['laser_ssr_length']})
    else:
        pulsedmasterlogic.set_analysis_settings({'method': 'sum', 'signal_start': 0.0,
                                                 'signal_end': qm_dict['signal_end']})

    #### Set up fast counter ###################
    ''' 
    countlength = fast counter histogram length, roughly equal to laser pulse duration
    num_of_points = n_steps_controlled_variable, or for just-SSR, meas_time/sequence_length
    '''

    # if 'num_of_tries' in qm_dict:
    #     qm_dict['measurement_time'] = qm_dict['sequence_length'] * qm_dict['num_of_tries']

    if 'ssr_normalise' in qm_dict and not qm_dict['ssr_normalise']:
        # regular, non-normalised SSR, with alternating CnNOTe and laser pulses
        qm_dict['countlength'] = qm_dict['laser_ssr_length'] + qm_dict['laser_delay']
        print('qm_dict[countlength] = {}'.format(qm_dict['countlength']))
        singleshotlogic.set_ssr_counter_settings({'num_of_points': qm_dict['num_of_points'],
                                                  'countlength': qm_dict['countlength'],
                                                  'bin_width_s': qm_dict['bin_width_s'],
                                                  'rollover': qm_dict['ssr_counter_rollover']})

    elif 'ssr_normalise' in qm_dict and qm_dict['ssr_normalise']:
        # normalised SSR. The CnNOTe gate on alternating nuclear states,
        # and we have 2 laser pulses for each fastcounter 1D timetrace
        qm_dict['countlength'] = 2*qm_dict['laser_ssr_length'] + qm_dict['wait_ssr_length'] +\
                                 qm_dict['mw_cnot_rabi_period']/2 + 2*qm_dict['laser_delay']+ 2*qm_dict['laser_safety']

        singleshotlogic.set_ssr_counter_settings({'num_of_points': qm_dict['num_of_points'],
                                                  'countlength': qm_dict['countlength'],
                                                  'bin_width_s': qm_dict['bin_width_s'],
                                                  'rollover': qm_dict['ssr_counter_rollover']})

    else:
        print('error: need to define qm_dict[normalise] parameter')
        return cause_an_error

    # time.sleep(0.1)
    return


def basic_ssr_measurement(qm_dict, rep_index=1, para_index=1):
    singleshotlogic.toggle_ssr_measurement(True)
    time.sleep(0.1)
    while not pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)

    user_terminated = control_measurement(qm_dict)

    singleshotlogic.toggle_ssr_measurement(False)

    while pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
    singleshotlogic.save_measurement(qm_dict['name']+'_'+str(rep_index)+'_'+str(para_index))
    return user_terminated

# def basic_ssr_measurement(qm_dict, rep_index=1, para_index=1):
#     singleshotlogic.toggle_ssr_measurement(True)
#     # pulsedmasterlogic.toggle_pulsed_measurement(True)
#     time.sleep(0.1)
#     while not pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
#     # singleshotlogic.log.error('measurement running')
#     # user_terminated = control_measurement(qm_dict)
#     time.sleep(5)
#     # singleshotlogic.log.error('user_terminated = {}'.format(user_terminated)0)
#     singleshotlogic.toggle_ssr_measurement(False)
#     # pulsedmasterlogic.toggle_pulsed_measurement(False)
#     # if singleshotlogic.ssr_mode == 'flip_probability':
#         # singleshotlogic.log.error('measurement ended, flip_probability = {}'.format(singleshotlogic.spin_flip_prob))
#     # elif singleshotlogic.ssr_mode == 'mapping':
#         # singleshotlogic.log.error('measurement ended, mapped_state = {}'.format(singleshotlogic.mapped_state))
#     # else:
#     #     singleshotlogic.log.error('measurement ended, no ssr_mode')
#     while pulsedmasterlogic.status_dict['measurement_running']: time.sleep(0.5)
#     # singleshotlogic.save_measurement(qm_dict['name']+'_'+str(rep_index)+'_'+str(para_index))
#     return False



###################################### SSR experiment, how it was done conventionally ###########################


def ssr_experiment_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 6, [1])
    return user_terminated

def ssr_experiment_mapping_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 7, [1])
    return user_terminated

def ssr_nuclear_experiment_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 5, [1,2])
    return user_terminated


# def ssr_guide(qm_dict, length, replace):
#
#     # initalize pulsed GUI
#     pulsedmeasurementlogic._initialize_data_arrays()
#     if qm_dict['ssr_mode'] == 'mapping':
#         pulsedmasterlogic.mapped_state_array = np.zeros(len(qm_dict['params']['controlled_variable']))
#     else:
#         pulsedmasterlogic.spin_flip_array = np.zeros(len(qm_dict['params']['controlled_variable']))
#         pulsedmasterlogic.spin_flip_error = np.zeros(len(qm_dict['params']['controlled_variable']))
#     # get the ensemble list
#     ensemble_list = sequencegeneratorlogic.get_sequence(qm_dict['name']).ensemble_list
#     # print('\nensemble_list = {}'.format(ensemble_list))
#     param = sequencegeneratorlogic.get_sequence(qm_dict['name']).sampling_information['step_parameters']
#     # print('param = {}'.format(param))
#     sequence_parameter_list = param[0:length]
#     sequence_length = 0.0
#     # print('ensemble_list = {}'.format(ensemble_list))
#     # print('ensemble_list[0] = {}'.format(ensemble_list[0]))
#     # print('ensemble_list[0][repetitions] = {}'.format(ensemble_list[0]['repetitions']))
#     # compute the length of the sequence for the first data point
#     for ii in range(length):
#         curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[ii]['ensemble'])
#         # print('curr_ensemble = {}'.format(curr_ensemble))
#         sequence_length += \
#             pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[ii]['repetitions']+1)
#     # set parameters
#     singleshotlogic.set_sequence_length(sequence_length)
#     qm_dict['sequence_length'] = sequence_length
#     print('sequence_length = {}'.format(sequence_length))
#
#     # intialise the SSR GUI and SSR-counter
#     set_up_ssr_measurement(qm_dict)
#
#     # # make sequence continuous
#     # sequence_parameter_list[-1][1]['go_to'] = 1
#     # play sequence once only
#     sequence_parameter_list[-1][1]['go_to'] = -1
#
#     # set_up_ssr_measurement(qm_dict)
#     user_terminated = False
#     for rr in range(qm_dict['total_repetitions']):
#         singleshotlogic.log.error('ssr_guide, rr={}'.format(rr))
#         for ii in range(len(qm_dict['params']['controlled_variable'])):
#             singleshotlogic.log.error('ssr_guide, ii={}'.format(ii))
#             # do an optical or frequency optimization
#             # check_optimization(qm_dict, rr, ii) #LOCALFIX Andrew 10/2/2019: commented out, as frequency optimisation doesn't work atm
#             # adapt sequence and sequence length
#             for entry in replace:
#                 sequence_parameter_list[entry] = param[entry + ii * length]
#                 if not rr==0 or not ii==0:
#                     if ii!=0:
#                         sequence_length = \
#                             update_sequence_length(sequence_length, ensemble_list, entry + (ii-1) * length,
#                                                    entry + ii * length)
#                     else:
#                         sequence_length = \
#                             update_sequence_length(sequence_length, ensemble_list, entry + (len(qm_dict['params']['controlled_variable']) - 1) * length,
#                                                    entry + ii * length)
#             qm_dict['sequence_length'] = sequence_length
#             # adapt fastcounter cycles
#             singleshotlogic.set_ssr_counter_settings({'num_of_points': qm_dict['num_of_points'],
#                                                       'countlength': qm_dict['countlength']})
#
#             awg.write_sequence(qm_dict['name'], sequence_parameter_list)
#             awg.load_sequence(qm_dict['name'])
#             # run just_ssr_measurement and break the loop if the user terminated
#             singleshotlogic.log.error('ssr_guide, start basic_ssr_measurement')
#             user_terminated = basic_ssr_measurement(qm_dict, rr, ii)
#             singleshotlogic.log.error('ssr_guide, finished basic_ssr_measurement, user_terminated = {}'.format(user_terminated))
#             if user_terminated: break
#             if qmeas['ssr_mode'] == 'mapping':
#                 update_mapped_state_result(singleshotlogic.mapped_state, rr, ii)
#             else:
#                 update_spin_flip_result(singleshotlogic.spin_flip_prob,
#                               singleshotlogic.spin_flip_error, rr, ii)
#         if user_terminated: break
#
#     # # set_up_ssr_measurement(qm_dict)
#     # user_terminated = False
#     #
#     # qm_dict['sequence_length'] = sequence_length
#     # # adapt fastcounter cycles
#     # singleshotlogic.set_ssr_counter_settings({'num_of_points': qm_dict['num_of_points'],
#     #                                               'countlength': qm_dict['countlength']})
#     # # awg.write_sequence(qm_dict['name'], sequence_parameter_list)
#     # # awg.load_sequence(qm_dict['name'])
#     # # run just_ssr_measurement and break the loop if the user terminated
#     # # singleshotlogic.log.error('ssr_guide, start basic_ssr_measurement')
#     # user_terminated = basic_ssr_measurement(qm_dict, 1, 1)
#     # # singleshotlogic.log.error(
#     # #     'ssr_guide, finished basic_ssr_measurement, user_terminated = {}'.format(user_terminated))
#     # if qmeas['ssr_mode'] == 'mapping':
#     #     update_mapped_state_result(singleshotlogic.mapped_state, 0, 0)
#     # else:
#     #     update_spin_flip_result(singleshotlogic.spin_flip_prob,
#     #                             singleshotlogic.spin_flip_error, 0, 0)
#
#
#
#     time.sleep(0.5)
#     return user_terminated

def ssr_guide(qm_dict, length, replace):

    # initalize pulsed GUI
    pulsedmeasurementlogic._initialize_data_arrays()
    if qm_dict['ssr_mode'] == 'mapping':
        pulsedmasterlogic.mapped_state_array = np.zeros(len(qm_dict['params']['controlled_variable']))
    else:
        pulsedmasterlogic.spin_flip_array = np.zeros(len(qm_dict['params']['controlled_variable']))
        pulsedmasterlogic.spin_flip_error = np.zeros(len(qm_dict['params']['controlled_variable']))

    # sequence_length = length*len(qm_dict['params']['controlled_variable'])
    sequence_length = qm_dict['sequence_length']
    # set parameters
    singleshotlogic.set_sequence_length(sequence_length)
    # qm_dict['sequence_length'] = sequence_length
    print('sequence_length = {}'.format(sequence_length))

    # intialise pulse extraction/analysis and SSR-counter
    set_up_ssr_measurement(qm_dict)

    # todo: should optimise frequency before measurement
    user_terminated = False
    # do an optical or frequency optimization
    # check_optimization(qm_dict, rr, ii) #LOCALFIX Andrew 10/2/2019: commented out, as frequency optimisation doesn't work atm

    user_terminated = basic_ssr_measurement(qm_dict)

    # todo: needs to be an array, not just a single value
    if qmeas['ssr_mode'] == 'mapping':
        update_mapped_state_result(singleshotlogic.mapped_state, 1, 1)
    else:
        update_spin_flip_result(singleshotlogic.spin_flip_prob,
                      singleshotlogic.spin_flip_error, 1, 1)

    time.sleep(0.5)
    return user_terminated

def update_spin_flip_result(spin_flip_tmp, error_spin_flip_tmp, rep, index):
    """
    send data to pulsed measurement gui for analysis
    :param spin_flip_tmp:
    :param error_spin_flip_tmp:
    :param rep:
    :param index:
    :return:
    """

    pulsedmasterlogic.spin_flip_array[index] = \
        (pulsedmasterlogic.spin_flip_array[index] * rep + spin_flip_tmp)/(rep + 1)

    pulsedmasterlogic.spin_flip_error[index] = np.sqrt(
        (rep/(rep+1)*pulsedmasterlogic.spin_flip_error[index])**2 +  (error_spin_flip_tmp/(rep+1))**2)

    pulsedmeasurementlogic.signal_data[1] = pulsedmasterlogic.spin_flip_array
    pulsedmeasurementlogic.measurement_error[1] = pulsedmasterlogic.spin_flip_error
    pulsedmeasurementlogic.sigMeasurementDataUpdated.emit()
    return

def update_mapped_state_result(mapped_state_tmp, rep, index):

    pulsedmasterlogic.mapped_state_array[index] = \
        (pulsedmasterlogic.mapped_state_array[index] * rep + mapped_state_tmp)/(rep + 1)

    pulsedmeasurementlogic.signal_data[1] = pulsedmasterlogic.mapped_state_array
    #pulsedmeasurementlogic.measurement_error[1] = pulsedmasterlogic.spin_flip_error
    pulsedmeasurementlogic.sigMeasurementDataUpdated.emit()
    return


# def update_sequence_length(sequence_length, ensemble_list, remove_index, added_index):
#     # substract the removed length
#     curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[remove_index]['ensemble'])
#     sequence_length -= \
#         pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[remove_index]['repetitions'] + 1)
#
#     # add the new length
#     curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[added_index]['ensemble'])
#     sequence_length += \
#         pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[remove_index]['repetitions'] + 1)
#     singleshotlogic.set_sequence_length(sequence_length)
#     return sequence_length



############################# Optimization function ####################################

def check_optimization(qm_dict, rep, index):
    ''' Checks if an optical or a frequency optimization should be done '''

    length = len(qm_dict['params']['controlled_variable'])
    if 'ssr_optical_optimize_interval' in qm_dict:
        if qm_dict['ssr_optical_optimize_interval'] != 0 and qm_dict['ssr_optical_optimize_interval'] != None:
            if (rep * length + index) % qm_dict['ssr_optical_optimize_interval'] == 0:
                optimize_position()

    if (index == 0 and rep == 0) and qm_dict['ssr_freq_initial_optimize']:
        optimize_ssr_frequency(qm_dict)
    if 'ssr_freq_optimize_interval' in qm_dict:
        if qm_dict['ssr_freq_optimize_interval'] != 0 and qm_dict['ssr_freq_optimize_interval'] != None:
            if not (index == 0 and rep == 0):
                if (rep * length + index) % qm_dict['ssr_freq_optimize_interval'] == 0:
                    optimize_ssr_frequency(qm_dict)
    return 0

def optimize_ssr_frequency(qm_dict):
    # store the current settings
    meas_settings, signal_data, meas_error = store_current_parameters()

    # get name of current measurement waveform on AWG, before overwriting for optimisation measurements
    awg = pulsedmasterlogic.sequencegeneratorlogic().pulsegenerator()
    awg_asset_name = awg.current_loaded_assets['AWG']
    awg_asset_type = awg.current_loaded_assets_type

    # Generate a new dictionary with the measurement parameters
    freq_opt = OrderedDict()
    freq_opt['name'] = 'Opt_freq2'
    for key in qm_dict:
        if key.startswith('ssr_optimize'):
            freq_opt[key.replace('ssr_optimize_', '')] = qm_dict[key]
    if not 'frequency_on_condition' in freq_opt or not freq_opt['frequency_on_condition']:
        # generate sequence, upload it, set the parameters and run optimization experiment
        do_experiment(experiment=freq_opt['experiment'], qm_dict=freq_opt, meas_type=conventional_measurement,
                  meas_info=add_conventional_information,
                  generate_new=freq_opt['generate_new'], save_tag='Freq_opt', load_tag='')
    else:
        do_experiment(experiment=freq_opt['experiment'], qm_dict=freq_opt, meas_type=perform_measurement_on_condition,
                      meas_info=add_conventional_information,
                      generate_new=freq_opt['generate_new'], save_tag='Freq_opt', load_tag='')
    # perform a final fit
    fit_data, fit_result = pulsedmeasurementlogic.do_fit(freq_opt['fit_method'])
    for key in freq_opt['update_parameters']:
        qm_dict[freq_opt['update_parameters'][key]] = fit_result.result_str_dict[key]['value']
    time.sleep(0.5)
    reset_parameters(meas_settings, signal_data, meas_error)
    update_ssr_sequence(qm_dict)
    return fit_result

def update_ssr_sequence(qm_dict):
    tmp_dict = qm_dict.copy()
    tmp_dict['sequence_mode'] = False
    tmp_dict['name'] = 'SSR'
    generate_sample_upload('singleshot_readout', tmp_dict)
    # time.sleep(0.5)
    # pulsedmeasurementlogic.fastcounter().set_gated(True)
    set_up_ssr_measurement(qm_dict)
    return

def store_current_parameters():
    # get the current measurement settings
    meas_settings = pulsedmasterlogic.measurement_settings
    # get the current measurement data
    signal_data = pulsedmeasurementlogic.signal_data
    # get the current measurement error
    meas_error = pulsedmeasurementlogic.measurement_error
    return meas_settings, signal_data, meas_error

def reset_parameters(meas_settings, signal_data, meas_error):
    # reset the measuremetn settings and reset the GUI
    pulsedmasterlogic.set_measurement_settings(meas_settings)
    # reset the data and error
    pulsedmeasurementlogic._initialize_data_arrays()
    pulsedmeasurementlogic.signal_data = signal_data
    pulsedmeasurementlogic.measurement_error = meas_error
    #remove the fit
    pulsedmasterlogic.do_fit('No Fit')
    return








