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
    qm_dict['sequence_mode'] = False
    qm_dict['charge_state_selection'] = False
    # for key in ssr.keys():
    #   if key not in qm_dict.keys():
    for key in ssr:
        if key not in qm_dict:
            qm_dict[key] = ssr[key]
    return qm_dict

def add_ssr_info(experiment, qm_dict):
    qm_dict['experiment'] = experiment
    qm_dict['gated'] = True
    qm_dict['sequence_mode'] = True
    qm_dict['charge_state_selection'] = False
    for key in ssr:
        if key not in qm_dict:
            qm_dict[key] = ssr[key]
    return qm_dict

def just_ssr_measurement(qm_dict):
    print('qm_dict[sequence_length] = {}'.format(qm_dict['sequence_length']))
    qm_dict['sequence_length'] = qm_dict['sequence_length'] * qm_dict['counts_per_readout']
    set_up_ssr_measurement(qm_dict)
    #time.sleep(2)
    basic_ssr_measurement(qm_dict)
    return


def set_up_ssr_measurement(qm_dict):
    #pulsedmeasurementlogic.fastcounter().change_sweep_mode(True)
    # time.sleep(0.01)
    singleshotlogic.set_parameters(qm_dict)
    # pulsedmasterlogic.set_extraction_settings({'method': 'conv_deriv'})
    # pulsedmasterlogic.set_analysis_settings({'method': 'sum', 'signal_start': 0.0,
    #                                          'signal_end': qm_dict['laser_length']})
    if 'ssr_normalise' in qm_dict and qm_dict['ssr_normalise']:
        pulsedmasterlogic.set_extraction_settings({'method': 'conv_deriv2'})
    else:
        pulsedmasterlogic.set_extraction_settings({'method': 'conv_deriv'})
    if not 'signal_end' in qm_dict:
        pulsedmasterlogic.set_analysis_settings({'method': 'sum', 'signal_start': 0.0,
                                             'signal_end': qm_dict['laser_length']})
    else:
        pulsedmasterlogic.set_analysis_settings({'method': 'sum', 'signal_start': 0.0,
                                                 'signal_end': qm_dict['signal_end']})
    if 'num_of_tries' in qm_dict:
        qm_dict['measurement_time'] = qm_dict['sequence_length'] * qm_dict['num_of_tries']
    # if not 'countlength' in qm_dict:
    qm_dict['countlength'] = int(1.03*qm_dict['measurement_time']/qm_dict['sequence_length']) # LOCALFIX Andrew: moved 1.03 inside int()
        # qm_dict['countlength'] = 1.03 * int(qm_dict['measurement_time'] / qm_dict['sequence_length'])
    print('qm_dict[countlength] = {}'.format(qm_dict['countlength']))

    print('qm_dict[counts_per_readout] = {}'.format(qm_dict['counts_per_readout']))
    print('qm_dict[delay_length] = {}'.format(qm_dict['delay_length']))
    if 'ssr_normalise' in qm_dict and qm_dict['ssr_normalise']:
        singleshotlogic.set_ssr_counter_settings({'counts_per_readout': 2*qm_dict['counts_per_readout'],
                                                  'countlength': qm_dict['countlength']})
        time.sleep(0.1) # todo: check if necessary - might just be linked to fastcomtec counter?
        pulsedmasterlogic.set_fast_counter_settings({'bin_width': qm_dict['bin_width'],
                                                     'record_length': qm_dict['params']['counting_length'],
                                                     'number_of_gates': qm_dict['countlength']})
    else:
        singleshotlogic.set_ssr_counter_settings({'counts_per_readout': qm_dict['counts_per_readout'],
                                                  'countlength': qm_dict['countlength']})
        time.sleep(0.1)  # todo: check if necessary - might just be linked to fastcomtec counter?
        pulsedmasterlogic.set_fast_counter_settings({'bin_width': qm_dict['bin_width'],
                                                     'record_length': qm_dict['params']['counting_length'],
                                                     'number_of_gates': qm_dict['countlength']})

    time.sleep(1)
    #pulsedmeasurementlogic.fastcounter().set_delay_start(qm_dict['delay_length'])
    #pulsedmeasurementlogic.fastcounter().change_save_mode(0)
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



###################################### SSR experiment, how it was done conventionally ###########################


def ssr_experiment_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 5, [1])
    return user_terminated

def ssr_experiment_mapping_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 6, [1])
    return user_terminated

def ssr_nuclear_experiment_individual(qm_dict):
    user_terminated = ssr_guide(qm_dict, 4, [1,2])
    return user_terminated


def ssr_guide(qm_dict, length, replace):
    # initalize pulsed GUI
    pulsedmeasurementlogic._initialize_data_arrays()
    if qm_dict['ssr_mode'] == 'mapping':
        pulsedmasterlogic.mapped_state_array = np.zeros(len(qm_dict['params']['controlled_variable']))
    else:
        pulsedmasterlogic.spin_flip_array = np.zeros(len(qm_dict['params']['controlled_variable']))
        pulsedmasterlogic.spin_flip_error = np.zeros(len(qm_dict['params']['controlled_variable']))
    # get the ensemble list
    ensemble_list = sequencegeneratorlogic.get_sequence(qm_dict['name']).ensemble_list
    param = sequencegeneratorlogic.get_sequence(qm_dict['name']).sampling_information['step_parameters']
    sequence_parameter_list = param[0:length]
    sequence_length = 0.0
    # print('ensemble_list = {}'.format(ensemble_list))
    # print('ensemble_list[0] = {}'.format(ensemble_list[0]))
    # print('ensemble_list[0][repetitions] = {}'.format(ensemble_list[0]['repetitions']))
    # compute the length of the sequence for the first data point
    for ii in range(length):
        curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[ii]['ensemble'])
        sequence_length += \
            pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[ii]['repetitions']+1)
    # set parameters
    singleshotlogic.set_sequence_length(sequence_length)
    qm_dict['sequence_length'] = sequence_length
    set_up_ssr_measurement(qm_dict)
    print('sequence_length = {}'.format(sequence_length))
    # make sequence continuous
    sequence_parameter_list[-1][1]['go_to'] = 1
    # set_up_ssr_measurement(qm_dict)
    user_terminated = False
    for rr in range(qm_dict['total_repetitions']):
        for ii in range(len(qm_dict['params']['controlled_variable'])):
            # do an optical or frequency optimization
            check_optimization(qm_dict, rr, ii)
            # adapt sequence and sequence length
            for entry in replace:
                sequence_parameter_list[entry] = param[entry + ii * length]
                if not rr==0 or not ii==0:
                    if ii!=0:
                        sequence_length = \
                            update_sequence_length(sequence_length, ensemble_list, entry + (ii-1) * length,
                                                   entry + ii * length)
                    else:
                        sequence_length = \
                            update_sequence_length(sequence_length, ensemble_list, entry + (len(qm_dict['params']['controlled_variable']) - 1) * length,
                                                   entry + ii * length)
            qm_dict['sequence_length'] = sequence_length
            # adapt fastcounter cycles
            singleshotlogic.set_ssr_counter_settings(
                    {'countlength': 1.1*int(qm_dict['measurement_time']/sequence_length)})
            awg.write_sequence(qm_dict['name'], sequence_parameter_list)
            awg.load_sequence(qm_dict['name'])
            # run just_ssr_measurement and break the loop if the user terminated
            user_terminated = basic_ssr_measurement(qm_dict,rr,ii)
            if user_terminated: break
            if qmeas['ssr_mode'] == 'mapping':
                update_mapped_state_result(singleshotlogic.mapped_state, rr, ii)
            else:
                update_spin_flip_result(singleshotlogic.spin_flip_prob,
                              singleshotlogic.spin_flip_error, rr, ii)
        if user_terminated: break
    time.sleep(0.5)
    return user_terminated


def update_sequence_length(sequence_length, ensemble_list, remove_index, added_index):
    # substract the removed length
    curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[remove_index]['ensemble'])
    sequence_length -= \
        pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[remove_index]['repetitions'] + 1)

    # add the new length
    curr_ensemble = sequencegeneratorlogic.get_ensemble(ensemble_list[added_index]['ensemble'])
    sequence_length += \
        pulsedmasterlogic.get_ensemble_info(curr_ensemble)[0] * (ensemble_list[remove_index]['repetitions'] + 1)
    singleshotlogic.set_sequence_length(sequence_length)
    return sequence_length


def update_spin_flip_result(spin_flip_tmp, error_spin_flip_tmp, rep, index):

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
    # change_fastcounter_mode(False)
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

def change_fastcounter_mode(gated=True):

    # pulsedmeasurementgui._delete_extraction_param_widgets()
    # pulsedmeasurementlogic.fastcounter().change_sweep_mode(gated)
    # pulsedmeasurementlogic.fastcounter().set_gated(gated)
    # pulsedmeasurementlogic.extraction_methods
    # # update the GUI
    # pulsedmeasurementgui._pe.extract_param_method_comboBox.clear()
    # pulsedmeasurementgui._pe.extract_param_method_comboBox.addItems(
    #     list(pulsedmasterlogic.extraction_methods))
    # time.sleep(0.01)
    return








