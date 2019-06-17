# -*- coding: utf-8 -*-

"""
This file contains the Qudi Predefined Methods for sequence generator

Qudi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qudi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qudi. If not, see <http://www.gnu.org/licenses/>.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at <https://github.com/Ulm-IQO/qudi/>
"""


import inspect
import importlib
import numpy as np
from logic.pulsed.pulse_objects import PulseBlock, PulseBlockEnsemble, PulseSequence
#from logic.pulsed.pulse_objects import PredefinedGeneratorBase
from logic.pulsed.predefined_generate_methods.helper_methods_setup3 import HelperMethods


"""
General Pulse Creation Procedure:
=================================
- Create at first each PulseBlockElement object
- add all PulseBlockElement object to a list and combine them to a
  PulseBlock object.
- Create all needed PulseBlock object with that idea, that means
  PulseBlockElement objects which are grouped to PulseBlock objects.
- Create from the PulseBlock objects a PulseBlockEnsemble object.
- If needed and if possible, combine the created PulseBlockEnsemble objects
  to the highest instance together in a PulseSequence object.
"""

class SSRPredefinedGeneratorS3(HelperMethods):
    """

    """


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # import os
        # cwd = os.getcwd()
        # print(cwd, os.listdir(cwd+'\logic\pulsed\predefined_generate_methods'))
        # module_name = r'\logic\pulsed\predefined_generate_methods\helper_methods_setup3.py'
        # print(module_name)
        # mod = importlib.import_module('{0}'.format(module_name))
        # importlib.reload(mod)
        # print(inspect.getmembers(mod,True))
        # helperclass = inspect.getmembers(mod,True)[0]
        # instance = helperclass()
        # self._helper_methods = dict()
        # for method_name, method_ref in inspect.getmembers(instance, inspect.ismethod):
        #     setattr(self, method_name, self.method_ref)
        #     self._helper_methods[method_name] = method_ref
        # print(self._helper_methods)

    non_normalised_safety = 1.4
    normalised_safety = 1.8


    # def generate_ajh_test_sequence(self, name='ajh_test_sequence', tau_start=1.0e-6, tau_step=1.0e-6, num_of_points=2,
    #                       laser_name='laser_wait', laser_length=0.4e-6, wait_length=1.0e-6,
    #                       rf_cnot_name='RF', rf_cnot_freq=2.0e6, rf_cnot_amp=0.2, rf_cnot_duration=3.0e-6, rf_cnot_phase=0,
    #                       ssr_name='SSR', mw_cnot_rabi_period=1.0e-6, mw_cnot_amplitude=2.0, mw_cnot_frequency=5.0e6,
    #                       mw_cnot_phase=0, mw_cnot_amplitude2=1.0, mw_cnot_frequency2=7.5e6,
    #                       mw_cnot_phase2=0, ssr_normalise=False, counts_per_readout=3, sync_gate_name='sync_gate',
    #                       rf_channel='a_ch2'):

    def generate_ajh_test_sequence(self, name='SSR', mw_cnot_rabi_period=2e-6, mw_cnot_amplitude=2.0,
        mw_cnot_frequency=5e6, mw_cnot_phase=0, mw_cnot_amplitude2=2.0,
        mw_cnot_frequency2=7.5e6, mw_cnot_phase2=0, ssr_normalise=False,
        counts_per_readout=3, laser_init_length=3.0e-6, wait_init_length=1.5e-6,
        wait_ssr_length=1.0e-6, n_reps = 100):
        """
        Initialise and perform SSR. Intended use is to create nuclear spin time traces for e.g. magnet alignment
        - currently only works for NV with a single qubit (nitrogen nuclear spin). To adapt for more qubits, would
        need to perform multi-frequency CNOT gate

        :param name:
        :param mw_cnot_rabi_period:
        :param mw_cnot_amplitude:
        :param mw_cnot_frequency:
        :param mw_cnot_phase:
        :param mw_cnot_amplitude2:
        :param mw_cnot_frequency2:
        :param mw_cnot_phase2:◙
        :param ssr_normalise:
        :param counts_per_readout:
        :param laser_init_length:
        :param wait_init_length:
        :param wait_ssr_length:
        :return:
        """

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # single_ssr_readout_time = laser_init_length + wait_init_length + \
        #                           counts_per_readout*(mw_cnot_rabi_period/2 + ssr['laser_length'] + ssr['wait_time'])
        single_ssr_readout_time = 1  # fixme: change to actual time
        time_array = np.array(range(n_reps))*single_ssr_readout_time

        ### create pulse elements ##############################

        # Add the laser initialization
        laser_name = 'laser_wait_init'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_laser_wait(name=laser_name, laser_length=laser_init_length, wait_length=wait_init_length,
                                     trigger=False)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        laser_wait_list = [laser_name, seq_param]

        # Add SSR 'sync' trigger
        sync_gate_name = 'sync_gate'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_trigger(name=sync_gate_name, tau=1e-7, digital_channel=self.sync_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        sync_list = [sync_gate_name, seq_param]

        # Add SSR readout mw/laser pulses and 'gate' trigger
        ssr_name = 'SSR'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_singleshot_readout(name=ssr_name, mw_cnot_rabi_period=mw_cnot_rabi_period,
                                             mw_cnot_amplitude=mw_cnot_amplitude,
                                             mw_cnot_frequency=mw_cnot_frequency,
                                             mw_cnot_phase=mw_cnot_phase,
                                             mw_cnot_amplitude2=mw_cnot_amplitude2,
                                             mw_cnot_frequency2=mw_cnot_frequency2,
                                             mw_cnot_phase2=mw_cnot_phase2,
                                             ssr_normalise=ssr_normalise)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({'repetitions': counts_per_readout - 1})
        ssr_list = [ssr_name, seq_param]

        ### Arrange pulse elements into a Sequence ##############################

        # arrange the individual blocks in the correct order
        element_list = list()
        for ii in range(n_reps):
            element_list.append(laser_wait_list.copy())
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous
        element_list = self._make_sequence_continous(element_list)
        sequence = PulseSequence(name=name, ensemble_list=element_list, rotating_frame=False)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=time_array, units=('s', ''), number_of_lasers=n_reps,
                                       labels=('Time', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences

    def generate_just_ssr(self, name='Just-SSR', mw_cnot_rabi_period=2e-6, mw_cnot_amplitude=2.0,
        mw_cnot_frequency=5e6, mw_cnot_phase=0, mw_cnot_amplitude2=2.0,
        mw_cnot_frequency2=7.5e6, mw_cnot_phase2=0, ssr_normalise=False,
        counts_per_readout=3, laser_init_length=3.0e-6, wait_init_length=1.5e-6,
        wait_ssr_length=1.0e-6, laser_ssr_length=500e-9):
        """
        Initialise and perform SSR. Intended use is to create nuclear spin time traces for e.g. magnet alignment
        - only a single SSR loop written to AWG. Need to loop playback and stop measurement using timer
        - currently only works for NV with a single qubit (nitrogen nuclear spin). To adapt for more qubits, would
        need to perform multi-frequency CNOT gate

        :param name:
        :param mw_cnot_rabi_period:
        :param mw_cnot_amplitude:
        :param mw_cnot_frequency:
        :param mw_cnot_phase:
        :param mw_cnot_amplitude2:
        :param mw_cnot_frequency2:
        :param mw_cnot_phase2:
        :param ssr_normalise:
        :param counts_per_readout:
        :param laser_init_length:
        :param wait_init_length:
        :param wait_ssr_length:
        :return:
        """
        n_reps=1

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        single_ssr_readout_time = 1  # todo: AJH 18/3/2019 - not sure if time_array is used at all
        time_array = np.array(range(n_reps))*single_ssr_readout_time

        ### create pulse elements ##############################

        # Add the laser initialization
        laser_name = 'laser_wait_init'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_laser_wait(name=laser_name, laser_length=laser_init_length, wait_length=wait_init_length,
                                     trigger=False)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        laser_wait_list = [laser_name, seq_param]

        # Add SSR 'sync' trigger
        sync_gate_name = 'sync_gate'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_trigger(name=sync_gate_name, tau=1e-7, digital_channel=self.sync_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        sync_list = [sync_gate_name, seq_param]

        # Add SSR readout mw/laser pulses and 'gate' trigger
        ssr_name = 'SSR'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_singleshot_readout(name=ssr_name, mw_cnot_rabi_period=mw_cnot_rabi_period,
                                             mw_cnot_amplitude=mw_cnot_amplitude,
                                             mw_cnot_frequency=mw_cnot_frequency,
                                             mw_cnot_phase=mw_cnot_phase,
                                             mw_cnot_amplitude2=mw_cnot_amplitude2,
                                             mw_cnot_frequency2=mw_cnot_frequency2,
                                             mw_cnot_phase2=mw_cnot_phase2,
                                             ssr_normalise=ssr_normalise,
                                             wait_ssr_length=wait_ssr_length,
                                             laser_ssr_length=laser_ssr_length)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({'repetitions': counts_per_readout - 1})
        ssr_list = [ssr_name, seq_param]

        ### Arrange pulse elements into a Sequence ##############################

        # arrange the individual blocks in the correct order
        element_list = list()
        for ii in range(n_reps):
            element_list.append(laser_wait_list.copy())
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous
        element_list = self._make_sequence_continous(element_list)
        sequence = PulseSequence(name=para_dict['name'], ensemble_list=element_list, rotating_frame=False)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=time_array, units=('s', ''), number_of_lasers=n_reps,
                                       labels=('Time', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences

    def generate_singleshot_readout(self, name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=0.1,
                                    mw_cnot_frequency=2.8e9, mw_cnot_phase = 0, mw_cnot_amplitude2=0.1,
                                    mw_cnot_frequency2=2.8e9, mw_cnot_phase2=0, ssr_normalise=True,
                                    wait_ssr_length=1.0e-6, laser_ssr_length=500e-9):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        ### prevent granularity problems
        mw_cnot_rabi_period = self._adjust_to_samplingrate(mw_cnot_rabi_period, 4)


        # get mw pi pulse block
        mw_pi_element = self._get_multiple_mw_element(length=mw_cnot_rabi_period/2,
                                                      increment=0.0,
                                                      amps=mw_cnot_amplitude,
                                                      freqs=mw_cnot_frequency,
                                                      phases=mw_cnot_phase)

        readout_element1 = self._get_readout_element(wait_time=wait_ssr_length, length=laser_ssr_length, trigger=True)
        block = PulseBlock(name=name)
        block.append(mw_pi_element)
        block.extend(readout_element1)

        if ssr_normalise:

            readout_element2 = self._get_readout_element(wait_time=wait_ssr_length, length=laser_ssr_length, trigger=False)

            # time_between_trigger = self.laser_length + self.wait_time + self.laser_delay
            # if time_between_trigger > self.laser_length * self.non_normalised_safety:
            #     wait = time_between_trigger - self.laser_length * self.non_normalised_safety
            #     extra_waiting_element = self._get_idle_element(length=wait*1.2, increment=0)
            mw_pi_element2 = self._get_multiple_mw_element(length=mw_cnot_rabi_period/2,
                                                           increment=0.0,
                                                           amps=mw_cnot_amplitude2,
                                                           freqs=mw_cnot_frequency2,
                                                           phases=mw_cnot_phase2)
            # waiting_element = self._get_idle_element(length=self.laser_length + 200e-9, increment=0)

            # if self.laser_length + self.wait_time + self.laser_delay > self.laser_length * self.non_normalised_safety:
            #     block.append(extra_waiting_element)
#
            block.append(mw_pi_element2)
            # block.append(waiting_element)
            block.extend(readout_element2)
        created_blocks.append(block)
        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=list(), controlled_variable = [0],
                                        counting_length = self.laser_length * self.normalised_safety if ssr_normalise
                                                         else laser_ssr_length * self.non_normalised_safety)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)

        return created_blocks, created_ensembles, created_sequences

############################################# SSR experiments ################################################

    def generate_ssr_rabi(self, name='SSR-Rabi', tau_start=1.0e-6, tau_step=1.0e-6, num_of_points=2,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9,
                          rf_cnot_name='RF', rf_cnot_freq=2.0e6, rf_cnot_amp=0.2, rf_cnot_duration=3.0e-6, rf_cnot_phase=0,
                          ssr_name='SSR', mw_cnot_rabi_period=1.0e-6, mw_cnot_amplitude=2.0, mw_cnot_frequency=5.0e6,
                          mw_cnot_phase=0, mw_cnot_amplitude2=1.0, mw_cnot_frequency2=7.5e6,
                          mw_cnot_phase2=0, ssr_normalise=False, counts_per_readout=3, sync_gate_name='sync_gate',
                          rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the Rabi pieces
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        para_list=list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name_tmp, tau=tau,
                                                 microwave_amplitude=self.microwave_amplitude,
                                                 microwave_frequency=self.microwave_frequency, microwave_phase=0.0)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._standard_ssr(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=tau_array, units=('s', ''), number_of_lasers=num_of_points,
                                       labels=('Tau', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences

    def generate_ssr_rabi_mapping(self, name='SSR-Rabi', tau_start=1.0e-9, tau_step=1.0e-9, num_of_points=50,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9,
                          rf_cnot_name='RF', rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_duration=100e-6, rf_cnot_phase=0,
                          ssr_name='SSR', mw_cnot_name='MW-CNOT', mw_cnot_rabi_period=20e-9,
                          mw_cnot_amplitude=1.0, mw_cnot_frequency=2.8e9,
                          mw_cnot_phase=0, mw_cnot_amplitude2=1.0, mw_cnot_frequency2=2.8e9,
                          mw_cnot_phase2=0, ssr_normalise=True, counts_per_readout=1000, sync_gate_name='sync_gate',
                          rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the Rabi pieces
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        para_list = list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name_tmp, tau=tau,
                                                 microwave_amplitude=self.microwave_amplitude,
                                                 microwave_frequency=self.microwave_frequency, microwave_phase=0.0)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._standard_ssr_mapping(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=tau_array, units=('s', ''), number_of_lasers=num_of_points,
                                       labels=('Tau', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences



    def generate_ssr_echo(self, name='SSR-Echo', tau_start=1.0e-9, tau_step=1.0e-9, num_of_points=50,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9,
                         rf_cnot_name='RF', rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_duration=100e-6, rf_cnot_phase=0,
                         ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=0.1, mw_cnot_frequency=2.8e9,
                         mw_cnot_phase=0, mw_cnot_amplitude2=0.1, mw_cnot_frequency2=2.8e9,
                         mw_cnot_phase2=0, ssr_normalise=True, counts_per_readout=1000, sync_gate_name='sync_gate',
                          rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the XY8 pieces
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        para_list=list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_echo_s3(name=name_tmp, tau=tau, microwave_amplitude=self.microwave_amplitude,
                                            microwave_frequency=self.microwave_frequency, rabi_period=self.rabi_period)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._standard_ssr(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=2*tau_array, units=('s', ''), number_of_lasers=num_of_points,
                                       labels=('Tau', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences


    def generate_ssr_xy8(self, name='SSR-XY8', tau_start=1.0e-9, tau_step=1.0e-9, num_of_points=50, xy8N=1, ylast=False,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9,
                         rf_cnot_name='RF', rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_duration=100e-6, rf_cnot_phase=0,
                         ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=0.1, mw_cnot_frequency=2.8e9,
                         mw_cnot_phase=0, mw_cnot_amplitude2=0.1, mw_cnot_frequency2=2.8e9,
                         mw_cnot_phase2=0, ssr_normalise=True, counts_per_readout=1000, sync_gate_name='sync_gate',
                          rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the XY8 pieces
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        para_list=list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_xy8_s3(name = name_tmp, tau=tau, microwave_amplitude=self.microwave_amplitude,
                                            microwave_frequency=self.microwave_frequency, rabi_period=self.rabi_period,
                                            xy8N = xy8N, ylast=ylast)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._standard_ssr(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=tau_array, units=('s', ''), number_of_lasers=num_of_points,
                                       labels=('Tau', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences


    def generate_ssr_xy8_Nsweep(self, name='SSR-XY8_Nsweep', tau=1.0e-6, N_start = 1, N_step=1.0e-9, num_of_points=50,
                                ylast=False, laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9,
                                rf_cnot_name='RF', rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_duration=100e-6,
                                rf_cnot_phase=0, ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=0.1,
                                mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=0.1,
                                mw_cnot_frequency2=2.8e9, mw_cnot_phase2=0, ssr_normalise=True, counts_per_readout=1000,
                                sync_gate_name='sync_gate', rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the XY8 pieces
        N_array = N_start + np.arange(num_of_points) * N_step
        para_list=list()
        for number, xy8N in enumerate(N_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_xy8_s3(name = name_tmp, tau=tau, microwave_amplitude=self.microwave_amplitude,
                                            microwave_frequency=self.microwave_frequency, rabi_period=self.rabi_period,
                                            xy8N = xy8N, ylast=ylast)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._standard_ssr(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=N_array*8*tau, units=('s', ''),
                                       number_of_lasers=num_of_points,
                                       labels=('Interaction time', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences



    ################################# Generate standard SSR sequence ###########################################

    def _standard_ssr(self, created_blocks, created_ensembles, para_list, para_dict):

        created_sequences = list()
        # generate initialization, rf control, and ssr_readout)
        created_blocks_tmp, created_ensembles_tmp, laser_wait_list, rf_list1, rf_list2, sync_list, ssr_list = \
            self._initalize_rf_ssr(laser_name=para_dict['laser_name'], laser_init_length=para_dict['laser_init_length'],
                                   wait_init_length=para_dict['wait_init_length'], rf_cnot_name=para_dict['rf_cnot_name'],
                                   rf_cnot_duration=para_dict['rf_cnot_duration'], rf_cnot_amp=para_dict['rf_cnot_amp'],
                                   rf_cnot_freq=para_dict['rf_cnot_freq'], rf_cnot_phase=para_dict['rf_cnot_phase'],
                                   ssr_name=para_dict['ssr_name'], mw_cnot_rabi_period=para_dict['mw_cnot_rabi_period'],
                                   mw_cnot_amplitude=para_dict['mw_cnot_amplitude'],
                                   mw_cnot_frequency=para_dict['mw_cnot_frequency'],
                                   mw_cnot_phase=para_dict['mw_cnot_phase'],
                                   mw_cnot_amplitude2=para_dict['mw_cnot_amplitude2'],
                                   mw_cnot_frequency2=para_dict['mw_cnot_frequency2'],
                                   mw_cnot_phase2=para_dict['mw_cnot_phase2'],
                                   ssr_normalise=para_dict['ssr_normalise'],
                                   counts_per_readout=para_dict['counts_per_readout'],
                                   sync_gate_name=para_dict['sync_gate_name'],
                                   rf_channel=para_dict['rf_channel'],
                                   wait_ssr_length=para_dict['wait_ssr_length'],
                                   laser_ssr_length=para_dict['laser_ssr_length'])
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp

        # bring the individual blocks in the correct order
        element_list = list()
        for ii in range(len(para_list)):
            element_list.append(laser_wait_list.copy())
            element_list.append(para_list[ii])
            element_list.append(rf_list1.copy())
            element_list.append(rf_list2.copy())
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous
        element_list = self._make_sequence_continous(element_list)
        # print('_standard_ssr: rf_list1 = {}'.format(rf_list1))
        # print('_standard_ssr: rf_list2 = {}'.format(rf_list2))
        # print('_standard_ssr: ssr_list = {}'.format(ssr_list))
        # print('_standard_ssr: element_list = {}'.format(element_list))
        sequence = PulseSequence(name=para_dict['name'], ensemble_list=element_list, rotating_frame=False)

        return created_blocks, created_ensembles, sequence


    def _standard_ssr_mapping(self, created_blocks, created_ensembles, para_list, para_dict):

        created_sequences = list()
        # generate initialization, rf control, and ssr_readout)
        created_blocks_tmp, created_ensembles_tmp, laser_wait_list, mw_cnot_list, rf_list1, rf_list2, sync_list, ssr_list = \
            self._initalize_rf_ssr_mapping(laser_name=para_dict['laser_name'], laser_init_length=para_dict['laser_init_length'],
                                   wait_init_length=para_dict['wait_init_length'], rf_cnot_name=para_dict['rf_cnot_name'],
                                   rf_cnot_duration=para_dict['rf_cnot_duration'], rf_cnot_amp=para_dict['rf_cnot_amp'],
                                   rf_cnot_freq=para_dict['rf_cnot_freq'], rf_cnot_phase=para_dict['rf_cnot_phase'],
                                   ssr_name=para_dict['ssr_name'], mw_cnot_name=para_dict['mw_cnot_name'],
                                   mw_cnot_rabi_period=para_dict['mw_cnot_rabi_period'],
                                   mw_cnot_amplitude=para_dict['mw_cnot_amplitude'],
                                   mw_cnot_frequency=para_dict['mw_cnot_frequency'],
                                   mw_cnot_phase=para_dict['mw_cnot_phase'],
                                   mw_cnot_amplitude2=para_dict['mw_cnot_amplitude2'],
                                   mw_cnot_frequency2=para_dict['mw_cnot_frequency2'],
                                   mw_cnot_phase2=para_dict['mw_cnot_phase2'],
                                   ssr_normalise=para_dict['ssr_normalise'],
                                   counts_per_readout=para_dict['counts_per_readout'],
                                   sync_gate_name=para_dict['sync_gate_name'],
                                   rf_channel=para_dict['rf_channel'],
                                   wait_ssr_length=para_dict['wait_ssr_length'],
                                   laser_ssr_length=para_dict['laser_ssr_length'])

        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp

        # bring the individual blocks in the correct order
        element_list = list()
        for ii in range(len(para_list)):
            element_list.append(laser_wait_list.copy())
            element_list.append(para_list[ii])
            element_list.append(mw_cnot_list.copy())
            element_list.append(rf_list1.copy())
            element_list.append(rf_list2.copy())
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous+
        element_list = self._make_sequence_continous(element_list)
        sequence = PulseSequence(name=para_dict['name'], ensemble_list=element_list, rotating_frame=False)

        return created_blocks, created_ensembles, sequence




    #################################### Nuclear control methods ###################################

    def generate_ssr_nuclear_odmr(self, name='Nuclear-ODMR', freq_start=1.0e6, freq_step=1.0e3, num_of_points=50,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9, initial_pi_pulse=False,
                                  rf_duration=1.0e-3, rf_amp=0.1, rf_phase=0,
                                  ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=1.0,
                                  mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=1.0,
                                  mw_cnot_frequency2=2.8e9,  mw_cnot_phase2=0, ssr_normalise=True,
                                  counts_per_readout=1000, sync_gate_name='sync_gate', rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the RF-Rabi pieces
        freq_array = freq_start + np.arange(num_of_points) * freq_step
        para_list = list()
        for number, freq in enumerate(freq_array):
            # name_tmp = name + '_' + str(number)
            # created_blocks_tmp, created_ensembles_tmp, list1, list2 = \
            #     self._chopped_rf_pulse(name = name_tmp, rf_duration=rf_duration, rf_amp=rf_amp,
            #                            rf_freq=freq, rf_phase=rf_phase, rf_channel=rf_channel)
            # created_blocks += created_blocks_tmp
            # created_ensembles += created_ensembles_tmp
            # para_list.append([list1, list2])
            # # print('\ncreated_blocks_tmp = {}'.format(created_blocks_tmp))
            # # print('created_ensembles_tmp = {}'.format(created_ensembles_tmp))
            # # print('list1 = {}'.format(list1))
            # # print('list2 = {}'.format(list2))

            # LOCALFIX Andrew 6/4/2019: chopped_rf_pulse seems to be causing spurious AWG output, so avoiding use
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name_tmp, tau=rf_duration,
                                              microwave_amplitude=rf_amp,
                                              microwave_frequency=freq, microwave_phase=rf_phase,
                                              channel=rf_channel
                                              )
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._nuclear_manipulation(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=freq_array, units=('Hz', ''), number_of_lasers=num_of_points,
                                       labels=('Frequency', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences

    def generate_ssr_nuclear_conditional_odmr(self, name='Nuclear-ODMR', freq_start=1.0e6, freq_step=1.0e3, num_of_points=50,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9, initial_pi_pulse=False,
                                  rf_duration=1.0e6, rf_amp=0.1, rf_phase=0,
                                  ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=1.0,
                                  mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=1.0,
                                  mw_cnot_frequency2=2.8e9,  mw_cnot_phase2=0, ssr_normalise=True,
                                  counts_per_readout=1000, sync_gate_name='sync_gate', rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the RF-Rabi pieces
        freq_array = freq_start + np.arange(num_of_points) * freq_step
        para_list = list()
        for number, freq in enumerate(freq_array):

            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name_tmp, tau=rf_duration,
                                              microwave_amplitude=rf_amp,
                                              microwave_frequency=freq, microwave_phase=rf_phase,
                                              channel=rf_channel
                                              )
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            seq_param = self._customize_seq_para({})
            para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._nuclear_manipulation_conditional(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=freq_array, units=('Hz', ''), number_of_lasers=num_of_points,
                                       labels=('Frequency', 'Spin flip probability'),
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)
        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences

    def generate_ssr_nuclear_rabi(self, name='Nuclear-Rabi', tau_start=1.0e-9, tau_step=1.0e-9, num_of_points=50,
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9, initial_pi_pulse=False,
                                  rf_freq=1.0e6, rf_amp=0.1, rf_phase=0,
                                  ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=1.0,
                                  mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=1.0,
                                  mw_cnot_frequency2=2.8e9, mw_cnot_phase2=0, ssr_normalise=True,
                                  counts_per_readout=1000, sync_gate_name='sync_gate', rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict = locals()

        # generate the RF-Rabi pieces
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        para_list=list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, list1, list2 = \
                self._chopped_rf_pulse(name = name_tmp, rf_duration=tau, rf_amp=rf_amp,
                                       rf_freq=rf_freq, rf_phase=rf_phase, rf_channel=rf_channel)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            para_list.append([list1, list2])

            # LOCALFIX Andrew 6/4/2019: chopped_rf_pulse seems to be causing spurious AWG output, so avoiding use
            # name_tmp = name + '_' + str(number)
            # created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            #     self.generate_single_mw_pulse(name=name_tmp, tau=tau,
            #                                   microwave_amplitude=rf_amp,
            #                                   microwave_frequency=rf_freq, microwave_phase=rf_phase,
            #                                   channel=rf_channel
            #                                   )
            # created_blocks += created_blocks_tmp
            # created_ensembles += created_ensembles_tmp
            # seq_param = self._customize_seq_para({})
            # para_list.append([name_tmp, seq_param])

        created_blocks, created_ensembles, sequence = \
            self._nuclear_manipulation(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating = False, laser_ignore_list = list(),
                                    controlled_variable = tau_array, units=('s', ''), number_of_lasers = num_of_points,
                                    counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)

        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences


    def generate_ssr_contrast(self, name='SSR-contrast',
                          laser_name='laser_wait', laser_init_length=3.0e-6, wait_init_length=1.0e-6,
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9, initial_pi_pulse=False,
                                  rf_cnot_duration=100e-6, rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_phase=0,
                                  ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=1.0,
                                  mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=1.0,
                                  mw_cnot_frequency2=2.8e9, mw_cnot_phase2=0, ssr_normalise=True,
                                  counts_per_readout=1000, sync_gate_name='sync_gate', rf_channel='a_ch2'):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        para_dict=locals()

        # generate the RF-Rabi pieces
        tau_array = [0, rf_cnot_duration]
        para_list=list()
        for number, tau in enumerate(tau_array):
            name_tmp = name + '_' + str(number)
            created_blocks_tmp, created_ensembles_tmp, list1, list2 = \
                self._chopped_rf_pulse(name = name_tmp, rf_duration=tau, rf_amp=rf_cnot_amp,
                                       rf_freq=rf_cnot_freq, rf_phase=rf_cnot_phase, rf_channel=rf_channel)
            created_blocks += created_blocks_tmp
            created_ensembles += created_ensembles_tmp
            para_list.append([list1, list2])

        created_blocks, created_ensembles, sequence = \
            self._nuclear_manipulation(created_blocks, created_ensembles, para_list, para_dict)

        self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                       controlled_variable=tau_array, units=('s', ''), number_of_lasers=2,
                                       counting_length=laser_ssr_length * self.normalised_safety if ssr_normalise
                                       else laser_ssr_length * self.non_normalised_safety)

        created_sequences.append(sequence)
        return created_blocks, created_ensembles, created_sequences


############################################# Standard nuclear manipulstion #####################################

    def _nuclear_manipulation(self, created_blocks, created_ensembles, para_list, para_dict):

        # generate initialization, rf control, and ssr_readout)
        created_blocks_tmp, created_ensembles_tmp, laser_wait_list, sync_list, ssr_list = \
            self._initialize_ssr(laser_name=para_dict['laser_name'], laser_init_length=para_dict['laser_init_length'],
                                   wait_init_length=para_dict['wait_init_length'],
                                   initial_pi_pulse=para_dict['initial_pi_pulse'],
                                   ssr_name=para_dict['ssr_name'],
                                   mw_cnot_rabi_period=para_dict['mw_cnot_rabi_period'],
                                   mw_cnot_amplitude=para_dict['mw_cnot_amplitude'],
                                   mw_cnot_frequency=para_dict['mw_cnot_frequency'],
                                   mw_cnot_phase=para_dict['mw_cnot_phase'],
                                   mw_cnot_amplitude2=para_dict['mw_cnot_amplitude2'],
                                   mw_cnot_frequency2=para_dict['mw_cnot_frequency2'],
                                   mw_cnot_phase2=para_dict['mw_cnot_phase2'],
                                   ssr_normalise=para_dict['ssr_normalise'],
                                   counts_per_readout=para_dict['counts_per_readout'],
                                   wait_ssr_length=para_dict['wait_ssr_length'],
                                   laser_ssr_length=para_dict['laser_ssr_length'])
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp

        # LOCALFIX Andrew 6/4/2019: chopped_rf_pulse seems to be causing spurious AWG output, so avoiding use
        # # bring the individual blocks in the correct order
        # element_list = list()
        # for ii in range(len(para_list)):
        #     element_list.append(laser_wait_list.copy())
        #     element_list.append(para_list[ii][0])
        #     element_list.append(para_list[ii][1])
        #     element_list.append(sync_list.copy())
        #     element_list.append(ssr_list.copy())
        # # make sequence continuous
        # element_list = self._make_sequence_continous(element_list)

        # bring the individual blocks in the correct order
        element_list = list()
        for ii in range(len(para_list)):
            element_list.append(laser_wait_list.copy())
            element_list.append(para_list[ii])
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous
        element_list = self._make_sequence_continous(element_list)

        sequence = PulseSequence(name=para_dict['name'], ensemble_list=element_list, rotating_frame=False)


        return created_blocks, created_ensembles, sequence

    def _nuclear_manipulation_conditional(self, created_blocks, created_ensembles, para_list, para_dict):

        # generate initialization, rf control, and ssr_readout)
        created_blocks_tmp, created_ensembles_tmp, laser_wait_list, sync_list, ssr_list = \
            self._initialize_ssr(laser_name=para_dict['laser_name'], laser_init_length=para_dict['laser_init_length'],
                                   wait_init_length=para_dict['wait_init_length'],
                                   initial_pi_pulse=para_dict['initial_pi_pulse'],
                                   ssr_name=para_dict['ssr_name'],
                                   mw_cnot_rabi_period=para_dict['mw_cnot_rabi_period'],
                                   mw_cnot_amplitude=para_dict['mw_cnot_amplitude'],
                                   mw_cnot_frequency=para_dict['mw_cnot_frequency'],
                                   mw_cnot_phase=para_dict['mw_cnot_phase'],
                                   mw_cnot_amplitude2=para_dict['mw_cnot_amplitude2'],
                                   mw_cnot_frequency2=para_dict['mw_cnot_frequency2'],
                                   mw_cnot_phase2=para_dict['mw_cnot_phase2'],
                                   ssr_normalise=para_dict['ssr_normalise'],
                                   counts_per_readout=para_dict['counts_per_readout'],
                                   wait_ssr_length=para_dict['wait_ssr_length'],
                                   laser_ssr_length=para_dict['laser_ssr_length'])
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp

        # Add electron CNOT
        # ToDo; Add multiple frequency element
        mw_cnot_name = 'mw_cnot'
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_multiple_mw_pulse(name=mw_cnot_name, tau=para_dict['mw_cnot_rabi_period'] / 2.0,
                                            microwave_amplitude=para_dict['mw_cnot_amplitude'],
                                            microwave_frequency=para_dict['mw_cnot_frequency'],
                                            microwave_phase=para_dict['mw_cnot_phase'])
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        mw_cnot_list = [mw_cnot_name, seq_param]

        # bring the individual blocks in the correct order
        element_list = list()
        for ii in range(len(para_list)):
            element_list.append(laser_wait_list.copy())
            element_list.append(mw_cnot_list.copy())
            element_list.append(para_list[ii])
            element_list.append(sync_list.copy())
            element_list.append(ssr_list.copy())
        # make sequence continuous
        element_list = self._make_sequence_continous(element_list)

        sequence = PulseSequence(name=para_dict['name'], ensemble_list=element_list, rotating_frame=False)


        return created_blocks, created_ensembles, sequence


############################################ Helper methods ##################################################

    def _initalize_rf_ssr(self, laser_name, laser_init_length, wait_init_length, rf_cnot_name, rf_cnot_duration,
                          rf_cnot_amp, rf_cnot_freq, rf_cnot_phase,
                          ssr_name, mw_cnot_rabi_period, mw_cnot_amplitude, mw_cnot_frequency, mw_cnot_phase,
                          mw_cnot_amplitude2, mw_cnot_frequency2, mw_cnot_phase2, ssr_normalise,
                          counts_per_readout, sync_gate_name='sync_gate', rf_channel='a_ch2',
                          wait_ssr_length=1.0e-6, laser_ssr_length=500e-9):

        created_blocks = list()
        created_ensembles = list()

        # Add the laser initialization
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_laser_wait(name=laser_name, laser_length=laser_init_length,
                                     wait_length=wait_init_length, trigger=False)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        laser_wait_list = [laser_name, seq_param]

        # Add RF pulse
        created_blocks_tmp, created_ensembles_tmp, rf_list1, rf_list2 = \
            self._chopped_rf_pulse(name=rf_cnot_name, rf_duration=rf_cnot_duration, rf_amp=rf_cnot_amp,
                                   rf_freq=rf_cnot_freq, rf_phase=rf_cnot_phase, rf_channel=rf_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        # # LOCALFIX Andrew 6/4/2019: chopped_rf_pulse seems to be causing spurious AWG output, so avoiding use
        # name_tmp = name + '_' + str(number)
        # created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
        #     self.generate_single_mw_pulse(name=name_tmp, tau=tau,
        #                                   microwave_amplitude=rf_amp,
        #                                   microwave_frequency=rf_freq, microwave_phase=rf_phase,
        #                                   channel=rf_channel
        #                                   )
        # created_blocks += created_blocks_tmp
        # created_ensembles += created_ensembles_tmp
        # seq_param = self._customize_seq_para({})
        # para_list.append([name_tmp, seq_param])

        # Add SSR 'sync' trigger
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_trigger(name=sync_gate_name, tau=1e-7, digital_channel=self.sync_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        sync_list = [sync_gate_name, seq_param]

        # Add SSR readout mw/laser pulses and 'gate' trigger
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_singleshot_readout(name=ssr_name, mw_cnot_rabi_period=mw_cnot_rabi_period,
                                         mw_cnot_amplitude=mw_cnot_amplitude,
                                         mw_cnot_frequency=mw_cnot_frequency,
                                         mw_cnot_phase=mw_cnot_phase,
                                         mw_cnot_amplitude2=mw_cnot_amplitude2,
                                         mw_cnot_frequency2=mw_cnot_frequency2,
                                         mw_cnot_phase2=mw_cnot_phase2,
                                         ssr_normalise=ssr_normalise,
                                         wait_ssr_length=wait_ssr_length,
                                         laser_ssr_length=laser_ssr_length)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({'repetitions': counts_per_readout-1})
        ssr_list = [ssr_name, seq_param]

        return created_blocks, created_ensembles, laser_wait_list, rf_list1, rf_list2, sync_list, ssr_list

    def _initalize_rf_ssr_mapping(self, laser_name, laser_init_length, wait_init_length, rf_cnot_name, rf_cnot_duration,
                                  rf_cnot_amp, rf_cnot_freq, rf_cnot_phase,
                                  ssr_name, mw_cnot_name, mw_cnot_rabi_period, mw_cnot_amplitude,
                                  mw_cnot_frequency, mw_cnot_phase,
                                  mw_cnot_amplitude2, mw_cnot_frequency2, mw_cnot_phase2, ssr_normalise,
                                  counts_per_readout, sync_gate_name='sync_gate', rf_channel='a_ch2',
                                  wait_ssr_length=1.0e-6, laser_ssr_length=500e-9):

        created_blocks = list()
        created_ensembles = list()

        # Add the laser initialization
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_laser_wait(name=laser_name, laser_length=laser_init_length, wait_length=wait_init_length)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        laser_wait_list = [laser_name, seq_param]

        # Add electron CNOT
        #ToDo; Add multiple frequency element
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_multiple_mw_pulse(name=mw_cnot_name, tau=mw_cnot_rabi_period/2.0,
                                          microwave_amplitude=mw_cnot_amplitude,
                                          microwave_frequency = mw_cnot_frequency,
                                          microwave_phase=mw_cnot_phase)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        mw_cnot_list = [mw_cnot_name, seq_param]

        # Add RF pulse
        created_blocks_tmp, created_ensembles_tmp, rf_list1, rf_list2 = \
            self._chopped_rf_pulse(name=rf_cnot_name, rf_duration=rf_cnot_duration, rf_amp=rf_cnot_amp,
                                   rf_freq=rf_cnot_freq, rf_phase=rf_cnot_phase, rf_channel=rf_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp

        # Add SSR 'sync' trigger
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_trigger(name=sync_gate_name, tau=1e-7, digital_channel=self.sync_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        sync_list = [sync_gate_name, seq_param]

        # Add SSR
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_singleshot_readout(name=ssr_name, mw_cnot_rabi_period=mw_cnot_rabi_period,
                                             mw_cnot_amplitude=mw_cnot_amplitude,
                                             mw_cnot_frequency=mw_cnot_frequency,
                                             mw_cnot_phase=mw_cnot_phase,
                                             mw_cnot_amplitude2=mw_cnot_amplitude2,
                                             mw_cnot_frequency2=mw_cnot_frequency2,
                                             mw_cnot_phase2=mw_cnot_phase2,
                                             ssr_normalise=ssr_normalise,
                                             wait_ssr_length=wait_ssr_length,
                                             laser_ssr_length=laser_ssr_length)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({'repetitions': counts_per_readout-1})
        ssr_list = [ssr_name, seq_param]

        return created_blocks, created_ensembles, laser_wait_list, mw_cnot_list, rf_list1, rf_list2, sync_list, ssr_list


    def _initialize_ssr(self, laser_name, laser_init_length, wait_init_length, initial_pi_pulse, ssr_name, mw_cnot_rabi_period,
                        mw_cnot_amplitude, mw_cnot_frequency, mw_cnot_phase, mw_cnot_amplitude2, mw_cnot_frequency2,
                        mw_cnot_phase2, ssr_normalise, counts_per_readout, sync_gate_name='sync_gate',
                        wait_ssr_length=1.0e-6, laser_ssr_length=500e-9):

        created_blocks = list()
        created_ensembles = list()
        if not initial_pi_pulse:
            # Add just laser initialization
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_laser_wait(name=laser_name, laser_length=laser_init_length, wait_length=wait_init_length)
        else:
            # Add an additional MW pi-pulse to initalize the NV into -1 or +1
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_laser_wait_pipulse(name=laser_name, laser_length=laser_init_length,  wait_length=wait_init_length)

        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        laser_wait_list = [laser_name, seq_param]

        # Add SSR 'sync' trigger
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_trigger(name=sync_gate_name, tau=1e-7, digital_channel=self.sync_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({})
        sync_list = [sync_gate_name, seq_param]

        # Add SSR
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_singleshot_readout(name=ssr_name, mw_cnot_rabi_period=mw_cnot_rabi_period,
                                             mw_cnot_amplitude=mw_cnot_amplitude,
                                             mw_cnot_frequency=mw_cnot_frequency,
                                             mw_cnot_phase=mw_cnot_phase,
                                             mw_cnot_amplitude2=mw_cnot_amplitude2,
                                             mw_cnot_frequency2=mw_cnot_frequency2,
                                             mw_cnot_phase2=mw_cnot_phase2,
                                             ssr_normalise=ssr_normalise,
                                             wait_ssr_length=wait_ssr_length,
                                             laser_ssr_length=laser_ssr_length)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        seq_param = self._customize_seq_para({'repetitions': counts_per_readout-1})
        ssr_list = [ssr_name, seq_param]

        return created_blocks, created_ensembles, laser_wait_list, sync_list, ssr_list


##########################################         Helper methods     ##########################################





    def generate_single_echo_s3(self, name='Hahn-Echo', rabi_period = 20e-9, tau=500e-9, microwave_amplitude=0.1,
                               microwave_frequency = 2.8e9):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        ### prevent granularity problems
        rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
        tau = self._adjust_to_samplingrate(tau, 4)

        # get pihalf element
        pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                              increment=0.0,
                                              amp=microwave_amplitude,
                                              freq=microwave_frequency,
                                              phase=0.0)
        # get pi elements
        pix_element = self._get_mw_element(length=rabi_period / 2,
                                           increment=0.0,
                                           amp=microwave_amplitude,
                                           freq=microwave_frequency,
                                           phase=0.0)

        tau_element = self._get_idle_element(length=tau - rabi_period / 2, increment= 0.0)

        # create XY8-N block element list
        block = PulseBlock(name=name)
        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        if (rabi_period+2*tau) * self.pulse_generator_settings['sample_rate'] < 4800:
            length_idle = 4800 / self.pulse_generator_settings['sample_rate'] - (rabi_period+2*tau)
            idle_element_extra = self._get_idle_element(length = length_idle, increment = 0.0)
            block.append(idle_element_extra)
        # actual XY8-N sequence
        block.append(pihalf_element)
        block.append(tau_element)
        block.append(pix_element)
        block.append(tau_element)
        block.append(pihalf_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences

    def generate_single_xy8_s3(self, name='XY8', rabi_period = 20e-9, tau=500e-9, microwave_amplitude=0.1,
                               microwave_frequency = 2.8e9, xy8N =1, ylast=False):

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        ### prevent granularity problems
        rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
        tau = self._adjust_to_samplingrate(tau, 4)

        # get pihalf element
        pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                              increment=0.0,
                                              amp=microwave_amplitude,
                                              freq=microwave_frequency,
                                              phase=0.0)
        if ylast:
            pihalf_y_element = self._get_mw_element(length=rabi_period / 4,
                                                    increment=0.0,
                                                    amp=microwave_amplitude,
                                                    freq=microwave_frequency,
                                                    phase=90.0)

        # get pi elements
        pix_element = self._get_mw_element(length=rabi_period / 2,
                                           increment=0.0,
                                           amp=microwave_amplitude,
                                           freq=microwave_frequency,
                                           phase=0.0)

        piy_element = self._get_mw_element(length=rabi_period / 2,
                                           increment=0.0,
                                           amp=microwave_amplitude,
                                           freq=microwave_frequency,
                                           phase=90.0)

        tauhalf_element = self._get_idle_element(length = tau/2-rabi_period/4, increment=0.0)
        tau_element = self._get_idle_element(length=tau - rabi_period / 2, increment= 0.0)

        # create XY8-N block element list
        block = PulseBlock(name=name)
        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        if (tau * 8 * xy8N) * self.pulse_generator_settings['sample_rate'] < 4800:
            length_idle = 4800 / self.pulse_generator_settings['sample_rate'] - (tau * 8 * xy8N)
            idle_element_extra = self._get_idle_element(length = length_idle, increment = 0.0)
            block.append(idle_element_extra)
        # actual XY8-N sequence
        block.append(pihalf_element)
        block.append(tauhalf_element)
        for n in range(xy8N):
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            if n != xy8N - 1:
                block.append(tau_element)
        block.append(tauhalf_element)
        if ylast:
            block.append(pihalf_y_element)
        else:
            block.append(pihalf_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_single_xy8_signal_s3(self, name='XY8', rabi_period = 20e-9, tau=500e-9, xy8N =1,
                                      signal_during_mw = False, lasty=False, signal_amplitude = 1,
                                      signal_frequency=1.0e6, signal_phase = 0.0):

        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        ### prevent granularity problems
        rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
        tau = self._adjust_to_samplingrate(tau, 4)

        if not signal_during_mw:

            # get pihalf element
            pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                                  increment=0.0,
                                                  amp=self.microwave_amplitude,
                                                  freq=self.microwave_frequency,
                                                  phase=0.0)
            if lasty:
                piyhalf_element = self._get_mw_element(length=rabi_period / 4,
                                                       increment=0.0,
                                                       amp=self.microwave_amplitude,
                                                       freq=self.microwave_frequency,
                                                       phase=90.0)

            # get pi elements
            pix_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=0.0)

            piy_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=90.0)

        else:
            # get pihalf element
            pihalf_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                           increment=0,
                                                           amps=[self.microwave_amplitude, signal_amplitude],
                                                           freqs=[self.microwave_frequency, signal_frequency],
                                                           phases=[0.0, signal_phase])
            if lasty:
                piyhalf_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                                increment=0,
                                                                amps=[self.microwave_amplitude, signal_amplitude],
                                                                freqs=[self.microwave_frequency, signal_frequency],
                                                                phases=[90.0, signal_phase])

            # get pi elements
            pix_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[0.0, signal_phase])

            piy_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[90.0, signal_phase])

        # get pure interaction elements
        tauhalf_element = self._get_mw_element(length=tau / 2.0 - rabi_period / 4,
                                               increment=0.0,
                                               amp=signal_amplitude,
                                               freq=signal_frequency,
                                               phase=signal_phase)

        tau_element = self._get_mw_element(length=tau - rabi_period / 2,
                                           increment=0.0,
                                           amp=signal_amplitude,
                                           freq=signal_frequency,
                                           phase=signal_phase)

        block = PulseBlock(name=name)
        if (tau + xy8N) * self.pulse_generator_settings['sample_rate'] < 4800:
            length_idle = 4800 / self.pulse_generator_settings['sample_rate'] - (tau + xy8N)
            idle_element_extra = self._get_idle_element(length_idle, 0.0)
            block.append(idle_element_extra)
        # actual XY8-N sequence
        block.append(pihalf_element)
        block.append(tauhalf_element)
        for n in range(xy8N):
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            if n != xy8N - 1:
                block.append(tau_element)
        block.append(tauhalf_element)
        if lasty:
            block.append(piyhalf_element)
        else:
            block.append(pihalf_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_single_xy4_signal_s3(self, name='XY4', rabi_period = 20e-9, tau=500e-9, xy4N =1,
                                      signal_during_mw = False, lasty=False, signal_amplitude = 1,
                                      signal_frequency=1.0e6, signal_phase = 0.0):

        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        ### prevent granularity problems
        rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
        tau = self._adjust_to_samplingrate(tau, 4)

        if not signal_during_mw:

            pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                                  increment=0.0,
                                                  amp=self.microwave_amplitude,
                                                  freq=self.microwave_frequency,
                                                  phase=0.0)

            if lasty:
                piyhalf_element = self._get_mw_element(length=rabi_period / 4,
                                                       increment=0.0,
                                                       amp=self.microwave_amplitude,
                                                       freq=self.microwave_frequency,
                                                       phase=90.0)

            # get pi elements
            pix_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=0.0)

            piy_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=90.0)

        else:
            pihalf_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                           increment=0,
                                                           amps=[self.microwave_amplitude, signal_amplitude],
                                                           freqs=[self.microwave_frequency, signal_frequency],
                                                           phases=[0.0, signal_phase])
            if lasty:
                piyhalf_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                                increment=0,
                                                                amps=[self.microwave_amplitude, signal_amplitude],
                                                                freqs=[self.microwave_frequency, signal_frequency],
                                                                phases=[90.0, signal_phase])

            pix_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[0.0, signal_phase])

            piy_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[90.0, signal_phase])

        # get pure interaction elements
        tauhalf_element = self._get_mw_element(length=tau / 2.0 - rabi_period / 4,
                                               increment=0.0,
                                               amp=signal_amplitude,
                                               freq=signal_frequency,
                                               phase=signal_phase)

        tau_element = self._get_mw_element(length=tau - rabi_period / 2,
                                           increment=0.0,
                                           amp=signal_amplitude,
                                           freq=signal_frequency,
                                           phase=signal_phase)

        block = PulseBlock(name=name)
        if (tau + xy4N) * self.pulse_generator_settings['sample_rate'] < 4800:
            length_idle = 4800 / self.pulse_generator_settings['sample_rate'] - (tau + xy4N)
            idle_element_extra = self._get_idle_element(length_idle, 0.0)
            block.append(idle_element_extra)
        # actual xy4-N sequence
        block.append(pihalf_element)
        block.append(tauhalf_element)
        for n in range(xy4N):
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            if n != xy4N - 1:
                block.append(tau_element)
            block.append(tauhalf_element)
        if lasty:
            block.append(piyhalf_element)
        else:
            block.append(pihalf_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_single_xy4_signal_adapted_readout_s3(self, name='XY4', rabi_period = 20e-9, tau=500e-9, xy4N =1,
                                                      signal_during_mw = False, lasty=False, signal_amplitude = 1,
                                                      signal_frequency=1.0e6, signal_phase = 0.0,
                                                      signal_amplitude_Hz = 1.0e3):

        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        ### prevent granularity problems
        rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
        tau = self._adjust_to_samplingrate(tau, 4)

        # compute the readout_axes:
        detuning = (signal_frequency - 1 / 2 / tau) * 2 * np.pi
        phase = 2 * signal_amplitude_Hz * (1 - np.cos(detuning * tau * 4 * xy4N)) / detuning / 2 / np.pi * 360

        if not signal_during_mw:

            pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                                  increment=0.0,
                                                  amp=self.microwave_amplitude,
                                                  freq=self.microwave_frequency,
                                                  phase=0.0)
            pihalf_readout_element = self._get_mw_element(length=rabi_period / 4,
                                                  increment=0.0,
                                                  amp=self.microwave_amplitude,
                                                  freq=self.microwave_frequency,
                                                  phase=phase)

            pix_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=0.0)

            piy_element = self._get_mw_element(length=rabi_period / 2,
                                               increment=0.0,
                                               amp=self.microwave_amplitude,
                                               freq=self.microwave_frequency,
                                               phase=90.0)

        else:
            pihalf_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                           increment=0,
                                                           amps=[self.microwave_amplitude, signal_amplitude],
                                                           freqs=[self.microwave_frequency, signal_frequency],
                                                           phases=[0.0, signal_phase])

            pihalf_readout_element = self._get_multiple_mw_element(length=rabi_period / 4,
                                                           increment=0,
                                                           amps=[self.microwave_amplitude, signal_amplitude],
                                                           freqs=[self.microwave_frequency, signal_frequency],
                                                           phases=[phase, signal_phase])

            pix_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[0.0, signal_phase])

            piy_element = self._get_multiple_mw_element(length=rabi_period / 2,
                                                        increment=0,
                                                        amps=[self.microwave_amplitude, signal_amplitude],
                                                        freqs=[self.microwave_frequency, signal_frequency],
                                                        phases=[90.0, signal_phase])

        # get pure interaction elements
        tauhalf_element = self._get_mw_element(length=tau / 2.0 - rabi_period / 4,
                                               increment=0.0,
                                               amp=signal_amplitude,
                                               freq=signal_frequency,
                                               phase=signal_phase)

        tau_element = self._get_mw_element(length=tau - rabi_period / 2,
                                           increment=0.0,
                                           amp=signal_amplitude,
                                           freq=signal_frequency,
                                           phase=signal_phase)

        block = PulseBlock(name=name)
        # additional time to fill up the waveform to 4800 samples if necessary
        if (tau+xy4N) * self.pulse_generator_settings['sample_rate'] < 4800:
            length_idle = 4800/self.pulse_generator_settings['sample_rate'] - (tau+xy4N)
            idle_element_extra = self._get_idle_element(length=length_idle, increment=0.0)
            block.append(idle_element_extra)
        # actual xy4-N sequence
        block.append(pihalf_element)
        block.append(tauhalf_element)
        for n in range(xy4N):
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            block.append(tau_element)
            block.append(pix_element)
            block.append(tau_element)
            block.append(piy_element)
            if n != xy4N - 1:
                block.append(tau_element)
        block.append(tauhalf_element)
        block.append(pihalf_readout_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences


######################################### Helper functions ########################################################



        ############################################# FIXME ##############################################################


def generate_single_T1_qdyne_s3(self, name='T1_qdyne', rabi_period=20e-9, tau=1e-6, seq_trig='d_ch1'):
    created_blocks = list()
    created_ensembles = list()
    created_sequences = list()
    ### prevent granularity problems
    rabi_period = self._adjust_to_samplingrate(rabi_period, 8)
    tau = self._adjust_to_samplingrate(tau, 4)

    # get pihalf element
    pihalf_element = self._get_mw_element(length=rabi_period / 4,
                                          increment=0.0,
                                          amp=self.microwave_amplitude,
                                          freq=self.microwave_frequency,
                                          phase=0.0)
    # get tau element
    tau_element = self._get_trigger_element(length=tau,
                                            increment=0.0,
                                            channels=seq_trig)
    # create single_T1_qdyne block element list
    block = PulseBlock(name=name)
    # block.append(waiting_element)
    block.append(tau_element)
    block.append(pihalf_element)

    created_blocks.append(block)

    # Create block ensemble
    block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
    block_ensemble.append((block.name, 0))

    # add metadata to invoke settings
    block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks,
                                                    alternating=False,
                                                    controlled_variable=np.array([0]))
    # append ensemble to created ensembles
    created_ensembles.append(block_ensemble)
    return created_blocks, created_ensembles, created_sequences


def generate_ssr_T1_Qdyne(self, name='SSR-T1qdyne', tau=1.0e-9, seq_trig='d_ch1',
                          laser_name='laser_wait', laser_length=1e-6, wait_length=1e-6,
                          rf_cnot_name='RF', rf_cnot_freq=1.0e6, rf_cnot_amp=0.1, rf_cnot_duration=100e-6, rf_cnot_phase=0,
                          ssr_name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=1.0, mw_cnot_frequency=2.8e9,
                          mw_cnot_phase=0, mw_cnot_amplitude2=1.0, mw_cnot_frequency2=2.8e9,
                          mw_cnot_phase2=0, ssr_normalise=True, counts_per_readout=1000):
    created_blocks = list()
    created_ensembles = list()
    created_sequences = list()
    para_dict = locals()

    # generate the Rabi pieces

    para_list = list()
    created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
        self.generate_single_T1_qdyne_s3(name=name, tau=tau, microwave_amplitude=self.microwave_amplitude,
                                         microwave_frequency=self.microwave_frequency, microwave_phase=0.0,
                                         seq_trig=seq_trig)
    created_blocks += created_blocks_tmp
    created_ensembles += created_ensembles_tmp
    seq_param = self._customize_seq_para({})
    para_list.append([name, seq_param])

    created_blocks, created_ensembles, sequence = \
        self._standard_ssr(created_blocks, created_ensembles, para_list, para_dict)

    self._add_metadata_to_settings(sequence, created_blocks=list(), alternating=False, laser_ignore_list=list(),
                                   controlled_variable=[0], units=('s', ''), number_of_lasers=1,
                                   counting_length=laser_length * self.normalised_safety if ssr_normalise else laser_length * 1.4)
    created_sequences.append(sequence)
    return created_blocks, created_ensembles, created_sequences












