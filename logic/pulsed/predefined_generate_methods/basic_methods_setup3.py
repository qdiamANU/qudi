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


import numpy as np
from logic.pulsed.pulse_objects import PulseBlock, PulseBlockEnsemble
from logic.pulsed.pulse_objects import PredefinedGeneratorBase


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


class BasicPredefinedGeneratorS3(PredefinedGeneratorBase):
    """

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def generate_ajh_test_function(self, name='ajh_test_function', mw_freq=1.0e6, tau_start=2.0e-6, tau_step=0.1e-6,
                                   num_of_points=5):
        """
        :param name:
        :param freq_list:
        :return:
        """

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get readout element
        readout_element = self._get_readout_element()

        # get tau array for measurement
        tau_array = tau_start + np.arange(num_of_points) * tau_step

        # create the mw element
        mw_element = self._get_mw_element(length=tau_start,
                                          increment=tau_step,
                                          amp=self.microwave_amplitude,
                                          freq=mw_freq,
                                          phase=0)

        # create second mw element
        mw_element2 = self._get_mw_element(length=tau_start*2,
                                          increment=tau_step*2,
                                          amp=self.microwave_amplitude,
                                          freq=mw_freq,
                                          phase=0)

        waiting_element = self._get_idle_element(length=self.wait_time,
                                                 increment=0)

        # Create block and append to created_blocks list
        test_block = PulseBlock(name=name)
        test_block.append(mw_element)
        test_block.append(waiting_element)
        test_block.append(mw_element2)
        test_block.extend(readout_element)
        created_blocks.append(test_block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)

        # Create and append sync trigger block if needed
        if self.sync_channel:
            sync_block = PulseBlock(name='sync_trigger')
            sync_block.append(self._get_sync_element())
            created_blocks.append(sync_block)
            block_ensemble.append((sync_block.name, 0))

        # append mw/readout block
        block_ensemble.append((test_block.name, num_of_points - 1))

        # add metadata to invoke settings later on
        block_ensemble.measurement_information['alternating'] = False
        block_ensemble.measurement_information['laser_ignore_list'] = list()
        block_ensemble.measurement_information['controlled_variable'] = tau_array
        block_ensemble.measurement_information['units'] = ('s', '')
        block_ensemble.measurement_information['labels'] = ('Tau', 'Signal')
        block_ensemble.measurement_information['number_of_lasers'] = num_of_points
        block_ensemble.measurement_information['counting_length'] = self._get_ensemble_count_length(
            ensemble=block_ensemble, created_blocks=created_blocks)

        # Append ensemble to created_ensembles list
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences

    # def generate_ajh_test_sequence(self, name='ajh_test_sequence', counts_per_readout = 5, mw_freq=1.0e6, tau_start=2.0e-6, tau_step=0.1e-6,
    #                                num_of_points=5):
    #     """
    #     :param name:
    #     :param freq_list:
    #     :return:
    #     """
    #
    #     created_blocks = list()
    #     created_ensembles = list()
    #     created_sequences = list()
    #
    #     laser_name = 'laser_initialisation'
    #     laser_length = 3.0e-6
    #     wait_length = self.wait_time
    #
    #     # Add the laser initialization
    #     created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
    #         self.generate_laser_wait(name=laser_name, laser_length=laser_length, wait_length=wait_length)
    #     created_blocks += created_blocks_tmp
    #     created_ensembles += created_ensembles_tmp
    #     seq_param = self._customize_seq_para({})
    #     laser_wait_list = [laser_name, seq_param]
    #
    #     # Add MW pulse
    #     mw_cnot_name = 'cNnotE_gate'
    #     mw_cnot_duration = 3.0e-6
    #     mw_cnot_amp = 2.0
    #     mw_cnot_freq = 1.0e6
    #     mw_cnot_phase = 0
    #     created_blocks_tmp, created_ensembles_tmp, rf_list1, rf_list2 = \
    #         self._chopped_rf_pulse(name=mw_cnot_name, rf_duration=mw_cnot_duration, rf_amp=mw_cnot_amp,
    #                                rf_freq=mw_cnot_freq, rf_phase=mw_cnot_phase)
    #     created_blocks += created_blocks_tmp
    #     created_ensembles += created_ensembles_tmp
    #     seq_param = self._customize_seq_para({})
    #     mw_cnot_list = [mw_cnot_name, seq_param]
    #
    #     # add conventional readout
    #     readout_name='readout_conventional'
    #     readout_element = self._get_readout_element()
    #     readout_block = PulseBlock(name=readout_name)
    #     readout_block.extend(readout_element)
    #     created_blocks.append(readout_block)
    #     readout_block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
    #     readout_block_ensemble.append((readout_block.name, 0))
    #     seq_param = self._customize_seq_para({'repetitions': counts_per_readout - 1})
    #     readout_list = [readout_name, seq_param]
    #
    #
    #     # bring the individual blocks in the correct order
    #     element_list = list()
    #     # para_list =
    #     # for ii in range(len(para_list)):
    #     #     element_list.append(laser_wait_list.copy())
    #     #     element_list.append(para_list[ii])
    #     #     element_list.append(mw_cnot_list.copy())
    #     #     element_list.append(readout_list.copy())
    #
    #     element_list.append(laser_wait_list.copy())
    #     element_list.append(mw_cnot_list.copy())
    #     element_list.append(readout_list.copy())
    #     # make sequence continous
    #     element_list = self._make_sequence_continous(element_list)
    #
    #     print('_standard_ssr: laser_wait_list = {}'.format(laser_wait_list))
    #     print('_standard_ssr: mw_cnot_list = {}'.format(mw_cnot_list))
    #     print('_standard_ssr: readout_list = {}'.format(readout_list))
    #     print('_standard_ssr: element_list = {}'.format(element_list))
    #
    #     sequence_name = 'test_sequence'
    #     sequence = PulseSequence(name=sequence_name, ensemble_list=element_list, rotating_frame=False)
    #
    #     return created_blocks, created_ensembles, sequence

    def generate_ramsey_from_list_s3(self, name='ramsey', tau_list='[1e-6, 2e-6]', alternating = True):
        """

        """

        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        tau_array = [n.strip() for n in tau_list]

        # get readout element
        readout_element = self._get_readout_element()

        # get pihalf element
        pihalf_element = self._get_mw_element(length=self.rabi_period / 4,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0.0)

        if alternating:
            if self.microwave_channel.startswith('a'):
                pi3half_element = self._get_mw_element(length=self.rabi_period / 4,
                                                       increment=0,
                                                       amp=self.microwave_amplitude,
                                                       freq=self.microwave_frequency,
                                                       phase=180)
            else:
                pi3half_element = self._get_mw_element(length=3 * self.rabi_period / 4,
                                                       increment=0,
                                                       amp=self.microwave_amplitude,
                                                       freq=self.microwave_frequency,
                                                       phase=0)

        # Create block and append to created_blocks list
        block = PulseBlock(name=name)
        for tau in tau_array:
            block.append(pihalf_element)
            tau_element = self._get_idle_element(length=tau, increment=0)
            block.append(tau_element)
            block.append(tau_element)
            block.append(pihalf_element)
            block.extend(readout_element)

            if alternating:
                block.append(pihalf_element)
                block.append(tau_element)
                block.append(pi3half_element)
                block.extend(readout_element)

        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences


    def generate_t1_s3(self, name='T1', tau_start=1.0e-6, tau_step=1.0e-6,
                    num_of_points=50, alternating = False):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        tau_array = tau_start + np.arange(num_of_points) * tau_step

        # get readout element
        readout_element = self._get_readout_element()

        if alternating: # get pi element
            pi_element = self._get_mw_element(length=self.rabi_period / 2,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0)

        tau_element = self._get_idle_element(length=tau_start, increment=tau_step)
        block = PulseBlock(name=name)
        block.append(tau_element)
        block.extend(readout_element)
        if alternating:
            block.append(pi_element)
            block.append(tau_element)
            block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, num_of_points - 1))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences


    def generate_t1_exponential_s3(self, name='T1_exp', tau_start=1.0e-6, tau_end=1.0e-6,
                    num_of_points=50, alternating=False):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        if tau_start == 0.0:
            tau_array = np.geomspace(1e-9, tau_end, num_of_points - 1)
            tau_array = np.insert(tau_array, 0, 0.0)
        else:
            tau_array = np.geomspace(tau_start, tau_end, num_of_points)

        # get readout element
        readout_element = self._get_readout_element()

        if alternating:  # get pi element
            pi_element = self._get_mw_element(length=self.rabi_period / 2,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0)
        block = PulseBlock(name=name)
        for tau in tau_array:
            tau_element = self._get_idle_element(length=tau, increment=0.0)
            block.append(tau_element)
            block.extend(readout_element)
            if alternating:
                block.append(pi_element)
                block.append(tau_element)
                block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences


    def generate_hahnecho_exp_s3(self, name='hahn_echo', tau_start=1.0e-6, tau_end=10.0e-6,
                          num_of_points=50, alternating=True):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        if tau_start == 0.0:
            tau_array = np.geomspace(1e-9, tau_end, num_of_points - 1)
            tau_array = np.insert(tau_array, 0, 0.0)
        else:
            tau_array = np.geomspace(tau_start, tau_end, num_of_points)

        # get readout element
        readout_element = self._get_readout_element()

        pihalf_element = self._get_mw_element(length=self.rabi_period / 4,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0)
        pi_element = self._get_mw_element(length=self.rabi_period / 2,
                                          increment=0,
                                          amp=self.microwave_amplitude,
                                          freq=self.microwave_frequency,
                                          phase=0)
        # Use a 180 deg phase shiftet pulse as 3pihalf pulse if microwave channel is analog
        if self.microwave_channel.startswith('a'):
            pi3half_element = self._get_mw_element(length=self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=180)
        else:
            pi3half_element = self._get_mw_element(length=3 * self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=0)


        # Create block and append to created_blocks list
        block = PulseBlock(name=name)
        for tau in tau_array:
            tau_element = self._get_idle_element(length=tau, increment=0.0)
            block.append(pihalf_element)
            block.append(tau_element)
            block.append(pi_element)
            block.append(tau_element)
            block.append(pihalf_element)
            block.extend(readout_element)
            if alternating:
                block.append(pihalf_element)
                block.append(tau_element)
                block.append(pi_element)
                block.append(tau_element)
                block.append(pi3half_element)
                block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, 0))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_hahnecho_s3(self, name='hahn_echo', tau_start=1.0e-6, tau_step=1.0e-6,
                          num_of_points=50, alternating=True):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        tau_array = tau_start + np.arange(num_of_points) * tau_step

        # get readout element
        readout_element = self._get_readout_element()

        pihalf_element = self._get_mw_element(length=self.rabi_period / 4,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0)
        pi_element = self._get_mw_element(length=self.rabi_period / 2,
                                          increment=0,
                                          amp=self.microwave_amplitude,
                                          freq=self.microwave_frequency,
                                          phase=0)
        # Use a 180 deg phase shiftet pulse as 3pihalf pulse if microwave channel is analog
        if self.microwave_channel.startswith('a'):
            pi3half_element = self._get_mw_element(length=self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=180)
        else:
            pi3half_element = self._get_mw_element(length=3 * self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=0)
        tau_element = self._get_idle_element(length=tau_start, increment=tau_step)

        # Create block and append to created_blocks list
        block = PulseBlock(name=name)
        block.append(pihalf_element)
        block.append(tau_element)
        block.append(pi_element)
        block.append(tau_element)
        block.append(pihalf_element)
        block.extend(readout_element)
        if alternating:
            block.append(pihalf_element)
            block.append(tau_element)
            block.append(pi_element)
            block.append(tau_element)
            block.append(pi3half_element)
            block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, num_of_points - 1))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_xy8_tau_s3(self, name='xy8_tau', tau_start=0.5e-6, tau_step=0.01e-6, num_of_points=50,
                         xy8_order=4, alternating=True):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get tau array for measurement ticks
        tau_array = tau_start + np.arange(num_of_points) * tau_step
        # calculate "real" start length of tau due to finite pi-pulse length
        real_start_tau = max(0, tau_start - self.rabi_period / 2)

        # get readout element
        readout_element = self._get_readout_element()
        pihalf_element = self._get_mw_element(length=self.rabi_period / 4,
                                              increment=0,
                                              amp=self.microwave_amplitude,
                                              freq=self.microwave_frequency,
                                              phase=0)
        # Use a 180 deg phase shiftet pulse as 3pihalf pulse if microwave channel is analog
        if self.microwave_channel.startswith('a'):
            pi3half_element = self._get_mw_element(length=self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=180)
        else:
            pi3half_element = self._get_mw_element(length=3 * self.rabi_period / 4,
                                                   increment=0,
                                                   amp=self.microwave_amplitude,
                                                   freq=self.microwave_frequency,
                                                   phase=0)
        pix_element = self._get_mw_element(length=self.rabi_period / 2,
                                           increment=0,
                                           amp=self.microwave_amplitude,
                                           freq=self.microwave_frequency,
                                           phase=0)
        piy_element = self._get_mw_element(length=self.rabi_period / 2,
                                           increment=0,
                                           amp=self.microwave_amplitude,
                                           freq=self.microwave_frequency,
                                           phase=90)
        tauhalf_element = self._get_idle_element(length=real_start_tau / 2, increment=tau_step / 2)
        tau_element = self._get_idle_element(length=real_start_tau, increment=tau_step)

        # Create block and append to created_blocks list
        block = PulseBlock(name=name)
        block.append(pihalf_element)
        block.append(tauhalf_element)
        for n in range(xy8_order):
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
            if n != xy8_order - 1:
                block.append(tau_element)
        block.append(tauhalf_element)
        block.append(pihalf_element)
        block.extend(readout_element)
        if alternating:
            block.append(pihalf_element)
            block.append(tauhalf_element)
            for n in range(xy8_order):
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
                if n != xy8_order - 1:
                    block.append(tau_element)
            block.append(tauhalf_element)
            block.append(pi3half_element)
            block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=True)
        block_ensemble.append((block.name, num_of_points - 1))

        # Create and append sync trigger block if needed
        created_blocks, block_ensemble = self._add_trigger(created_blocks, block_ensemble)
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=alternating,
                                                        controlled_variable=tau_array)

        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences




    def generate_rabi_DTG(self, name='Rabi_dtg', tau_start=1e-9, tau_step=1e-9, num_of_points=50,
                          seq_trig='d_ch1', shorten_record=300e-9):
        """
        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        tau_array = tau_start + np.arange(num_of_points) * tau_step

        # get MW element (here just DC trigger)
        mw_element = self._get_trigger_element(length=tau_start, increment=tau_step,  channels=seq_trig)
        # get readout element
        readout_element = self._get_readout_element()
        block = PulseBlock(name=name)
        # Create element list for Rabi PulseBlock
        block.append(mw_element)
        block.extend(readout_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, num_of_points-1))

        # add metadata to invoke settings
        sequence_length = self._get_ensemble_count_length(ensemble=block_ensemble, created_blocks=created_blocks)
        block_ensemble = \
            self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False,
                                                        controlled_variable=tau_array,
                                                        counting_length=sequence_length-shorten_record)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences


    def generate_multifreq_rabi(self, name='multifreq_rabi', amps='[1, 1, 1]', freqs='[10e6, 11e6, 12e6]', phases='[0, 0, 0]',
                                tau_start=10.0e-9, tau_step=10.0e-9, num_of_points=50):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # amps_list = list(map(float, amps.split(',')))
        # freqs_list = list(map(float, freqs.split(',')))
        # phases_list = list(map(float, phases.split(',')))
        amps_list = eval(amps)
        freqs_list = eval(freqs)
        phases_list = eval(phases)


        # get tau array for measurement ticks
        tau_array = tau_start + np.arange(num_of_points) * tau_step

        # create the laser_mw element
        mw_element = self._get_multiple_mw_element(length=tau_start,
                                          increment=tau_step,
                                          amps=amps_list,
                                          freqs=freqs_list,
                                          phases=phases_list)
        readout_element = self._get_readout_element(wait_time=self.wait_time,
                                                    length=self.laser_length,
                                                    trigger=True)
        # waiting_element = self._get_idle_element(length=self.wait_time,
        #                                          increment=0)
        # laser_element = self._get_laser_gate_element(length=self.laser_length,
        #                                              increment=0)
        # delay_element = self._get_delay_gate_element()

        # Create block and append to created_blocks list
        rabi_block = PulseBlock(name=name)
        rabi_block.append(mw_element)
        # rabi_block.append(laser_element)
        # rabi_block.append(delay_element)
        # rabi_block.append(waiting_element)
        rabi_block.extend(readout_element)
        created_blocks.append(rabi_block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((rabi_block.name, num_of_points - 1))

        # Create and append sync trigger block if needed
        if self.sync_channel:
            sync_block = PulseBlock(name='sync_trigger')
            sync_block.append(self._get_sync_element())
            created_blocks.append(sync_block)
            block_ensemble.append((sync_block.name, 0))

        # add metadata to invoke settings later on
        block_ensemble.measurement_information['alternating'] = False
        block_ensemble.measurement_information['laser_ignore_list'] = list()
        block_ensemble.measurement_information['controlled_variable'] = tau_array
        block_ensemble.measurement_information['units'] = ('s', '')
        block_ensemble.measurement_information['labels'] = ('Tau', 'Signal')
        block_ensemble.measurement_information['number_of_lasers'] = num_of_points
        block_ensemble.measurement_information['counting_length'] = self._get_ensemble_count_length(
            ensemble=block_ensemble, created_blocks=created_blocks)

        # Append ensemble to created_ensembles list
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences
















