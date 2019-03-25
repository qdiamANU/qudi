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
from logic.pulsed.pulse_objects import PulseBlock, PulseBlockEnsemble, PulseSequence
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

class HelperMethods(PredefinedGeneratorBase):
    """

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # todo: the minimum number of samples for a sequence waveform is currently hard-coded into this file.
        # Defining a single parameter here is better than just using a number throughout the class, as before
        # But the class should really be re-written so that all of this is dealt with in the hardware file
        # the minimum waveform length for the Spectrum AWG is 384 samples. Any required idle time is already
        # applied in the AWG hardware file, when writing the waveform for each pulse ensemble
        self.minimum_number_of_samples = 384

##########################################         Laser Helper methods     ##########################################


    def generate_laser(self, name='laser', tau=1e-6):
        """ Generates an element where the orange laser is turned on.

        @param str name: Name of the PulseBlockEnsemble
        @param float tau: duration in seconds


        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()

        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau, 2)
        ### fill up to minimum length of self.minimum_number_of_samples if necessary
        if tau * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            tau = self._adjust_to_samplingrate(self.minimum_number_of_samples/self.pulse_generator_settings['sample_rate'], 2)
            self.log.warning('Laser duration adusted to minimum sequence length of {0} for samplingrate '
                             '{1}'.format(tau, self.pulse_generator_settings['sample_rate']))

        # get the laser element
        laser_element = self._get_trigger_element(length=tau, increment=0.0, channels=self.laser_channel)

        block = PulseBlock(name=name)
        block.append(laser_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks,
                                                        alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, list()

    def generate_laser_wait(self, name='laser_wait', laser_length=500e-9, wait_length = 1e-6, trigger=True):
        """ Generates Laser pulse and waiting (idle) time.

        @param str name: Name of the PulseBlockEnsemble
        @param float length: laser duration in seconds
        @param float amp: In case of analogue laser channel this value will be the laser on voltage.

        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()

        ### prevent granularity problems
        laser_length = self._adjust_to_samplingrate(laser_length, 2)
        ### fill up to minimum length of self.minimum_number_of_samples if necessary
        if (laser_length + wait_length) * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            wait_length = self._adjust_to_samplingrate(
                self.minimum_number_of_samples / self.pulse_generator_settings['sample_rate'] - laser_length, 2)
        else:
            wait_length = self._adjust_to_samplingrate(wait_length, 2)

        # create the laser element
        if trigger:
            laser_element = self._get_laser_gate_element(length=laser_length, increment=0)
        else:
            laser_element = self._get_laser_element(length=laser_length, increment=0)
        waiting_element = self._get_idle_element(length=wait_length, increment=0.0)
        # Create the element list
        block = PulseBlock(name=name)
        block.append(laser_element)
        block.append(waiting_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, list()


    def generate_laser_wait_pipulse(self, name='laser_wait_pi', laser_length=500e-9, wait_length = 1e-6):
        """ Generates Laser pulse and waiting (idle) time.

        @param str name: Name of the PulseBlockEnsemble
        @param float length: laser duration in seconds
        @param float amp: In case of analogue laser channel this value will be the laser on voltage.

        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        # create the laser element
        laser_element = self._get_laser_gate_element(length=laser_length, increment=0)
        waiting_element = self._get_idle_element(length=wait_length, increment=0.0)
        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(self.rabi_period/2.0, 4)
        mw_element = self._get_mw_element(length=tau,
                                          increment=0.0,
                                          amp=self.microwave_amplitude,
                                          freq=self.microwave_frequency,
                                          phase=0.0)
        # Create the element list
        block = PulseBlock(name=name)
        block.extend([laser_element, waiting_element, mw_element])
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences



    def generate_pipulse_trigger_laser_wait(self, name='laser_wait_pi', tau = 1e-6, laser_length=500e-9,
                                             wait_length = 1e-6):
        """ Generates Laser pulse and waiting (idle) time.

        @param str name: Name of the PulseBlockEnsemble
        @param float length: laser duration in seconds
        @param float amp: In case of analogue laser channel this value will be the laser on voltage.

        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        # get readout element
        laser_element = self._get_laser_gate_element(length=laser_length, increment=0)
        waiting_element = self._get_idle_element(length=wait_length, increment=0.0)
        seqtrig_element = self._get_sync_element()

        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau, 4)
        mw_element = self._get_mw_element(length=tau,
                                          increment=0.0,
                                          amp=self.microwave_amplitude,
                                          freq=self.microwave_frequency,
                                          phase=0.0)

        block = PulseBlock(name=name)
        block.extend([mw_element, seqtrig_element, laser_element, waiting_element])
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences

################################### idle helper methods ######################################


    def generate_idle_s3(self, name='idle', tau=1e-6):
        """ Generates waiting (idle) time.

        @param str name: Name of the PulseBlockEnsemble
        @param float tau: duration in seconds


        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()

        ### fill up to minimum length of self.minimum_number_of_samples if necessary
        if tau * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            needed_extra_time = self.minimum_number_of_samples / self.pulse_generator_settings['sample_rate'] - tau
        else:
            needed_extra_time = 0

        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau + needed_extra_time, 4)

        # get the idle element
        idle_element = self._get_idle_element(length=tau, increment=0.0)
        # Create the element list
        block = PulseBlock(name=name)
        block.append(idle_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, list()

################################### Trigger helper methods ######################################

    def generate_trigger(self, name='trigger', tau=1e-6, digital_channel='d_ch1'):
        """ Generates an element with a digital trigger

        @param str name: Name of the PulseBlockEnsemble
        @param float tau: duration in seconds


        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()

        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau, 2)
        ### fill up to minimum length of self.minimum_number_of_samples if necessary
        if tau * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            tau = self._adjust_to_samplingrate(self.minimum_number_of_samples/self.pulse_generator_settings['sample_rate'], 2)
            self.log.warning('Trigger duration adjusted to minimum sequence length of {0} for samplingrate '
                             '{1}'.format(tau, self.pulse_generator_settings['sample_rate']))

        # get the trigger element
        trigger_element = self._get_trigger_element(length=tau, increment=0.0, channels=digital_channel)

        block = PulseBlock(name=name)
        block.append(trigger_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks,
                                                        alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, list()



    def generate_trigger_wait(self, name='trigger_wait', trigger_length=1e-6, wait_length=1e-6, digital_channel = 'd_ch1'):
        """ Generates an element where the orange laser is turned on.

        @param str name: Name of the PulseBlockEnsemble
        @param float tau: duration in seconds


        @return object: the generated PulseBlockEnsemble object.
        """
        created_blocks = list()
        created_ensembles = list()

        ### prevent granularity problems
        trigger_length = self._adjust_to_samplingrate(trigger_length, 2)
        ### fill up to minimum length of self.minimum_number_of_samples if necessary
        if (trigger_length+wait_length) * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            idle_length = self._adjust_to_samplingrate(self.minimum_number_of_samples/self.pulse_generator_settings['sample_rate']-trigger_length, 2)
        else:
            idle_length = self._adjust_to_samplingrate(wait_length, 2)

        # get the trigger and idle elements
        trigger_element = self._get_trigger_element(length=trigger_length, increment=0.0, channels=digital_channel)
        idle_element = self._get_idle_element(length=idle_length, increment=0.0)

        block = PulseBlock(name=name)
        block.append(trigger_element)
        block.append(idle_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks,
                                                        alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, list()

    ################################## sequence parts ######################################

    def generate_singleshot_readout_s3(self, name='SSR', mw_cnot_rabi_period=20e-9, mw_cnot_amplitude=0.1,
                                    mw_cnot_frequency=2.8e9, mw_cnot_phase=0, mw_cnot_amplitude2=0.1,
                                    mw_cnot_frequency2=2.8e9, mw_cnot_phase2=0, ssr_normalise=True):
        """

        """
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()

        ### prevent granularity problems
        mw_cnot_rabi_period = self._adjust_to_samplingrate(mw_cnot_rabi_period, 4)

        # get mw pi pulse block
        mw_pi_element = self._get_multiple_mw_element(length=mw_cnot_rabi_period / 2,
                                                      increment=0.0,
                                                      amps=mw_cnot_amplitude,
                                                      freqs=mw_cnot_frequency,
                                                      phases=mw_cnot_phase)

        trigger_element = self._get_sync_element()

        readout_element = self._get_readout_element()
        block = PulseBlock(name=name)
        block.append(mw_pi_element)
        block.append(trigger_element)
        block.extend(readout_element)

        if ssr_normalise:
            time_between_trigger = self.laser_length + self.wait_time + self.laser_delay
            if time_between_trigger > self.laser_length * 1.4:
                wait = time_between_trigger - self.laser_length * 1.4
                extra_waiting_element = self._get_idle_element(length=wait * 1.2, increment=0)
            mw_pi_element2 = self._get_multiple_mw_element(length=mw_cnot_rabi_period / 2,
                                                           increment=0.0,
                                                           amps=mw_cnot_amplitude2,
                                                           freqs=mw_cnot_frequency2,
                                                           phases=mw_cnot_phase2)
            waiting_element = self._get_idle_element(length=self.laser_length + 200e-9, increment=0)

            if self.laser_length + self.wait_time + self.laser_delay > self.laser_length * 1.4:
                block.append(extra_waiting_element)
                #
            block.append(mw_pi_element2)
            block.append(trigger_element)
            block.append(waiting_element)
            block.extend(readout_element)
        created_blocks.append(block)
        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=list(),
                                                        controlled_variable=[0],
                                                        counting_length=self.laser_length * 1.8 if ssr_normalise
                                                                  else self.laser_length * 1.4)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)

        return created_blocks, created_ensembles, created_sequences





################################## MW + RF pulse methods ######################################


    def generate_single_mw_pulse(self, name='MW_pulse', tau=1e-6, microwave_amplitude=1.0,
                                    microwave_frequency = 1e6, microwave_phase=0.0, channel='none'):

        # In sequence mode there is a minimum waveform length of self.minimum_number_of_samples sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau, 4)

        mw_element = self._get_mw_element(length=tau,
                                          increment=0.0,
                                          amp=microwave_amplitude,
                                          freq=microwave_frequency,
                                          phase=microwave_phase,
                                          channel=channel)

        # Create PulseBlock object
        block = PulseBlock(name=name)
        if tau * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            length_idle = self.minimum_number_of_samples/self.pulse_generator_settings['sample_rate'] -tau
            idle_element = self._get_idle_element(length = length_idle, increment= 0.0)
            block.append(idle_element)

        block.append(mw_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences

    def generate_multiple_mw_pulse(self, name='MW_pulse', tau=1e-6, microwave_amplitude=1.0,
                                   microwave_frequency=1e6, microwave_phase=0.0):

        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        created_blocks = list()
        created_ensembles = list()
        created_sequences = list()
        ### prevent granularity problems
        tau = self._adjust_to_samplingrate(tau, 4)

        mw_element = self._get_multiple_mw_element(length=tau,
                                                   increment=0.0,
                                                   amps=microwave_amplitude,
                                                   freqs=microwave_frequency,
                                                   phases=microwave_phase)

        # Create PulseBlock object
        block = PulseBlock(name=name)
        if tau * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            length_idle = self.minimum_number_of_samples / self.pulse_generator_settings['sample_rate'] - tau
            idle_element = self._get_idle_element(length=length_idle, increment=0.0)
            block.append(idle_element)

        block.append(mw_element)
        created_blocks.append(block)

        # Create block ensemble
        block_ensemble = PulseBlockEnsemble(name=name, rotating_frame=False)
        block_ensemble.append((block.name, 0))
        # add metadata to invoke settings
        block_ensemble = self._add_metadata_to_settings(block_ensemble, created_blocks=created_blocks, alternating=False)
        # append ensemble to created ensembles
        created_ensembles.append(block_ensemble)
        return created_blocks, created_ensembles, created_sequences

    def generate_single_xy8_s3(self, name='XY8', rabi_period=20e-9, tau=500e-9, microwave_amplitude=0.1,
                               microwave_frequency=2.8e9, xy8N=1, ylast=False):

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

        tauhalf_element = self._get_idle_element(length=tau / 2 - rabi_period / 4, increment=0.0)
        tau_element = self._get_idle_element(length=tau - rabi_period / 2, increment=0.0)

        # create XY8-N block element list
        block = PulseBlock(name=name)
        # In sequence mode there is a minimum waveform length of 4800 sample. If the pulse is to short add an
        # extra idle time before the pulse to take that into account
        if (tau * 8 * xy8N) * self.pulse_generator_settings['sample_rate'] < self.minimum_number_of_samples:
            length_idle = 4800 / self.pulse_generator_settings['sample_rate'] - (tau * 8 * xy8N)
            idle_element_extra = self._get_idle_element(length=length_idle, increment=0.0)
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


################################################## RF methods ###################################


    def _chopped_rf_pulse(self, name, rf_duration, rf_amp, rf_freq, rf_phase, rf_channel):
        """
        The function of chopped_rf_pulse is to reduce the data upload burden to the AWG, by 'chopping' the rf
        pulse into a number of periods

        generates two waveforms - rf pulse 1 which is an integer number of rf periods (or idle if pulse is shorter than
        one period), and rf pulse 2, which is the remaining fractional period.
        The output is then two sequence elements, where rf pulse 1 is repeated N times, and rf pulse 2 is played once.
        The number of repetitions, N, is chosen to match the desired total rf pulse duration

        Unfortunately, this technique doesn't work well with Spectrum AWGs (at least how the ANU one is set up as of
        25/2/2019), as there is a ca. 10ns idle period between repetitions of a waveform

        :param name:
        :param rf_duration:
        :param rf_amp:
        :param rf_freq:
        :param rf_phase:
        :param rf_channel:
        :return: created_blocks, created_ensembles, list1, list2
        list 1 = [name+'1', seq_param] - details for rf pulse 1
        list 1 = [name+'2', seq_param2] - details for rf pulse 2
        """


        rf_period = 1.0/rf_freq
        n_periods = np.floor(rf_duration / rf_period)
        # print('\n n_periods = {}'.format(n_periods))

        if n_periods == 0:
            # If there is not more than 1 period just make rf pulse 1 an idle period with minimal length
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_idle_s3(name=name + '1',
                                      tau=self.minimum_number_of_samples / self.pulse_generator_settings['sample_rate'])
            n_repetitions = 0
            rf_pulse2_duration = rf_duration
        else:
            # min. number of rf periods limited by min. number of samples for a waveform in an AWG sequence
            rf_pulse1_n_periods = np.ceil((self.minimum_number_of_samples / self.sample_rate) / rf_period)
            # rf_pulse1_duration multiplicity has to be two to avoid granularity problems with AWG
            rf_pulse1_duration = self._adjust_to_samplingrate(rf_pulse1_n_periods * rf_period, 2)
            # print('rf_pulse1_n_periods = {}, rf_pulse1_duration = {}'.format(rf_pulse1_n_periods, rf_pulse1_duration))

            # calculate duration for rf pulse 2
            rf_pulse2_duration = rf_duration % rf_pulse1_duration
            rf_pulse2_duration = self._adjust_to_samplingrate(rf_pulse2_duration + rf_period, 2)
            # extra ps length helps to prevent granularity problems (from Ulm - untested with ANU hardware):
            rf_pulse2_duration += 1e-12
            # print('rf_pulse2_duration = {}'.format(rf_pulse2_duration))

            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name+'1', tau=rf_pulse1_duration,
                                              microwave_amplitude=rf_amp, microwave_frequency=rf_freq,
                                              microwave_phase=rf_phase, channel=rf_channel)
            n_repetitions = int(np.floor(rf_duration / rf_pulse1_duration)-1)
            # print('n_repetitions = {}'.format(n_repetitions))

        created_blocks = created_blocks_tmp
        created_ensembles = created_ensembles_tmp
        # add sequence parameters
        seq_param = self._customize_seq_para({'repetitions': n_repetitions})
        list1 = [name+'1', seq_param]

        # generate second part of rf pulse
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_single_mw_pulse(name=name+'2', tau=rf_pulse2_duration, microwave_amplitude=rf_amp,
                                             microwave_frequency=rf_freq, microwave_phase=rf_phase, channel=rf_channel)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        # add sequence parameters
        seq_param2 = self._customize_seq_para({})
        list2 = [name+'2', seq_param2]

        return created_blocks, created_ensembles, list1, list2


    def _analyse_rf_pulse(self, rf_duration, rf_freq, number_of_rf_periods_per_waveform):
        rf_dict={'rf_duration': rf_duration}
        rf_dict['rf_freq'] = rf_freq
        period = 1.0/rf_freq

        # multiplicity has to be two to avoid granularity problems  with AWG
        rf_dict['period'] = self._adjust_to_samplingrate(period, 2)
        rf_dict['number_periods'] = np.floor(rf_duration / period)
        if rf_dict['number_periods'] > 0:
            remainder = rf_duration % rf_dict['period']
            remainder = self._adjust_to_samplingrate(remainder + rf_dict['period'], 2)
        else:
            remainder = self._adjust_to_samplingrate(rf_duration, 2)
        # this helps to prevent granularity problems
        rf_dict['remainder'] = remainder + 1e-12
        return rf_dict

################################### Sequence options ##########################################


    def _customize_seq_para(self, seq_para_dict):
        if 'event_trigger' not in seq_para_dict:
            seq_para_dict['event_trigger'] = 'OFF'
        if 'event_jump_to' not in seq_para_dict:
            seq_para_dict['event_jump_to'] = 0
        if 'wait_for' not in seq_para_dict:
            seq_para_dict['wait_for'] = 'OFF'
        if 'repetitions' not in seq_para_dict:
            seq_para_dict['repetitions'] = 0
        if 'go_to' not in seq_para_dict:
            seq_para_dict['go_to'] = 0
        return seq_para_dict

    def _make_sequence_continous(self, element_list):
        # change the goto key of the last list in element_list to 1
        tmp_dict = dict(element_list[-1][1])
        tmp_dict['go_to'] = 1
        element_list[-1][1] = tmp_dict
        return element_list


