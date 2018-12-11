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
        # the minimum waveform length for the Spectrum AWG is only 384 samples. Any required idle time is already
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



    def generate_laser_wait(self, name='laser_wait', laser_length=500e-9, wait_length = 1e-6):
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
        laser_element = self._get_laser_gate_element(length=laser_length, increment=0)
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
            self.log.warning('Trigger duration adjusted to minimum sequence length of {0} for samplingrate '
                             '{1}'.format(tau, self.pulse_generator_settings['sample_rate']))

        # get the trigger element
        laser_element = self._get_trigger_element(length=tau, increment=0.0, channels=digital_channel)

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

        # get the laser and idle elements
        laser_element = self._get_trigger_element(length=trigger_length, increment=0.0, channels=digital_channel)
        idle_element = self._get_idle_element(length=idle_length, increment=0.0)

        block = PulseBlock(name=name)
        block.append(laser_element)
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





################################## MW + RF pulse methods ######################################


    def generate_single_mw_pulse(self, name='MW_pulse', tau=1e-6, microwave_amplitude=1.0,
                                    microwave_frequency = 1e6, microwave_phase=0.0):

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
                                          phase=microwave_phase)

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



    def _chopped_rf_pulse(self, name, rf_duration, rf_amp, rf_freq, rf_phase):
        # analyse the rf pulse
        rf_dict = self._analyse_rf_pulse(rf_duration, rf_freq)
        # generate first part of rf pulse
        if rf_dict['number_periods'] > 0:
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_single_mw_pulse(name=name+'1', tau=rf_dict['period'], microwave_amplitude=rf_amp,
                                                 microwave_frequency=rf_freq, microwave_phase=rf_phase)
        else:
            # If there is not more than 1 period just makes this an idle with minimal length
            created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
                self.generate_idle_s3(name=name + '1', tau=self.minimum_number_of_samples/self.pulse_generator_settings['sample_rate'])
            # set it to 1 so it is repeated at least once
            rf_dict['number_periods'] = 1
        created_blocks = created_blocks_tmp
        created_ensembles = created_ensembles_tmp
        # add sequence parameters
        seq_param = self._customize_seq_para({'repetitions': rf_dict['number_periods']-1})
        list1 = [name+'1', seq_param]

        # generate second part of rf pulse
        created_blocks_tmp, created_ensembles_tmp, created_sequences_tmp = \
            self.generate_single_mw_pulse(name=name+'2', tau=rf_dict['remainder'], microwave_amplitude=rf_amp,
                                             microwave_frequency=rf_freq, microwave_phase=rf_phase)
        created_blocks += created_blocks_tmp
        created_ensembles += created_ensembles_tmp
        # add sequence parameters
        seq_param2 = self._customize_seq_para({})
        list2 = [name+'2', seq_param2]

        return created_blocks, created_ensembles, list1, list2


    def _analyse_rf_pulse(self, rf_duration, rf_freq):
        rf_dict={'rf_duration': rf_duration}
        rf_dict['rf_freq'] = rf_freq
        period = 1.0/rf_freq
        # multiplicity has to be two to avoid granularity problems  with AWG
        rf_dict['period'] = self._adjust_to_samplingrate(period, 2)
        rf_dict['number_periods'] = int(rf_duration / period) - 1
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


