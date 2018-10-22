# -*- coding: utf-8 -*-

"""
This file contains the Qudi hardware module for the M4i6631x8 series
AWG from Spectrum Instruments.

The M4i6631x8 has two 16-bit analogue channels. Bit from these analogue
channels can be sacrificed to create up to 3 digital output channels. The
assignment of the bits is flexible. We choose:
analog channel 1: bit 0-13: a_ch1
                  bit   15: d_ch1
                  bit   14: d_ch2
analog channel 2: bit 0-14: a_ch2
                  bit   15: d_ch3

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

import time
import numpy as np
import ctypes
from collections import OrderedDict
from thirdparty.spectrum_instruments.pyspcm import *

from core.module import Base, StatusVar, ConfigOption
from interface.pulser_interface import PulserInterface, PulserConstraints


class AWGSpectrumM4i6631x8(Base, PulserInterface):
    """ Dummy class for  PulseInterface

    Be careful in adjusting the method names in that class, since some of them
    are also connected to the mwsourceinterface (to give the AWG the possibility
    to act like a microwave source).
    """
    _modclass = 'awgSpectrumM4i6631x8'
    _modtype = 'hardware'

    _pulsed_file_dir = ConfigOption('pulsed_file_dir', False, missing='warn')

    activation_config = StatusVar(default=None)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self.log.info('Pulser __init__')

        self.connected = False
        self.sample_rate = 1.25e9

        # Deactivate all channels at first:
        self.channel_states = {'a_ch1': False, 'a_ch2': False,
                               'd_ch1': False, 'd_ch2': False, 'd_ch3': False}

        # for each analog channel one value
        self.amplitude_dict = {'a_ch1': 4.0, 'a_ch2': 4.0}
        self.offset_dict = {}

        # for each digital channel one value
        self.digital_high_dict = {'d_ch1': 3.3, 'd_ch2': 3.3, 'd_ch3': 3.3}
        self.digital_low_dict = {'d_ch1': 0.0, 'd_ch2': 0.0, 'd_ch3': 0.0}

        self.waveform_set = set()
        self.sequence_dict = dict()

        self.current_loaded_assets = dict()

        self.use_sequencer = True
        #self.interleave = False

        self.current_status = 0    # that means off, not running.

        # FIXME -> initialisation leaves digital output states high

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        # open card
        hCard = spcm_hOpen(create_string_buffer(b'/dev/spcm0'))
        self._hCard = hCard
        if hCard == None:
            self.log.warning('no card found...')
            return 1

        # prepare the card in a defined state
        self.reset()
        # FIXME -> initialisation leaves digital output states high

        # read type, function and sn and check for A/D card
        lCardType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCITYP, byref(lCardType))
        lSerialNumber = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCISERIALNO, byref(lSerialNumber))
        lFncType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_FNCTYPE, byref(lFncType))

        spcm_dwSetParam_i32(hCard, SPC_TIMEOUT, 25000)  # timeout 15 s

        self._active_channels = {'a_ch1': True,
                                'a_ch2': True,
                                'd_ch1': True,
                                'd_ch2': True,
                                'd_ch3': True, }
        self.set_active_channels(self._active_channels)
        self.set_analog_level({'a_ch1': 4.0, 'a_ch2': 4.0})
        self.pulser_off()
        self._loaded_asset = ''

        self.connected = True

        # dummy pulser ##############

        #self.channel_states = {'a_ch1': False, 'a_ch2': False, 'a_ch3': False,
        #                       'd_ch1': False, 'd_ch2': False, 'd_ch3': False, 'd_ch4': False,
        #                       'd_ch5': False, 'd_ch6': False, 'd_ch7': False, 'd_ch8': False}

        if self.activation_config is None:
            self.activation_config = self.get_constraints().activation_config['config_a12d123']
        elif self.activation_config not in self.get_constraints().activation_config.values():
            self.activation_config = self.get_constraints().activation_config['config_a12d123']

        #for chnl in self.activation_config:
        #    self.channel_states[chnl] = True

        ##################################



    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """

        if self._hCard != None:
            spcm_vClose(self._hCard)

        self.connected = False
        return

    def get_constraints(self):
        """
        Retrieve the hardware constrains from the Pulsing device.

        @return constraints object: object with pulser constraints as attributes.

        Provides all the constraints (e.g. sample_rate, amplitude, total_length_bins,
        channel_config, ...) related to the pulse generator hardware to the caller.

            SEE PulserConstraints CLASS IN pulser_interface.py FOR AVAILABLE CONSTRAINTS!!!

        If you are not sure about the meaning, look in other hardware files to get an impression.
        If still additional constraints are needed, then they have to be added to the
        PulserConstraints class.

        Each scalar parameter is an ScalarConstraints object defined in cor.util.interfaces.
        Essentially it contains min/max values as well as min step size, default value and unit of
        the parameter.

        PulserConstraints.activation_config differs, since it contain the channel
        configuration/activation information of the form:
            {<descriptor_str>: <channel_set>,
             <descriptor_str>: <channel_set>,
             ...}

        If the constraints cannot be set in the pulsing hardware (e.g. because it might have no
        sequence mode) just leave it out so that the default is used (only zeros).
        """
        constraints = PulserConstraints()

        constraints.sample_rate.min = 50.0e6
        constraints.sample_rate.max = 1.25e9
        constraints.sample_rate.step = 1.0e6
        constraints.sample_rate.default = 1.25e9

        constraints.a_ch_amplitude.min = 0.08 * 2
        constraints.a_ch_amplitude.max = 2.0 * 2
        constraints.a_ch_amplitude.step = 0.001 * 2
        constraints.a_ch_amplitude.default = 2.0 * 2

        constraints.a_ch_offset.min = 0.0
        constraints.a_ch_offset.max = 0.0
        constraints.a_ch_offset.step = 0.0
        constraints.a_ch_offset.default = 0.0

        constraints.d_ch_low.min = 0.0
        constraints.d_ch_low.max = 0.0
        constraints.d_ch_low.step = 0.00
        constraints.d_ch_low.default = 0.0

        constraints.d_ch_high.min = 3.3
        constraints.d_ch_high.max = 3.3
        constraints.d_ch_high.step = 0.00
        constraints.d_ch_high.default = 3.3

        constraints.waveform_length.min = 80
        constraints.waveform_length.max = 64800000
        constraints.waveform_length.step = 1
        constraints.waveform_length.default = 80

        constraints.waveform_num.min = 1
        constraints.waveform_num.max = 32000
        constraints.waveform_num.step = 1
        constraints.waveform_num.default = 1

        constraints.sequence_num.min = 1
        constraints.sequence_num.max = 4096
        constraints.sequence_num.step = 1
        constraints.sequence_num.default = 1

        constraints.subsequence_num.min = 1
        constraints.subsequence_num.max = 4000
        constraints.subsequence_num.step = 1
        constraints.subsequence_num.default = 1

        # If sequencer mode is available then these should be specified
        constraints.repetitions.min = 0
        constraints.repetitions.max = 65539
        constraints.repetitions.step = 1
        constraints.repetitions.default = 0

        #constraints.event_triggers = ['A', 'B']
        #constraints.flags = ['A', 'B', 'C', 'D']
        # Device has only one trigger and no flags
        constraints.event_triggers = ['ON']
        constraints.flags = list()

        constraints.sequence_steps.min = 0
        constraints.sequence_steps.max = 4096
        constraints.sequence_steps.step = 1
        constraints.sequence_steps.default = 0

        # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # channels. Here all possible channel configurations are stated, where only the generic
        # names should be used. The names for the different configurations can be customary chosen.
        activation_config = OrderedDict()
        activation_config['config_a12d123'] = {'a_ch1', 'a_ch2', 'd_ch1', 'd_ch2', 'd_ch3'}
        activation_config['config_a12d1'] = {'a_ch1', 'a_ch2', 'd_ch1'}
        activation_config['config_a12d13'] = {'a_ch1', 'a_ch2', 'd_ch1', 'd_ch3'}
        # Usage of channel 1 only:
        activation_config['config_a1d12'] = {'a_ch1', 'd_ch1', 'd_ch2'}
        # Usage of channel 2 only:
        activation_config['config_a2d3'] = {'a_ch2', 'd_ch3'}
        # Usage of only digital channels:
        activation_config['config_d123'] = {'d_ch1', 'd_ch2', 'd_ch3'}
        # Usage of only one analog channel:
        activation_config['config_a1'] = {'a_ch1'}
        activation_config['config_a2'] = {'a_ch2'}
        # Usage of only the analog channels:
        activation_config['config_a12'] = {'a_ch1', 'a_ch2'}
        constraints.activation_config = activation_config

        return constraints

    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:stopped, -1:error, 1:running)
        """
        #if self.current_status == 0:
        #    self.current_status = 1
        #    self.log.info('PulserDummy: Switch on the Output.')
        #   time.sleep(1)
        #    return 0
        #else:
        #    return -1

        if self.current_status == 1:
            self.log.info('Pulser already on!')
        elif self.current_status == -1:
            self.log.error('Pulser in error state, cannot switch output on!')
        elif self.current_status == 0:
            self.current_status = 1
            self.log.info('Pulser: Switch on the output.')

            # FIXME -> this assumes we always want both AWG channels activated, which is inconsistent with the full range of possible AWG configurations
            if self._spcm_dwGetParam_i32(SPC_ENABLEOUT0) == 0:
                spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 1)  # enable channel 1: a_ch1, d_ch1, d_ch2
            if self._spcm_dwGetParam_i32(SPC_ENABLEOUT1) == 0:
                spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 1)  # enable channel 2: a_ch2, d_ch3

            if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
                spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

        if self._read_out_error():
            return -1
        else:
            return 0

    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:stopped, -1:error, 1:running)
        """
        if self.current_status == 0:
            self.log.info('Pulser already off!')
        elif self.current_status == -1:
            self.log.error('Pulser in error state, cannot switch output off!')
        elif self.current_status == 1:
            self.current_status = 0
            self.log.info('Pulser: Switch off the output.')

            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 0)  # disable channel 1: a_ch1, d_ch1, d_ch2
            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 0)  # disable channel 2: a_ch2, d_ch3

        if self._read_out_error():
            return -1
        else:
            return 0

    def write_waveform(self, name, analog_samples, digital_samples, is_first_chunk, is_last_chunk,
                       total_number_of_samples):
        """
        Write a new waveform or append samples to an already existing waveform on the device memory.
        The flags is_first_chunk and is_last_chunk can be used as indicator if a new waveform should
        be created or if the write process to a waveform should be terminated.

        NOTE: All sample arrays in analog_samples and digital_samples must be of equal length!

        @param str name: the name of the waveform to be created/append to
        @param dict analog_samples: keys are the generic analog channel names (i.e. 'a_ch1') and
                                    values are 1D numpy arrays of type float32 containing the
                                    voltage samples.
        @param dict digital_samples: keys are the generic digital channel names (i.e. 'd_ch1') and
                                     values are 1D numpy arrays of type bool containing the marker
                                     states.
        @param bool is_first_chunk: Flag indicating if it is the first chunk to write.
                                    If True this method will create a new empty wavveform.
                                    If False the samples are appended to the existing waveform.
        @param bool is_last_chunk:  Flag indicating if it is the last chunk to write.
                                    Some devices may need to know when to close the appending wfm.
        @param int total_number_of_samples: The number of sample points for the entire waveform
                                            (not only the currently written chunk)

        @return (int, list): Number of samples written (-1 indicates failed process) and list of
                             created waveform names
        """
        # TODO -> change from PulserDummy
        waveforms = list()

        # Sanity checks
        if len(analog_samples) > 0:
            number_of_samples = len(analog_samples[list(analog_samples)[0]])
        elif len(digital_samples) > 0:
            number_of_samples = len(digital_samples[list(digital_samples)[0]])
        else:
            self.log.error('No analog or digital samples passed to write_waveform method in dummy '
                           'pulser.')
            return -1, waveforms

        for chnl, samples in analog_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of sample arrays for different channels in dummy '
                               'pulser.')
                return -1, waveforms
        for chnl, samples in digital_samples.items():
            if len(samples) != number_of_samples:
                self.log.error('Unequal length of sample arrays for different channels in dummy '
                               'pulser.')
                return -1, waveforms

        # Determine if only digital samples are active. In that case each channel will get a
        # waveform. Otherwise only the analog channels will have a waveform with digital channel
        # samples included (as it is the case in Tektronix and Keysight AWGs).
        # Simulate a 1Gbit/s transfer speed. Assume each analog waveform sample is 5 bytes large
        # (4 byte float and 1 byte marker bitmask). Assume each digital waveform sample is 1 byte.
        if len(analog_samples) > 0:
            for chnl in analog_samples:
                waveforms.append(name + chnl[1:])
                time.sleep(number_of_samples * 5 * 8 / 1024 ** 3)
        else:
            for chnl in digital_samples:
                waveforms.append(name + chnl[1:])
                time.sleep(number_of_samples * 8 / 1024 ** 3)

        self.waveform_set.update(waveforms)

        self.log.info('Waveforms with nametag "{0}" directly written on dummy pulser.'.format(name))
        return number_of_samples, waveforms

    def write_sequence(self, name, sequence_parameter_list):
        """
        Write a new sequence on the device memory.

        @param name: str, the name of the waveform to be created/append to
        @param sequence_parameter_list: list, contains the parameters for each sequence step and
                                        the according waveform names.

        @return: int, number of sequence steps written (-1 indicates failed process)
        """
        # TODO -> change from PulserDummy
        # Check if all waveforms are present on virtual device memory
        for waveform_tuple, param_dict in sequence_parameter_list:
            for waveform in waveform_tuple:
                if waveform not in self.waveform_set:
                    self.log.error('Failed to create sequence "{0}" due to waveform "{1}" not '
                                   'present in device memory.'.format(name, waveform))
                    return -1

        if name in self.sequence_dict:
            del self.sequence_dict[name]

        self.sequence_dict[name] = len(sequence_parameter_list[0][0])
        time.sleep(1)

        self.log.info('Sequence with name "{0}" directly written on dummy pulser.'.format(name))
        return len(sequence_parameter_list)

    def get_waveform_names(self):
        """ Retrieve the names of all uploaded waveforms on the device.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """
        # TODO -> change from PulserDummy
        return list(self.waveform_set)

    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """
        # TODO -> change from PulserDummy
        return list(self.sequence_dict)

    def delete_waveform(self, waveform_name):
        """ Delete the waveform with name "waveform_name" from the device memory.

        @param str waveform_name: The name of the waveform to be deleted
                                  Optionally a list of waveform names can be passed.

        @return list: a list of deleted waveform names.
        """
        # TODO -> change from PulserDummy
        if isinstance(waveform_name, str):
            waveform_name = [waveform_name]

        deleted_waveforms = list()
        for waveform in waveform_name:
            if waveform in self.waveform_set:
                self.waveform_set.remove(waveform)
                deleted_waveforms.append(waveform)

        return deleted_waveforms

    def delete_sequence(self, sequence_name):
        """ Delete the sequence with name "sequence_name" from the device memory.

        @param str sequence_name: The name of the sequence to be deleted
                                  Optionally a list of sequence names can be passed.

        @return list: a list of deleted sequence names.
        """
        # TODO -> change from PulserDummy
        if isinstance(sequence_name, str):
            sequence_name = [sequence_name]

        deleted_sequences = list()
        for sequence in sequence_name:
            if sequence in self.sequence_dict:
                del self.sequence_dict[sequence]
                deleted_sequences.append(sequence)

        return deleted_sequences

    def load_waveform(self, load_dict):
        """ Loads a waveform to the specified channel of the pulsing device.
        For devices that have a workspace (i.e. AWG) this will load the waveform from the device
        workspace into the channel.
        For a device without mass memory this will make the waveform/pattern that has been
        previously written with self.write_waveform ready to play.

        @param load_dict:  dict|list, a dictionary with keys being one of the available channel
                                      index and values being the name of the already written
                                      waveform to load into the channel.
                                      Examples:   {1: rabi_ch1, 2: rabi_ch2} or
                                                  {1: rabi_ch2, 2: rabi_ch1}
                                      If just a list of waveform names if given, the channel
                                      association will be invoked from the channel
                                      suffix '_ch1', '_ch2' etc.

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """
        # TODO -> change from PulserDummy
        if isinstance(load_dict, list):
            new_dict = dict()
            for waveform in load_dict:
                channel = int(waveform.rsplit('_ch', 1)[1])
                new_dict[channel] = waveform
            load_dict = new_dict

        # Determine if the device is purely digital and get all active channels
        analog_channels = [chnl for chnl in self.activation_config if chnl.startswith('a')]
        digital_channels = [chnl for chnl in self.activation_config if chnl.startswith('d')]
        pure_digital = len(analog_channels) == 0

        # Check if waveforms are present in virtual dummy device memory and specified channels are
        # active. Create new load dict.
        new_loaded_assets = dict()
        for channel, waveform in load_dict.items():
            if waveform not in self.waveform_set:
                self.log.error('Loading failed. Waveform "{0}" not found on device memory.'
                               ''.format(waveform))
                return self.current_loaded_assets
            if pure_digital:
                if 'd_ch{0:d}'.format(channel) not in digital_channels:
                    self.log.error('Loading failed. Digital channel {0:d} not active.'
                                   ''.format(channel))
                    return self.current_loaded_assets
            else:
                if 'a_ch{0:d}'.format(channel) not in analog_channels:
                    self.log.error('Loading failed. Analog channel {0:d} not active.'
                                   ''.format(channel))
                    return self.current_loaded_assets
            new_loaded_assets[channel] = waveform
        self.current_loaded_assets = new_loaded_assets
        return self.get_loaded_assets()

    def load_sequence(self, sequence_name):
        """ Loads a sequence to the channels of the device in order to be ready for playback.
        For devices that have a workspace (i.e. AWG) this will load the sequence from the device
        workspace into the channels.

        @param sequence_name:  str, name of the sequence to load

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """
        # TODO -> change from PulserDummy
        if sequence_name not in self.sequence_dict:
            self.log.error('Sequence loading failed. No sequence with name "{0}" found on device '
                           'memory.'.format(sequence_name))
            return self.get_loaded_assets()

        # Determine if the device is purely digital and get all active channels
        analog_channels = sorted(chnl for chnl in self.activation_config if chnl.startswith('a'))
        digital_channels = sorted(chnl for chnl in self.activation_config if chnl.startswith('d'))
        pure_digital = len(analog_channels) == 0

        if pure_digital and len(digital_channels) != self.sequence_dict[sequence_name]:
            self.log.error('Sequence loading failed. Number of active digital channels ({0:d}) does'
                           ' not match the number of tracks in the sequence ({1:d}).'
                           ''.format(len(digital_channels), self.sequence_dict[sequence_name]))
            return self.get_loaded_assets()
        if not pure_digital and len(analog_channels) != self.sequence_dict[sequence_name]:
            self.log.error('Sequence loading failed. Number of active analog channels ({0:d}) does'
                           ' not match the number of tracks in the sequence ({1:d}).'
                           ''.format(len(analog_channels), self.sequence_dict[sequence_name]))
            return self.get_loaded_assets()

        new_loaded_assets = dict()
        if pure_digital:
            for track_index, chnl in enumerate(digital_channels):
                chnl_num = int(chnl.split('ch')[1])
                new_loaded_assets[chnl_num] = '{0}_{1:d}'.format(sequence_name, track_index)
        else:
            for track_index, chnl in enumerate(analog_channels):
                chnl_num = int(chnl.split('ch')[1])
                new_loaded_assets[chnl_num] = '{0}_{1:d}'.format(sequence_name, track_index)

        self.current_loaded_assets = new_loaded_assets
        return self.get_loaded_assets()

    def get_loaded_assets(self):
        """
        Retrieve the currently loaded asset names for each active channel of the device.
        The returned dictionary will have the channel numbers as keys.
        In case of loaded waveforms the dictionary values will be the waveform names.
        In case of a loaded sequence the values will be the sequence name appended by a suffix
        representing the track loaded to the respective channel (i.e. '<sequence_name>_1').

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel,
                             string describing the asset type ('waveform' or 'sequence')
        """
        # TODO -> change from PulserDummy
        # Determine if it's a waveform or a sequence
        asset_type = None
        for asset_name in self.current_loaded_assets.values():
            if 'ch' in asset_name.rsplit('_', 1)[1]:
                current_type = 'waveform'
            else:
                current_type = 'sequence'

            if asset_type is None or asset_type == current_type:
                asset_type = current_type
            else:
                self.log.error('Unable to determine loaded asset type. Mixed naming convention '
                               'assets loaded (waveform and sequence tracks).')
                return dict(), ''

        return self.current_loaded_assets, asset_type

    def clear_all(self):
        """ Clears all loaded waveform from the pulse generators RAM.

        @return int: error code (0:OK, -1:error)

        Unused for digital pulse generators without storage capability
        (PulseBlaster, FPGA).
        """
        # TODO -> change from PulserDummy
        self.current_loaded_assets = dict()
        self.waveform_set = set()
        self.sequence_dict = dict()
        return 0

    def get_status(self):
        """ Retrieves the status of the pulsing hardware

        @return (int, dict): inter value of the current status with the
                             corresponding dictionary containing status
                             description for all the possible status variables
                             of the pulse generator hardware
        @return (int, dict): tuple with an integer value of the current status and a corresponding
                             dictionary containing status description for all the possible status
                             variables of the pulse generator hardware.

                            _spcm_dwGetParam_i32(SPC_M2STATUS):
                             1h Acquisition modes only: the pretrigger area has been filled.
                             2h The first trigger has been detected.
                             4h The card has finished its run and is ready.
                             8h Multi/ABA/Gated acquisition of M4i/M4x/M2p only:
                                the pretrigger area of one segment has been filled.)
        """
        # FIXME -> finish this!!
        status = self._spcm_dwGetParam_i32(SPC_M2STATUS)
        self.log.info('M4i6631x8 status: {0:x} (see manual for more information)'.format(status))
        return 0, {0: ''}
        status_dic = {}

        status_dic[-1] = 'Failed Request or Communication'
        status_dic[0] = 'Device has stopped, but can receive commands.'
        status_dic[1] = 'Device is active and running.'
        # All the other status messages should have higher integer values
        # than 1.
        return self.current_status, status_dic

    def get_sample_rate(self):
        """ Get the sample rate of the pulse generator hardware

        @return float: The current sample rate of the device (in Hz)

        Do not return a saved sample rate in a class variable, but instead
        retrieve the current sample rate directly from the device.
        """
        # TODO -> check this works
        sample_rate = self._spcm_dwGetParam_i32(SPC_SAMPLERATE)
        return float(sample_rate)

    def set_sample_rate(self, sample_rate):
        """ Set the sample rate of the pulse generator hardware

        @param float sample_rate: The sampling rate to be set (in Hz)

        @return float: the sample rate returned from the device.

        Note: After setting the sampling rate of the device, retrieve it again
              for obtaining the actual set value and use that information for
              further processing.
        """
        # TODO -> check this works
        constraint = self.get_constraints().sample_rate
        if sample_rate > constraint.max:
            self.sample_rate = constraint.max
        elif sample_rate < constraint.min:
            self.sample_rate = constraint.min
        else:
            self.sample_rate = sample_rate

        spcm_dwSetParam_i32(self._hCard, SPC_SAMPLERATE, int32(int(sample_rate)))
        self._read_out_error()
        return self.get_sample_rate()

    def get_analog_level(self, amplitude=None, offset=None):
        """ Retrieve the analog amplitude and offset of the provided channels.

        @param list amplitude: optional, if a specific amplitude value (in Volt
                               peak to peak, i.e. the full amplitude) of a
                               channel is desired.
        @param list offset: optional, if a specific high value (in Volt) of a
                            channel is desired.

        @return dict: with keys being the generic string channel names and items
                      being the values for those channels. Amplitude is always
                      denoted in Volt-peak-to-peak and Offset in (absolute)
                      Voltage.

        Note: Do not return a saved amplitude and/or offset value but instead
              retrieve the current amplitude and/or offset directly from the
              device.

        If no entries provided then the levels of all channels where simply
        returned. If no analog channels provided, return just an empty dict.
        Example of a possible input:
            amplitude = ['a_ch1','a_ch4'], offset =[1,3]
        to obtain the amplitude of channel 1 and 4 and the offset
            {'a_ch1': -0.5, 'a_ch4': 2.0} {'a_ch1': 0.0, 'a_ch3':-0.75}
        since no high request was performed.

        The major difference to digital signals is that analog signals are
        always oscillating or changing signals, otherwise you can use just
        digital output. In contrast to digital output levels, analog output
        levels are defined by an amplitude (here total signal span, denoted in
        Voltage peak to peak) and an offset (a value around which the signal
        oscillates, denoted by an (absolute) voltage).

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        # TODO -> check this works

        ampl_ch1 = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_AMP0, byref(ampl_ch1))  # channel0 amplitude in mV

        ampl_ch2 = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_AMP1, byref(ampl_ch2))  # channel1 amplitude in mV

        ampl_dict = {'a_ch1': float(ampl_ch1.value * 2 / 1000), 'a_ch2': float(ampl_ch2.value * 2 / 1000)}

        if amplitude is not None:
            ampl_dict = {key: ampl_dict[key] for key in amplitude}
        if offset is not None:
            offset_dict = {key: 0.0 for key in offset}
        else:
            offset_dict = {'a_ch1': 0.0, 'a_ch2': 0.0}

        self._read_out_error()
        return ampl_dict, offset_dict

    def set_analog_level(self, amplitude=None, offset=None):
        """ Set amplitude and/or offset value of the provided analog channel.

        @param dict amplitude: dictionary, with key being the channel and items
                               being the amplitude values (in Volt peak to peak,
                               i.e. the full amplitude) for the desired channel.
        @param dict offset: dictionary, with key being the channel and items
                            being the offset values (in absolute volt) for the
                            desired channel.

        @return (dict, dict): tuple of two dicts with the actual set values for
                              amplitude and offset.

        If nothing is passed then the command will return two empty dicts.

        Note: After setting the analog and/or offset of the device, retrieve
              them again for obtaining the actual set value(s) and use that
              information for further processing.

        The major difference to digital signals is that analog signals are
        always oscillating or changing signals, otherwise you can use just
        digital output. In contrast to digital output levels, analog output
        levels are defined by an amplitude (here total signal span, denoted in
        Voltage peak to peak) and an offset (a value around which the signal
        oscillates, denoted by an (absolute) voltage).

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        # TODO -> check if this works
        # TODO -> change to be robust against out-of-bounds inputs. Base off tektronix_awg70k
        # if sample_rate > constraint.max:
        #    self.sample_rate = constraint.max
        # elif sample_rate < constraint.min:
        #    self.sample_rate = constraint.min
        # else:
        #    self.sample_rate = sample_rate


        #for a_ch, amp in amplitude.items():
        #    self.amplitude_dict[a_ch] = amp
#
        #for a_ch, off in offset.items():
        #    self.offset_dict[a_ch] = off


        if amplitude is None:
            amplitude = dict()
        if offset is not None:
            self.log.warning('Setting analog offset values is not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')

            #if self.sample_rate > constraint.max:
            #    self.sample_rate = constraint.max
            #elif sample_rate < constraint.min:
            #    self.sample_rate = constraint.min
            #else:
            #    self.sample_rate = sample_rate

        constraints = self.get_constraints()
        if 'a_ch1' in amplitude.keys():
            if (amplitude['a_ch1'] >= constraints.a_ch_amplitude.min) and\
                    (amplitude['a_ch1'] <= constraints.a_ch_amplitude.max):
                # channel0 amplitude in mV
                spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(int(amplitude['a_ch1']*1000/2)))
            else:
                self.log.warning('Amplitude voltage level for the analog output channel 1 is out of constraints for the'
                                 'Spectrum M4i AWG series!\nMethod call will be ignored.')
        if 'a_ch2' in amplitude.keys():
            if (amplitude['a_ch2'] >= constraints.a_ch_amplitude.min) and\
                    (amplitude['a_ch2'] <= constraints.a_ch_amplitude.max):
                # channel0 amplitude in mV
                spcm_dwSetParam_i32(self._hCard, SPC_AMP1, int32(int(amplitude['a_ch2']*1000/2)))
            else:
                self.log.warning('Amplitude voltage level for the analog output channel 2 is out of constraints for the'
                                 'Spectrum M4i AWG series!\nMethod call will be ignored.')

        self._read_out_error()
        return self.get_analog_level()

    def get_digital_level(self, low=None, high=None):
        """ Retrieve the digital low and high level of the provided channels.

        @param list low: optional, if a specific low value (in Volt) of a
                         channel is desired.
        @param list high: optional, if a specific high value (in Volt) of a
                          channel is desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel
                               number and items being the values for those
                               channels. Both low and high value of a channel is
                               denoted in (absolute) Voltage.

        Note: Do not return a saved low and/or high value but instead retrieve
              the current low and/or high value directly from the device.

        If no entries provided then the levels of all channels where simply
        returned. If no digital channels provided, return just an empty dict.

        Example of a possible input:
            low = ['d_ch1', 'd_ch4']
        to obtain the low voltage values of digital channel 1 an 4. A possible
        answer might be
            {'d_ch1': -0.5, 'd_ch4': 2.0} {}
        since no high request was performed.

        The major difference to analog signals is that digital signals are
        either ON or OFF, whereas analog channels have a varying amplitude
        range. In contrast to analog output levels, digital output levels are
        defined by a voltage, which corresponds to the ON status and a voltage
        which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        # TODO -> check this works

        constraints = self.get_constraints()
        voltage_low = constraints.d_ch_low.default
        voltage_high = constraints.d_ch_high.default
        low_dict = {'d_ch1': voltage_low, 'd_ch2': voltage_low, 'd_ch3': voltage_low}
        high_dict = {'d_ch1': voltage_high, 'd_ch2': voltage_high, 'd_ch3': voltage_high}

        if low is not None:
            low_dict = {key: low_dict[key] for key in low}
        if high is not None:
            high_dict = {key: high_dict[key] for key in high}

        return low_dict, high_dict

    def set_digital_level(self, low=None, high=None):
        """ Set low and/or high value of the provided digital channel.

        @param dict low: dictionary, with key being the channel and items being
                         the low values (in volt) for the desired channel.
        @param dict high: dictionary, with key being the channel and items being
                         the high values (in volt) for the desired channel.

        @return (dict, dict): tuple of two dicts where first dict denotes the
                              current low value and the second dict the high
                              value.

        If nothing is passed then the command will return two empty dicts.

        Note: After setting the high and/or low values of the device, retrieve
              them again for obtaining the actual set value(s) and use that
              information for further processing.

        The major difference to analog signals is that digital signals are
        either ON or OFF, whereas analog channels have a varying amplitude
        range. In contrast to analog output levels, digital output levels are
        defined by a voltage, which corresponds to the ON status and a voltage
        which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between
        (amplitude, offset) and (value high, value low)!
        """
        self.log.warning('Setting the digital voltage levels is not available for the Spectrum M4i AWG series!\n'
                         'Method call will be ignored.')
        return self.get_digital_level()

    def get_active_channels(self, ch=None):
        """ Get the active channels of the pulse generator hardware.

        @param list ch: optional, if specific analog or digital channels are
                        needed to be asked without obtaining all the channels.

        @return dict:  where keys denoting the channel number and items boolean
                       expressions whether channel are active or not.

        Example for an possible input (order is not important):
            ch = ['a_ch2', 'd_ch2', 'a_ch1', 'd_ch5', 'd_ch1']
        then the output might look like
            {'a_ch2': True, 'd_ch2': False, 'a_ch1': False, 'd_ch5': True, 'd_ch1': False}

        If no parameters are passed to this method all channels will be asked
        for their setting.
        """
        # TODO -> check this works
        #if ch is None:
        #    ch = []
#
        #active_ch = {}
#
        #if not ch:
        #    active_ch = self.channel_states
#
        #else:
        #    for channel in ch:
        #        active_ch[channel] = self.channel_states[channel]
#
        #return active_ch

        # a_ch1 and a_ch2
        a_channel = self._spcm_dwGetParam_i32(SPC_CHENABLE)

        self._active_channels['a_ch1'] = bool(a_channel & CHANNEL0)
        self._active_channels['a_ch2'] = bool(a_channel & CHANNEL1)

        # digital channels
        # FIXME -> configfile
        # current implementation:
        # analog channel 1: bit 0-13: a_ch1
        #                   bit   15: d_ch1
        #                   bit   14: d_ch2
        # analog channel 2: bit 0-14: a_ch2
        #                   bit   15: d_ch3

        # d_ch1
        x0_mode = self._spcm_dwGetParam_i32(SPCM_X0_MODE)
        self._active_channels['d_ch1'] \
            = bool(x0_mode & SPCM_XMODE_DIGOUT) and bool(x0_mode & SPCM_XMODE_DIGOUTSRC_CH0)

        # d_ch2
        x1_mode = self._spcm_dwGetParam_i32(SPCM_X1_MODE)
        self._active_channels['d_ch2'] \
            = bool(x1_mode & SPCM_XMODE_DIGOUT) and bool(x1_mode & SPCM_XMODE_DIGOUTSRC_CH0)

        # d_ch3
        x2_mode = self._spcm_dwGetParam_i32(SPCM_X2_MODE)
        self._active_channels['d_ch3'] \
            = bool(x2_mode & SPCM_XMODE_DIGOUT) and bool(x2_mode & SPCM_XMODE_DIGOUTSRC_CH1)

        if ch is None:
            return self._active_channels
        else:
            return {key: self._active_channels[key] for key in ch.keys()}

    def set_active_channels(self, ch=None):
        """
        Set the active/inactive channels for the pulse generator hardware.
        The state of ALL available analog and digital channels will be returned
        (True: active, False: inactive).
        The actually set and returned channel activation must be part of the available
        activation_configs in the constraints.
        You can also activate/deactivate subsets of available channels but the resulting
        activation_config must still be valid according to the constraints.
        If the resulting set of active channels can not be found in the available
        activation_configs, the channel states must remain unchanged.

        @param dict ch: dictionary with keys being the analog or digital string generic names for
                        the channels (i.e. 'd_ch1', 'a_ch2') with items being a boolean value.
                        True: Activate channel, False: Deactivate channel

        @return dict: with the actual set values for ALL active analog and digital channels

        If nothing is passed then the command will simply return the unchanged current state.

        Note: After setting the active channels of the device, use the returned dict for further
              processing.

        Example for possible input:
            ch={'a_ch2': True, 'd_ch1': False, 'd_ch3': True, 'd_ch4': True}
        to activate analog channel 2 digital channel 3 and 4 and to deactivate
        digital channel 1. All other available channels will remain unchanged.
        """
        # TODO -> check this works
        #if ch is None:
        #    ch = {}
        #old_activation = self.channel_states.copy()
        #for channel in ch:
        #    self.channel_states[channel] = ch[channel]
#
        #active_channel_set = {chnl for chnl, is_active in self.channel_states.items() if is_active}
        #if active_channel_set not in self.get_constraints().activation_config.values():
        #    self.log.error('Channel activation to be set not found in constraints.\n'
        #                   'Channel activation unchanged.')
        #    self.channel_states = old_activation
        #else:
        #    self.activation_config = active_channel_set
#
        #return self.get_active_channels(ch=list(ch))

        if ch is None:
            return self.get_active_channels()

        for channel in ch:
            if channel in self._active_channels:
                self._active_channels[channel] = ch[channel]
            else:
                self.log.error('Trying to (de)activate channel "{0}". This channel is not present '
                               'in AWG. Setting channels aborted.'.format(channel))

        # analog channels
        channel0 = CHANNEL0 * (self._active_channels['a_ch1'] or self._active_channels['d_ch1'] or
                               self._active_channels['d_ch2'])
        channel1 = CHANNEL1 * (self._active_channels['a_ch2'] or self._active_channels['d_ch3'])
        spcm_dwSetParam_i32(self._hCard, SPC_CHENABLE, channel0 | channel1)

        # digital channels
        # FIXME -> configfile
        # current implementation:
        # analog channel 1: bit 0-13: a_ch1
        #                   bit   15: d_ch1
        #                   bit   14: d_ch2
        # analog channel 2: bit 0-14: a_ch2
        #                   bit   15: d_ch3

        # d_ch1
        if self._active_channels['d_ch1']:
            spcm_dwSetParam_i32(
                self._hCard, SPCM_X0_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT15)
        else:
            spcm_dwSetParam_i32(self._hCard, SPCM_X0_MODE, SPCM_XMODE_DISABLE)
        # d_ch2
        if self._active_channels['d_ch2']:
            spcm_dwSetParam_i32(
                self._hCard, SPCM_X1_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT14)
        else:
            spcm_dwSetParam_i32(self._hCard, SPCM_X1_MODE, SPCM_XMODE_DISABLE)
        # d_ch3
        if self._active_channels['d_ch3']:
            spcm_dwSetParam_i32(
                self._hCard, SPCM_X2_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH1 | SPCM_XMODE_DIGOUTSRC_BIT15)
        else:
            spcm_dwSetParam_i32(self._hCard, SPCM_X2_MODE, SPCM_XMODE_DISABLE)

        self._active_channels = self.get_active_channels()

        self._read_out_error()
        return self._active_channels

    def get_interleave(self):
        """ Check whether Interleave is ON or OFF in AWG.

        @return bool: True: ON, False: OFF

        Unused for pulse generator hardware other than an AWG.
        """
        self.log.info('Interleave mode not available for the Spectrum M4i AWG series!')
        return False

    def set_interleave(self, state=False):
        """ Turns the interleave of an AWG on or off.

        @param bool state: The state the interleave should be set to
                           (True: ON, False: OFF)

        @return bool: actual interleave status (True: ON, False: OFF)

        Note: After setting the interleave of the device, retrieve the
              interleave again and use that information for further processing.

        Unused for pulse generator hardware other than an AWG.
        """
        if state:
            self.log.warning('Interleave mode not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')
        return False

    def write(self, command):
        """ Sends a command string to the device.

        @param string command: string containing the command

        @return int: error code (0:OK, -1:error)
        """
        # TODO -> change from PulserDummy
        self.log.info('It is so nice that you talk to me and told me "{0}"; '
                      'as a dummy it is very dull out here! :) '.format(command))
        return 0

    def query(self, question):
        """ Asks the device a 'question' and receive and return an answer from it.

        @param string question: string containing the command

        @return string: the answer of the device to the 'question' in a string
        """
        # TODO -> change from PulserDummy
        self.log.info('Dude, I\'m a dummy! Your question \'{0}\' is way too '
                      'complicated for me :D !'.format(question))
        return 'I am a dummy!'

    def reset(self):
        """ Reset the device.

        @return int: error code (0:OK, -1:error)
        """
        #self.__init__()
        # FIXME -> reset leaves digital output states high
        self._read_out_error()
        if spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_RESET):
            self._read_out_error()
            self.log.error('AWG reset error')
            return -1
        else:
            self.log.info('AWG reset successful')
            self.connected = True
            return 0

    def has_sequence_mode(self):
        """ Asks the pulse generator whether sequence mode exists.

        @return: bool, True for yes, False for no.
        """
        available_modes = int32(0)
        spcm_dwSetParam_i32(self._hCard, SPC_AVAILCARDMODES, byref(available_modes))
        self._read_out_error()
        if available_modes.value & SPC_REP_STD_SEQUENCE:
            return True
        else:
            return False

    def force_trigger(self, wait=False):
        """
        Sends a trigger to the awg.

        @param wait bool: waits until the card has completed the current run (until timeout)

        @return int: error code (0:OK, -1:error)
        """
        if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
            return -1
        if wait:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER | M2CMD_CARD_WAITREADY)
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER)

        if self._read_out_error():
            return -1
        else:
            return 0

    def _read_out_error(self):
        """checks the error state and prints it out. Errors must be read out before the AWG
        can accept further commands"""
        spcm_dwGetParam_i32
        errortext = (ctypes.c_char*ERRORTEXTLEN)()
        err = spcm_dwGetErrorInfo_i32(self._hCard, None, None, errortext)
        if err:
            self.log.warning(errortext.value)
        return err

    def _spcm_dwGetParam_i32(self, lRegister):
        """
        Reads an internal register.

        @param lRegister (int): register, that should be read out

        @return int: (-1: error , else parameter value)
        """
        value = int32(0)
        spcm_dwGetParam_i32(self._hCard, lRegister, byref(value))
        if self._read_out_error():
           return -1
        else:
            return value.value

    def _spcm_dwSetParam_i32(self, lRegister, plValue):
        """
        Sets an internal register.

        @param lRegister (int): register, that should be set
               plValue   (int): value of the register

        @return int: (-1: error , else parameter value)
        """
        spcm_dwSetParam_i32(self._hCard, lRegister, int32(plValue))
        if self._read_out_error():
           return -1
        else:
            return 0

    def _create_combined_buffer_data(self, llMemSamples, channel_data):
        """
        digital channels are encoded in analog samples for synchonous readout
        -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)

        @param llMemSamples: number of samples for all channels
               channel_data: dictionary with channel keys and channel data
                             e.g. {'a_ch1': np.array([...]), 'd_ch1': np.array([])}

        @return (ptr16,  int64): pointer to buffer data, buffersize in bytes
        """

        # FIXME -> configfile
        # current implementation:
        # analog channel 1: bit 0-13: a_ch1
        #                   bit   15: d_ch1
        #                   bit   14: d_ch2
        # analog channel 2: bit 0-14: a_ch2
        #                   bit   15: d_ch3

        active_channels = self._active_channels

        for key in channel_data.keys():
            if channel_data[key].size != llMemSamples.value:
                self.log.error('Channel data for channel {} has the wrong size!'.format(key))
                return -1

        for key in active_channels.keys():
            if key not in channel_data.keys():
                channel_data[key] = np.zeros(llMemSamples.value, dtype=np.dtype(np.int16))

        if active_channels['d_ch2']:
            a_ch1_data = ((channel_data['a_ch1'].astype(np.uint16) >> 2) |
                          (channel_data['d_ch1'].astype(bool) << 15) |
                          (channel_data['d_ch2'].astype(bool) << 14))
        elif active_channels['d_ch1']:
            a_ch1_data = ((channel_data['a_ch1'].astype(np.unint16) >> 1) |
                          (channel_data['d_ch1'].astype(bool) << 15))

        if active_channels['d_ch3']:
            a_ch2_data = ((channel_data['a_ch2'].astype(np.uint16) >> 1) |
                          (channel_data['d_ch3'].astype(np.uint16) << 15))

        if active_channels['a_ch1'] and active_channels['a_ch2']:
            combined_channel_data = np.vstack((a_ch1_data, a_ch2_data)).ravel('F').astype(np.int16)
        elif active_channels['a_ch1']:
            combined_channel_data = a_ch1_data.astype(np.int16)
        elif active_channels['a_ch2']:
            combined_channel_data = a_ch2_data.astype(np.int16)
        else:
            self.log.error('No channel activated.')
            return -1

        # setup software buffer
        chcount = self._spcm_dwGetParam_i32(SPC_CHCOUNT)
        lBytesPerSample = self._spcm_dwGetParam_i32(SPC_MIINST_BYTESPERSAMPLE)

        qwBufferSize = uint64(llMemSamples.value * lBytesPerSample * chcount)
        pvBuffer = create_string_buffer(qwBufferSize.value)

        # calculate the data
        pnBuffer = cast(pvBuffer, ptr16)
        for i in range(0, llMemSamples.value * chcount, 1):
            pnBuffer[i] = int16(combined_channel_data[i])

        return pvBuffer, qwBufferSize

    def _set_trigger(self, software_trigger = False):
        """
        Sets up the trigger setting for the awg.

        @param software_trigger bool: enables (True) or disables (False) the software trigger

        @return int: error code (0:OK, -1:error)
        """
        # Setting up the trigger
        if software_trigger:
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)

        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, SPC_TMASK_NONE)

        if self._read_out_error():
            return -1
        else:
            return 0

    def set_digital_async_on(self,channel):
        """
        Turns on one of the Multi Purpose I/O Lines as a digital output

        @param channel int: select channel 0,1,2

        @return int: error code (0:OK, -1:error)
        """
        if channel == 0:
            spcm_dwSetParam_i32(self._hCard, SPCM_X0_MODE, SPCM_XMODE_ASYNCOUT)
        elif channel == 1:
            spcm_dwSetParam_i32(self._hCard, SPCM_X1_MODE, SPCM_XMODE_ASYNCOUT)
        elif channel == 2:
            spcm_dwSetParam_i32(self._hCard, SPCM_X2_MODE, SPCM_XMODE_ASYNCOUT)


        if self._read_out_error():
            return -1
        else:
            return 0

    def set_digital_async_off(self, channel):
        """
        Turns off one of the Multi Purpose I/O Lines as a digital output

        @param channel int: select channel 0,1,2

        @return int: error code (0:OK, -1:error)
        """



        if self._read_out_error():
            return -1
        else:
            return 0
    def _testfunction(self, channel, amplitude, offset, numberofsamples=int64(KILO_B(64))):
        """we try to produce a sine wave"""

        # setup the mode
        #qwChEnable = uint64(channel+1)
        llMemSamples = numberofsamples
        self.log.info(llMemSamples)
        llLoops = int64(0)  # loop continuously
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SINGLERESTART) #continuous
        spcm_dwSetParam_i64(self._hCard, SPC_CHENABLE, CHANNEL0)                  #channel enabled 1
        spcm_dwSetParam_i64(self._hCard, SPC_MEMSIZE, llMemSamples)               #replay lenth
        spcm_dwSetParam_i64(self._hCard, SPC_LOOPS, llLoops)                      #number of repetitions

        spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 1)                       #enable analog output 1

        lSetChannels = int32(0)                                                 #number of activated channels
        spcm_dwGetParam_i32(self._hCard, SPC_CHCOUNT, byref(lSetChannels))
        lBytesPerSample = int32(0)                                              #bytes per sample
        spcm_dwGetParam_i32(self._hCard, SPC_MIINST_BYTESPERSAMPLE, byref(lBytesPerSample))

        # setup the trigger mode
        # (SW trigger, no output)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, 0)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, 0)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, 0)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, 0)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, 0)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, 0)

        spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(1000)) #channel0 amplitude in mV
        spcm_dwSetParam_i32(self._hCard, SPC_AMP1, int32(1000)) #channel1 amplitude in mV

        # setup software buffer
        qwBufferSize = uint64(llMemSamples.value * lBytesPerSample.value * lSetChannels.value)
        pvBuffer = create_string_buffer(qwBufferSize.value)

        self.log.info(1)
        # calculate the data
        pnBuffer = cast(pvBuffer, ptr16)
        for i in range(0, llMemSamples.value, 1):
            if channel == 0:
                pnBuffer[i] = i
            elif channel == 1:
                pnBuffer[i] = int((np.sin(i/int(llMemSamples.value)*10e6/1.25e9 * 2*np.pi) + offset)*amplitude)
            else:
                print('error')

        # we define the buffer for transfer and start the DMA transfer
        print("Starting the DMA transfer and waiting until data is in board memory\n")
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        print("... data has been transferred to board memory\n")

        # We'll start and wait until the card has finished or until a timeout occurs
        spcm_dwSetParam_i32(self._hCard, SPC_TIMEOUT, 15000)
        print("Starting the card and waiting for ready interrupt\n(continuous and single restart will have timeout)\n")
        dwError = spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START)
        #dwError = spcm_dwSetParam_i32(self._hCard, SPC_M2CMD,
        #                             M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY)
        #if dwError == ERR_TIMEOUT:
        #    spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)
        #    print("error: ", dwError)

    def _testfunction2(self, frequency_1 = 100e6, frequency_2 = 200e6):

        self.pulser_off()
        # Setup of channel enable, output conditioning as well as trigger setup not shown for simplicity

        #lBytesPerSample = int32(0)

        # Read out used bytes per sample
        #spcm_dwGetParam_i32(self._hCard, SPC_MIINST_BYTESPERSAMPLE, byref(lBytesPerSample));

        # setup the trigger mode -> no trigger
        self._read_out_error()
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE) # (SW trigger, no output)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, 0)

        # Setting up the card mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)  # enable sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, 2)  # Divide on - board mem in two parts
        # spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_STARTSTEP, 0) # Step#0 is the first step after card start

        # Setting up the data memory and transfer data
        llMemSamples = int64(KILO_B(64))
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 0)  # set current configuration switch to segment 0
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, llMemSamples)  # define size of current segment 0

        sample_rate = self.get_sample_rate()
        self._read_out_error()

        a_ch1_signal = (np.sin(np.arange(llMemSamples.value) * (frequency_1/sample_rate) * 2 * np.pi) * (2 ** 15 - 1)).astype(
            dtype=np.int16)
        a_ch2_signal = np.zeros(llMemSamples.value).astype(dtype=np.int16)
        d_ch1_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)
        d_ch2_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)
        d_ch3_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)

        pvBuffer, qwBufferSize = self._create_combined_buffer_data(llMemSamples, {
            'a_ch1': a_ch1_signal, 'a_ch2': a_ch2_signal,
            'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal, 'd_ch3': d_ch3_signal})

        # it is assumed, that the Buffer memory has been allocated and is already filled with valid data
        llMemSamples = int64(KILO_B(64))
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

        # Setting up the data memory and transfer data
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 1) # set current configuration switch to segment 1
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, llMemSamples) # define size of current segment 1

        a_ch1_signal = np.zeros(llMemSamples.value).astype(dtype=np.int16)
        a_ch2_signal = (np.sin(np.arange(llMemSamples.value) * (frequency_2 / sample_rate) * 2 * np.pi) * (2 ** 15 - 1)).astype(
            dtype=np.int16)
        d_ch1_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)
        d_ch2_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)
        d_ch3_signal = np.zeros(llMemSamples.value).astype(dtype=np.bool)

        pvBuffer, qwBufferSize = self._create_combined_buffer_data(llMemSamples, {
            'a_ch1': a_ch1_signal, 'a_ch2': a_ch2_signal,
            'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal, 'd_ch3': d_ch3_signal})
        # it is assumed, that the Buffer memory has been allocated and is already filled with valid data
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

        # Setting up the sequence memory(Only two steps used here as an example)
        lStep = 0 # current step is Step  # 0
        llSegment = 0 # associated with memory
        llLoop = 1 # Pattern will be repeated 10 times
        llNext = 1 # Next step is Step  # 1
        llCondition = SPCSEQ_ENDLOOPONTRIG # Unconditionally leave current step

        # combine all the parameters to one int64 bit value
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        lStep = 1 # current step is Step  # 1
        llSegment = 1 # associated with memory segment 1
        llLoop = 1 # Pattern will be repeated once before condition is checked
        llNext = 0; # Next step is Step  # 0
        llCondition = SPCSEQ_ENDLOOPONTRIG # Repeat current step until a trigger has occurred
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        # Start the card
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

        # ... wait here or do something else ...

        #Stop the card
        # spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP);