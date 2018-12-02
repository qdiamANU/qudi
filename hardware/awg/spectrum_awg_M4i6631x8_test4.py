# -*- coding: utf-8 -*-

"""
This file contains the Qudi hardware module for the M4i6631x8 series
AWG from Spectrum Instruments.

Written by Andrew Horsley and Sven Bodenstedt, 2018

###############################################################

The M4i6631x8 has two 16-bit analogue channels. Bit from these analogue
channels can be sacrificed to create up to 3 digital output channels. The
assignment of the bits is flexible. We choose:
analog channel 1: bit 0-13: a_ch1
                  bit   15: d_ch1
                  bit   14: d_ch2
analog channel 2: bit 0-14: a_ch2
                  bit   15: d_ch3

###############################################################

Unlike many AWGs, the Spectrum AWGs have no onboard file system, and only
the currently 'loaded' waveforms are actually stored on the device. To integrate with
Qudi operations, we use self.waveform_dict and self.sequence_dict

self.waveform_dict is particularly important, as this is where waveforms written using
self.write_waveform() are stored, and then subsequently recalled in self.load_waveform()
or self.write_sequence()

self.sequence_dict() only includes the name of the currently loaded sequence (as a
dictionary key). The corresponding item is an empty string, ''

###############################################################

# todo: currently only tested for use with both AWG output channels. Will probably break if one channel is deactivated


###############################################################
###############################################################
###############################################################

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
import ctypes
from collections import OrderedDict
from thirdparty.spectrum_instruments.pyspcm import *
import time

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

        self.connected = False
        self.interleave = False  # no interleave mode available for Spectrum AWGs

        # Deactivate all channels at first:
        # todo: this doesn't do anything at the moment
        self.channel_states = {'a_ch1': False, 'a_ch2': False,
                               'd_ch1': False, 'd_ch2': False, 'd_ch3': False}
#

        ## for each analog channel one value
        #self.amplitude_dict = {'a_ch1': 4.0, 'a_ch2': 4.0}
        #self.offset_dict = {}
#
        ## for each digital channel one value
        #self.digital_high_dict = {'d_ch1': 3.3, 'd_ch2': 3.3, 'd_ch3': 3.3}
        #self.digital_low_dict = {'d_ch1': 0.0, 'd_ch2': 0.0, 'd_ch3': 0.0}

        self.waveform_names = list()
        self.waveform_dict = dict()
        self.sequence_names = list()
        self.sequence_dict = dict()

        self.current_loaded_assets = dict()
        self.current_loaded_assets_type = ''

        self.blank_stored_waveform = StoredWaveform()
        #self.use_sequencer = True

        self.current_status = 0    # current_status = 0 means off, not running.

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """
        # open card
        hCard = spcm_hOpen(create_string_buffer(b'/dev/spcm0'))
        self._hCard = hCard
        if self._hCard == None:
            self.log.warning('no card found...')
            return 1

        self.connected = True

        time.sleep(3)

        # prepare the card in a defined state
        self.reset()

        # read out awg card parameters
        self.memory_size_bytes_per_channel = self._spcm_dwGetParam_i64(SPC_PCIMEMSIZE) / self._spcm_dwGetParam_i64(SPC_MIINST_CHPERMODULE)
        self.memory_size_samples_per_channel = self.memory_size_bytes_per_channel / self._spcm_dwGetParam_i64(SPC_MIINST_BYTESPERSAMPLE)
        self.max_sequence_segments = self._spcm_dwGetParam_i64(SPC_SEQMODE_AVAILMAXSEGMENT)  # max number of segments
        self.max_sequence_steps = self._spcm_dwGetParam_i64(SPC_SEQMODE_AVAILMAXSTEPS)  # max number of sequence steps
        self.max_sequence_loops = self._spcm_dwGetParam_i64(SPC_SEQMODE_AVAILMAXLOOP)  # max number of loops for a given step

        # todo: not sure this is doing anything at the moment
        if self.activation_config is None:
            self.activation_config = self.get_constraints().activation_config['config_a12d123']
        elif self.activation_config not in self.get_constraints().activation_config.values():
            self.activation_config = self.get_constraints().activation_config['config_a12d123']

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

        constraints.a_ch_amplitude.min = 0.08 * 2  # factor of 2 because qudi uses peak-to-peak amplitude
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
        if self.current_status == 1:
            self.log.info('Pulser already on!')
            return 0
        elif self.current_status == -1:
            self.log.error('Pulser in error state, cannot switch output on!')
            return -1
        elif self.current_status == 0:
            self.current_status = 1
            self.log.info('Pulser: Switch on the output.')

            # check which AWG channels are required for the current active channels
            channel0 = CHANNEL0 * (self._active_channels['a_ch1'] or self._active_channels['d_ch1'] or
                                   self._active_channels['d_ch2'])
            channel1 = CHANNEL1 * (self._active_channels['a_ch2'] or self._active_channels['d_ch3'])
            #print('pulser_on: channel0={}, channel1={}'.format(channel0, channel1))
            if channel0 == CHANNEL0:
                if self._spcm_dwGetParam_i32(SPC_ENABLEOUT0) == 0:
                    self._spcm_dwSetParam_i32(SPC_ENABLEOUT0, 1)  # enable channel 1 output: a_ch1, d_ch1, d_ch2
            if channel1 == CHANNEL1:
                if self._spcm_dwGetParam_i32(SPC_ENABLEOUT1) == 0:
                    self._spcm_dwSetParam_i32(SPC_ENABLEOUT1, 1)  # enable channel 2 output: a_ch2, d_ch3

            # start AWG card and enable trigger
            if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
                spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

            #self.force_trigger()

        if self.read_out_error():
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
            return -1
        elif self.current_status == 1:
            self.current_status = 0
            self.log.info('Pulser: Switch off the output.')

            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 0)  # disable channel 1 output: a_ch1, d_ch1, d_ch2
            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 0)  # disable channel 2 output: a_ch2, d_ch3

        if self.read_out_error():
            return -1
        else:
            return 0

    def write_waveform(self, name, analog_samples, digital_samples, is_first_chunk, is_last_chunk,
                       total_number_of_samples):
        """
        Write a new waveform to the PC workspace.

        NOTES:
        -- Unlike many AWGs, with an onboard file system, only the active waveform/sequence
        is stored on the Spectrum AWG.
        -- The key output of write_waveform is a StoredWaveform instance, which is stored in
        the dictionary self.waveform_dict[name] = stored_waveform
        -- All sample arrays in analog_samples and digital_samples must be of equal length!
        -- Digital channels are encoded in analog samples for synchronous readout
                -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)

        @param str name: the name of the waveform to be created/append to
        @param dict analog_samples: keys are the generic analog channel names (i.e. 'a_ch1') and
                                    values are 1D numpy arrays of type float32 containing the
                                    voltage samples.
        @param dict digital_samples: keys are the generic digital channel names (i.e. 'd_ch1') and
                                     values are 1D numpy arrays of type bool containing the marker
                                     states.
        @param bool is_first_chunk: [NOT USED FOR SPECTRUM AWG]
                                    Flag indicating if it is the first chunk to write.
                                    If True this method will create a new empty waveform.
                                    If False the samples are appended to the existing waveform.
        @param bool is_last_chunk:  [NOT USED FOR SPECTRUM AWG]
                                    Flag indicating if it is the last chunk to write.
                                    Some devices may need to know when to close the appending wfm.
        @param int total_number_of_samples: The number of sample points for the entire waveform
                                            (not only the currently written chunk)

        @return (int, list): Number of samples written (-1 indicates failed process) and list of
                             created waveform names
        """

        self.pulser_off()

        #print('\nWrite Waveform\n')
        # check if waveform name already exists in stored waveforms. If so, delete previous waveform instance
        # todo: check if previous waveform identical - if so, leave unchanged and return
        if name in self.waveform_names:
            self.delete_waveform(name)
            # del self.waveform_dict[name]
            self.log.info('Overwriting waveform: {}'.format(name))

        # fixme: current waveform name should really only be defined after load_waveform()
        # need to check whether qudi requires it to be defined already though

        ###### Sanity checks and adjustment to data to 16-bit format ##########################

        # perform sanity check on input name, and update stored waveform names
        if isinstance(name, str):
            self.waveform_names.append(name)
            current_waveform_name = [name]
        elif isinstance(name, list):  # this check appears unnecessary, as pulsed GUI sends 'name' as a str
            if len(name) != 1:
                self.log.error('Cannot write waveform: too many waveform names given')
            else:
                self.waveform_names.append(name)
                current_waveform_name = name
        else:
            self.log.error('Cannot write waveform: waveform name is not list or str.')
        #print('\ncurrent_waveform_name = {}, type = {}'.format(current_waveform_name, type(current_waveform_name)))
        #print('self.waveform_names = {}, type = {}'.format(self.waveform_names, type(self.waveform_names)))

        # analogue channel sanity checks and adjustment to 16-bit format:
        for key in analog_samples:
            # perform sanity check on waveform amplitude: should be bounded by +/-1
            if max(analog_samples[key]) > 1 or min(analog_samples[key]) < -1:
                self.log.error('Cannot write waveform: waveform amplitude ({} to {}) for {} is out of bounds (+/-1)'.
                               format(min(analog_samples[key]), max(analog_samples[key]), key))
                return -1, current_waveform_name

            # perform sanity check on number of samples: should match total_number_of_samples
            len_samples = len(analog_samples[key])
            if len_samples != total_number_of_samples:
                self.log.error(
                    'Cannot write waveform: number of samples for {} ({}) does not match expected value ({})'.
                    format(key, len_samples, total_number_of_samples))
                return -1, current_waveform_name

            # adjust analogue amplitude to match expected 16-bit input for Spectrum AWG
            analog_samples[key] = (analog_samples[key]*(2 ** 15 - 1)).astype(dtype=np.int16)

        # digital channel sanity check
        for key in digital_samples:

            # perform sanity check on number of samples: should match total_number_of_samples
            len_samples = len(digital_samples[key])
            if len_samples != total_number_of_samples:
                self.log.error('Cannot write waveform: number of samples for {} ({}) does not match expected value ({})'.
                               format(key, len_samples, total_number_of_samples))
                return -1, current_waveform_name

        number_of_samples_writen = total_number_of_samples

        ############# end sanity checks #####################################################

        # combine analogue and digital sample dictionaries into single dictionary
        channel_data = {**analog_samples, **digital_samples}
        #print(channel_data)

        number_of_samples = int64(total_number_of_samples)
        _, padded_number_of_samples = self._padded_number_of_samples(number_of_samples)

        # create a sample buffer (containing all channels) in a format suitable to upload to the AWG
        pvBuffer, qwBufferSize = self._create_combined_buffer_data(number_of_samples, channel_data)

        # create a stored_waveform object, and add to dictionary of waveforms
        stored_waveform = StoredWaveform()  # HAVE TO create a new instance of class, otherwise we'll overwrite / corrupt previously stored waveforms
        stored_waveform.name = name
        stored_waveform.waveform_buffer = pvBuffer
        stored_waveform.n_samples = padded_number_of_samples.value
        stored_waveform.buffersize = qwBufferSize
        self.waveform_dict[name] = stored_waveform

        if self.read_out_error():
            self.log.error('Error in writing waveform')
            return -1, current_waveform_name
        else:
            return number_of_samples_writen, current_waveform_name

    def write_sequence(self, name, sequence_parameter_list, print_time=True):
        """
        Write a new sequence to the AWG.
        NOTE: unlike many AWGs, with an onboard file system, only the active
        waveform/sequence is stored on the Spectrum AWG. Due to the structure
        of the Spectrum AWG memory, we need to upload the sequence with
        write_sequence(), instead of the standard way, where write_sequence()
        creates a sequence file that is subsequently loaded to the AWG using load_sequence().
        Here, load_sequence() just sends a final initialisation command to the AWG

        @param name: str, the name of the waveform to be created/append to
        @param sequence_parameter_list: list, contains the parameters for each sequence step and
                                        the according waveform names.

        sequence_parameter_list = list( segment_1, segment_2, ...)
        segment_1 = list( ('name',), parameters_dict )
        parameters_dict: 'repetitions': 2               # Number of times to repeat segment
                        'go_to': -1                     # unused
                        'event_jump_to': -1             # Next segment to play after this one
                        'event_trigger': 'ON' / 'OFF'   # unused
                        'wait_for': 'ON' / 'OFF'        # unused
                        'flag_trigger': 'ON' / 'OFF'    # unused
                        'flag_high': 'ON' / 'OFF'       # unused
        todo: understand intended usage of parameters_dict items

        @return: int, number of sequence steps written (-1 indicates failed process)
        """
        start_time = time.time()
        print('write_sequence: name = {}'.format(name))
        #print('sequence_parameter_list = {}, \ntype = {}'.format(sequence_parameter_list, type(sequence_parameter_list)))

        # Check if all waveforms are present in PC memory
        for waveform_tuple, param_dict in sequence_parameter_list:
            for waveform in waveform_tuple:
                if waveform not in self.waveform_dict:
                    self.log.error('Failed to create sequence "{0}" due to waveform "{1}" not '
                                   'present in device memory.'.format(name, waveform))
                    return -1

        # overwrite the sequence dictionary:
        self.sequence_dict = dict()
        self.sequence_dict[name] = ''
        self.sequence_names = list({name})
        self.current_loaded_assets_type = 'sequence'

        # desired minimum number of memory segments for sequence
        n_segments = len(sequence_parameter_list)
        # AWG memory must be divided into 2^n segments
        n_actual_segments = int(pow(2, np.ceil(np.log2(n_segments))))
        # Set up the card in sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)
        # break the memory up into segments
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, n_actual_segments)
        # Set step #0 as the first step replayed after card start
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_STARTSTEP, 0)

        # Loop through and upload the segments:
        i = 0
        total_segments = len(sequence_parameter_list)
        for waveform_tuple, param_dict in sequence_parameter_list:
            for waveform_name in waveform_tuple:

                print('writing segment {}/{}, {}'.format(i+1, total_segments, waveform_name))

                waveform = self.waveform_dict[waveform_name]

                # set memory segment to write to
                spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_WRITESEGMENT, i)
                # define memory size of segment
                spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_SEGMENTSIZE, waveform.n_samples)
                # upload segment data
                spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, waveform.waveform_buffer, 0,
                                       waveform.buffersize)
                spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

                # cycle sample rate, to avoid a bug in AWG that can corrupt the output data
                self._cycle_sample_rate()

                # Define how the segment is replayed in the sequence
                lStep = i  # current step is Step  # i
                llSegment = i  # associated with memory segment #i
                llLoop = param_dict['repetitions']  # Pattern will be repeated this many times
                if i == n_segments - 1:  # todo: this may result in unexpected output for sequences defined non-sequentially
                    llNext = 0
                    if param_dict['wait_for'] == 'ON':
                        llCondition = SPCSEQ_ENDLOOPONTRIG
                        llLoop=1
                    else:
                        #llCondition = SPCSEQ_END  # End of sequence
                        llCondition = SPCSEQ_ENDLOOPALWAYS  # Unconditionally leave current step -> continuously loop sequence
                else:
                    llNext = param_dict['event_jump_to']  # Next step
                    if param_dict['wait_for'] == 'ON':
                        llCondition = SPCSEQ_ENDLOOPONTRIG
                        llLoop = 1
                    else:
                        llCondition = SPCSEQ_ENDLOOPALWAYS  # Unconditionally leave current step
                llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
                spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

                i += 1
        counter = i
        end_time = time.time()

        if self.read_out_error() | counter != len(sequence_parameter_list):
            self.log.error('Failed to upload sequence correctly')
            return -1
        else:
            if print_time:
                self.log.info('finished writing sequence, write time = {} s'.format(end_time-start_time))
            print('finished writing sequence, write time = {} s'.format(end_time-start_time))
            return len(sequence_parameter_list)

    def get_waveform_names(self):
        """ Retrieve the names of all 'uploaded' waveforms in the PC memory.
        NOTE: unlike many AWGs, with an onboard file system, only the active
        waveform/sequence is stored on the Spectrum AWG.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """

        return self.waveform_names

    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """

        return self.sequence_names

    def delete_waveform(self, waveform_del):
        """ Delete the waveform with name "waveform_name" from the PC memory.
        NOTE: unlike many AWGs, with an onboard file system, only the active
        waveform/sequence is stored on the Spectrum AWG.

        @param str waveform_del: The name of the waveform to be deleted
                                  Optionally a list of waveform names can be passed.

        @return list: a list of deleted waveform names.
        """
        #todo: update once waveform name structure is finalised

        # check if input is a string - if so, convert to list
        if isinstance(waveform_del, str):
            waveform_del = list([waveform_del])

        # scan through list of waveforms to delete. Remove from both self.waveform_names and self.waveform_dict
        for ind in range(len(waveform_del)):
            if waveform_del[ind] not in self.waveform_names:
                self.log.warning('[{}] not in waveform_names! Unable to delete.'.format(waveform_del[ind]))
            else:
                while waveform_del[ind] in self.waveform_names:
                    self.waveform_names.remove(waveform_del[ind])
            while waveform_del[ind] in self.waveform_dict:
                del self.waveform_dict[waveform_del[ind]]

        return waveform_del

    def delete_sequence(self, sequence_name):
        """ Delete the sequence with name "sequence_name" from the PC memory.
        NOTE: unlike many AWGs, with an onboard file system, only the active
        waveform/sequence is stored on the Spectrum AWG.

        @param str sequence_name: The name of the sequence to be deleted
                                  Optionally a list of sequence names can be passed.

        @return list: a list of deleted sequence names.
        """
        self.log.warning('Cannot delete sequences from Spectrum AWG (can only overwrite). No action taken.')
        return list(sequence_name)

    def load_waveform(self, load_dict, print_time=True):
        """ Loads a waveform to the specified channel of the pulsing device.
        For devices that have a workspace (i.e. AWG) this will load the waveform from the device
        workspace into the channel.
        For a device without mass memory this will make the waveform/pattern that has been
        previously written with self.write_waveform ready to play.

        # todo: update description of load_dict
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
        #print('\nload_dict = {}, type = {}\n'.format(load_dict, type(load_dict)))

        # If the AWG was previously in sequence mode, need to reset before running in waveform mode
        # todo: find a way around this issue that doesn't require resetting the AWG - probably have to contact manufacturer
        if self.current_loaded_assets_type != 'waveform':
            self.reset_for_waveform_sequence_switch()

        start_time = time.time()

        load_key = ''.join(load_dict)
        waveform = self.waveform_dict[load_key]

        # set up AWG parameters
        llLoops = 0  # 0 -> loop continuously
        self._spcm_dwSetParam_i32(SPC_CARDMODE, SPC_REP_STD_SINGLERESTART)  # The programmed memory is repeated once after each single trigger event.
        self._spcm_dwSetParam_i32(SPC_LOOPS, llLoops)  # number of repetitions
        self._spcm_dwSetParam_i64(SPC_MEMSIZE, waveform.n_samples)  # number of samples per channel
        self.set_software_trigger(software_trigger=True)  # set AWG to software trigger

        # upload the sample buffer to the AWG
        self.log.info("Starting waveform transfer to AWG and waiting until data is in AWG memory")
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, waveform.waveform_buffer, 0, waveform.buffersize)
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        end_time = time.time()
        self.log.info('finished loading waveform, load time = {} s'.format(end_time - start_time))

        self.log.info(
            "Starting the card and waiting for ready interrupt\n(continuous and single restart will have timeout)")
        #spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START)

        self.current_loaded_assets_type = 'waveform'
        self.current_loaded_assets['AWG'] = load_key

        if print_time:
            print('finished loading waveform, load time = {} s'.format(end_time - start_time))

        return  # todo: include return dictionary

    def load_sequence(self, sequence_name):
        """ Loads a sequence to the channels of the device in order to be ready for playback.
        For devices that have a workspace (i.e. AWG) this will load the sequence from the device
        workspace into the channels.

        @param sequence_name:  str, name of the sequence to load

        @return (dict, str): Dictionary with keys being the channel number and values being the
                             respective asset loaded into the channel, string describing the asset
                             type ('waveform' or 'sequence')
        """

        if sequence_name not in self.sequence_dict:
            self.log.error('Failed loading sequence "{}".'.format(sequence_name),
                           ' \n Expected load call for sequence "{}",'.format(next(iter(self.sequence_dict))),
                           'which should be prepared currently on AWG')
            return self.get_loaded_assets()

        self.log.info(
            "Starting the card and waiting for ready interrupt\n(continuous and single restart will have timeout)")
        #spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

        self.current_loaded_assets['AWG'] = sequence_name
        self.current_loaded_assets_type = 'sequence'

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

        print('get_loaded_assets: {}'.format(self.current_loaded_assets))

        asset_type = self.current_loaded_assets_type

        return self.current_loaded_assets, asset_type

    def clear_all(self):
        """ Clears all loaded waveform from the pulse generators RAM.
        @return int: error code (0:OK, -1:error)

        Unused for digital pulse generators without storage capability
        (PulseBlaster, FPGA).

        For the Spectrum AWG, we achieve this by uploading a list of
        zeros to each channel.

        This function is called by the clear button in the GUI
        """
        # TODO check this works, once write_waveform has been implemented
        self.log.warning('clear all does nothing at the moment')

        # upload zeros to all channels
        total_number_of_samples = 32
        analog_samples = {"a_ch1": np.zeros(total_number_of_samples),
                          "a_ch2": np.zeros(total_number_of_samples)}
        digital_samples = {"d_ch1": np.zeros(total_number_of_samples),
                          "d_ch2": np.zeros(total_number_of_samples),
                          "d_ch3": np.zeros(total_number_of_samples)}
        #self.write_waveform('cleared_AWG', analog_samples, digital_samples, True, True, total_number_of_samples)

        self.current_loaded_assets = dict()  # still needed??

        if self.read_out_error():
            return -1
        else:
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

                             100h The next data block as defined in the notify size is available.
                                    It is at least the amount of data available but it also can be more data.
                             200h The data transfer has completed. This status information will only occur
                                    if the notify size is set to zero.
                             400h The data transfer had on overrun (acquisition) or underrun (replay) while doing FIFO transfer.
                             800h An internal error occurred while doing data transfer
        """
        # TODO: this isn't as well-implemented as it could be - would ideally derive active/ready status from AWG query

        status = self._spcm_dwGetParam_i32(SPC_M2STATUS)
        self.log.info('M4i6631x8 status: {0:x} (see manual for more information)'.format(status))
        error = self.read_out_error()

        # change current_status if error detected:
        if error or status >= 800:
            self.current_status = -1

        # otherwise, leave status unchanged (0 or 1)

        # Define the status dictionary for this AWG:
        status_dic = {}
        status_dic[1] = 'Device is active and running.'
        status_dic[0] = 'Device has stopped, but can receive commands.'
        status_dic[-1] = 'Device error'
        # All the other status messages should have higher integer values than 1.

        return self.current_status, status_dic

    def get_sample_rate(self):
        """ Get the sample rate of the pulse generator hardware

        @return float: The current sample rate of the device (in Hz)

        Do not return a saved sample rate in a class variable, but instead
        retrieve the current sample rate directly from the device.
        """
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
        constraint = self.get_constraints().sample_rate
        if sample_rate > constraint.max:
            self.sample_rate = constraint.max
        elif sample_rate < constraint.min:
            self.sample_rate = constraint.min
        else:
            self.sample_rate = sample_rate

        spcm_dwSetParam_i32(self._hCard, SPC_SAMPLERATE, int32(int(sample_rate)))
        self.read_out_error()
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

        self.read_out_error()
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
        self.log.warning('set_analogue_level')
        constraints = self.get_constraints()
        # print('\namplitude constraints: {} to {}'.format(constraints.a_ch_amplitude.min, constraints.a_ch_amplitude.max))
        # print('\namplitude a_ch1 cmd = {}\n'.format(amplitude['a_ch1']))
        if 'a_ch1' in amplitude.keys():
            if (amplitude['a_ch1'] >= constraints.a_ch_amplitude.min) and\
                    (amplitude['a_ch1'] <= constraints.a_ch_amplitude.max):
                # channel0 amplitude in mV
                spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(int(amplitude['a_ch1']*1000/2))) #todo: why the factor of 2?
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

        self.read_out_error()
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
        # todo: no testing/dev done with AWG in any mode other than with all analogue and digital channels activated
        # analogue channels:
        a_channel = self._spcm_dwGetParam_i32(SPC_CHENABLE)

        self._active_channels['a_ch1'] = bool(a_channel & CHANNEL0)
        self._active_channels['a_ch2'] = bool(a_channel & CHANNEL1)

        # digital channels:

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
        # TODO - integrate with activation_config
        # todo: no testing/dev done with AWG in any mode other than with all analogue and digital channels activated
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

        self.read_out_error()
        return self._active_channels

    def get_interleave(self):
        """ Check whether Interleave is ON or OFF in AWG.
        Interleave consists of using two or more AWG channels working at a nominal sample rates to generate a signal
        as if it were created by a higher sample rate device

        @return bool: True: ON, False: OFF

        Unused for pulse generator hardware other than an AWG.
        """
        self.log.info('Interleave mode not available for the Spectrum M4i AWG series!')
        return False

    def set_interleave(self, state=False):
        """ Turns the interleave of an AWG on or off.
        Interleave consists of using two or more AWG channels working at a nominal sample rates to generate a signal
        as if it were created by a higher sample rate device

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
        With the Spectrum AWG, this command is used to set an
        internal register

        @param string command: string containing the command

        @return int: error code (0:OK, -1:error)


        @param lRegister (int): register, that should be set
               plValue   (int): value of the register

        @return int: (-1: error , else parameter value)
        """
        # TODO -> finish this
        # need to work out how to go from str command input to give
        # the required lRegister and plValue

        #lRegister
        #plValue
        #temp = _spcm_dwSetParam_i32(self, lRegister, plValue)
        #return temp

        self.log.warning('AWG write command not yet implemented')
        return -1

    def query(self, question):
        """ Asks the device a 'question' and receive and return an answer from it.

        @param string question: string containing the command. For the Spectrum AWG,
        this is converted to an integer corresponding to the register that should be
        read out. I.e. the input should be an integer-like string.

        @return string: the answer of the device to the 'question' in a string.
        """
        # TODO -> check this works

        temp = _spcm_dwGetParam_i32(self, int(question))
        if temp == -1:
            return 'error'
        else:
            return str(temp)

    def reset(self):
        """ Reset the device, and prepare in a default state:

        All channels activated (a_ch1, a_ch2, d_ch1, d_ch2, d_ch3)
        Output for all channels deactivated
        Analogue amplitudes: 4 V peak-to-peak
        AWG timeout:         25 s

        @return int: error code (0:OK, -1:error)
        """
        #self.__init__()

        self.read_out_error()
        reset_outcome = self._spcm_dwSetParam_i64(SPC_M2CMD, M2CMD_CARD_RESET)
        self.current_status = 0

        # disable output on all channels
        self._spcm_dwSetParam_i64(SPC_ENABLEOUT0, 0)
        self._spcm_dwSetParam_i64(SPC_ENABLEOUT1, 0)

        # set sample rate to max value, 1.25 GS/s
        #self.set_sample_rate(1.25e9)
        spcm_dwSetParam_i32(self._hCard, SPC_SAMPLERATE, int32(int(1.25e9)))
        self.sample_rate = int(1.25e9)

        # set AWG timeout
        awg_timeout = 25000  # timeout in ms
        spcm_dwSetParam_i32(self._hCard, SPC_TIMEOUT, awg_timeout)

        # set active channels
        spcm_dwSetParam_i32(self._hCard, SPC_CHENABLE, CHANNEL0 | CHANNEL1)
        spcm_dwSetParam_i32(self._hCard, SPCM_X0_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT15)
        spcm_dwSetParam_i32(self._hCard, SPCM_X1_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT14)
        spcm_dwSetParam_i32(self._hCard, SPCM_X2_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH1 | SPCM_XMODE_DIGOUTSRC_BIT15)
        self._active_channels = {'a_ch1': True,
                                 'a_ch2': True,
                                 'd_ch1': True,
                                 'd_ch2': True,
                                 'd_ch3': True}

        # set analogue amplitudes
        # input is zero-peak in mV, so input=2000 gives a peak-to-peak amplitude of 4V
        spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(2000))  # 2000 -> set channel 0 amplitude to 4V P2P
        spcm_dwSetParam_i32(self._hCard, SPC_AMP1, int32(2000))  # 2000 -> set channel 1 amplitude to 4V P2P

        error = self.read_out_error()

        if reset_outcome == -1 or error:
            self.log.error('AWG reset error')
            return -1
        else:
            self.log.info('AWG reset successful')
            self.connected = True
            return 0

    def reset_for_waveform_sequence_switch(self):
        """ Reset the device, and prepare in a default state:

        All channels activated (a_ch1, a_ch2, d_ch1, d_ch2, d_ch3)
        Output for all channels deactivated
        Analogue amplitudes: 4 V peak-to-peak
        AWG timeout:         25 s

        @return int: error code (0:OK, -1:error)
        """

        # record current device settings:
        ampl_dict, offset_dict = self.get_analog_level()
        sample_rate = self.get_sample_rate()
        awg_timeout = self._spcm_dwGetParam_i32(SPC_TIMEOUT)

        self.read_out_error()
        reset_outcome = self._spcm_dwSetParam_i64(SPC_M2CMD, M2CMD_CARD_RESET)
        self.current_status = 0

        # return to previous settings
        print(ampl_dict)
        self.set_analog_level(ampl_dict)
        self.set_sample_rate(sample_rate)
        self._spcm_dwSetParam_i32(SPC_TIMEOUT, awg_timeout)

        # disable output on all channels
        self._spcm_dwSetParam_i64(SPC_ENABLEOUT0, 0)
        self._spcm_dwSetParam_i64(SPC_ENABLEOUT1, 0)

        # set active channels
        spcm_dwSetParam_i32(self._hCard, SPC_CHENABLE, CHANNEL0 | CHANNEL1)
        spcm_dwSetParam_i32(self._hCard, SPCM_X0_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT15)
        spcm_dwSetParam_i32(self._hCard, SPCM_X1_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT14)
        spcm_dwSetParam_i32(self._hCard, SPCM_X2_MODE, SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH1 | SPCM_XMODE_DIGOUTSRC_BIT15)
        self._active_channels = {'a_ch1': True,
                                 'a_ch2': True,
                                 'd_ch1': True,
                                 'd_ch2': True,
                                 'd_ch3': True}

        error = self.read_out_error()

        if reset_outcome == -1 or error:
            self.log.error('AWG reset error')
            return -1
        else:
            self.log.info('AWG reset for waveform-sequence switch successful')
            self.connected = True
            return 0

    def has_sequence_mode(self):
        """ Asks the pulse generator whether sequence mode exists.

        @return: bool, True for yes, False for no.
        """
        available_modes = int32(0)
        spcm_dwSetParam_i32(self._hCard, SPC_AVAILCARDMODES, byref(available_modes))
        self.read_out_error()
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
        # not sure if the following check is needed, or why it was written
        #if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
        #    return -1

        if wait:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER | M2CMD_CARD_WAITREADY)
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER)

        if self.read_out_error():
            return -1
        else:
            return 0

    def laser_on(self, laser_channel='d_ch1'):
        """ AWG output results in cw laser, with everything else off"""
        # fixme: laser is off for ca. 25.6 ns at end of every n_samples
        # with n_samples = 32*10000, the waveform length is ca. 256 us, so the 'off' duty cycle is 1e-4

        if self.get_status()[0] > 0:
            self.log.error('Can´t load a waveform, because pulser running. Switch off the pulser and try again.')

            return 0

        else:
            # create waveforms, assuming laser channel is d_ch1
            n_samples = 32 * 10000
            a_ch1_signal = np.zeros(n_samples)
            a_ch2_signal = np.zeros(n_samples)
            if laser_channel == 'd_ch1':
                d_ch1_signal = np.ones(n_samples).astype(dtype=np.bool)
                d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
                d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)
            elif laser_channel == 'd_ch2':
                d_ch1_signal = np.zeros(n_samples).astype(dtype=np.bool)
                d_ch2_signal = np.ones(n_samples).astype(dtype=np.bool)
                d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)
            elif laser_channel == 'd_ch3':
                d_ch1_signal = np.zeros(n_samples).astype(dtype=np.bool)
                d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
                d_ch3_signal = np.ones(n_samples).astype(dtype=np.bool)
            else:
                self.log.error('Invalid laser channel!')
                return -1

            # combine analogue and digital sample dictionaries into two dictionaries
            analog_samples = {'a_ch1': a_ch1_signal,
                              'a_ch2': a_ch2_signal}
            digital_samples = {'d_ch1': d_ch1_signal,
                               'd_ch2': d_ch2_signal,
                               'd_ch3': d_ch3_signal}

            # save waveform in awg.waveform_dict[name]:
            name = 'laser_on'
            self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

            # upload the data to the AWG
            self.load_waveform(name, print_time=False)

            # turn on pulser
            self.pulser_on()

            return

    def laser_off(self):
        """ Turns off laser_on() """

        self.pulser_off()

        return

    def all_off(self):
        """ Writes zeros to all AWG channels and turns AWG output chennels off.
        This is a kill switch: function does not check if any measurement is running """

        if self.current_status == 1:
            self.pulser_off()

        # create waveforms, assuming laser channel is d_ch1
        n_samples = 32 * 1000
        a_ch1_signal = np.zeros(n_samples)
        a_ch2_signal = np.zeros(n_samples)
        d_ch1_signal = np.zeros(n_samples).astype(dtype=np.bool)
        d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
        d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)

        # combine analogue and digital sample dictionaries into two dictionaries
        analog_samples = {'a_ch1': a_ch1_signal,
                          'a_ch2': a_ch2_signal}
        digital_samples = {'d_ch1': d_ch1_signal,
                           'd_ch2': d_ch2_signal,
                           'd_ch3': d_ch3_signal}

        # save waveform in awg.waveform_dict[name]:
        name = 'all_off'
        self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

        # upload the data to the AWG
        self.load_waveform(name, print_time=False)

        # DON'T turn on pulser
        self.pulser_off()

        return

    def write_cw_odmr_waveform(self, freq_start, freq_stop, freq_step, name):
        """
        a_ch1: microwave        (frequency scan, on cw)
        a_ch2: rf               (off)
        d_ch1: green laser      (on cw)
        d_ch2: red laser        (unused)
        d_ch3: counter trigger  (pulse at start of each frequency step)

         """

        self._frequency_list = list(np.arange(freq_start, freq_stop, freq_step))
        n_freq_steps = len(self._frequency_list)

        n_samples_per_step = 400 * 32  # with a 1.25 GS/s sampling rate, this equates to 10.2us per frequency step
        n_samples = n_samples_per_step * n_freq_steps

        sample_rate = self.get_sample_rate()
        trigger_length_seconds = 20e-9  # trigger length in seconds
        trigger_length_samples = int(trigger_length_seconds * sample_rate)

        a_ch1_signal = np.array([])
        a_ch2_signal = np.array([])
        d_ch1_signal = np.array([])
        d_ch2_signal = np.array([])
        d_ch3_signal = np.array([])

        ##### Generate Data and Write Waveforms ########################

        digital_on = np.ones(n_samples_per_step).astype(dtype=np.bool)
        digital_off = np.zeros(n_samples_per_step).astype(dtype=np.bool)
        trigger_signal = digital_off
        trigger_signal[0:trigger_length_samples] = True

        for i in range(n_freq_steps):
            # print('\ni={}'.format(i))

            mw_signal = np.sin(np.arange(n_samples_per_step) * (self._frequency_list[i] / sample_rate) * 2 * np.pi)

            a_ch1_signal = np.concatenate((a_ch1_signal, mw_signal))
            a_ch2_signal = np.concatenate((a_ch2_signal, np.zeros(n_samples_per_step)))
            d_ch1_signal = np.concatenate((d_ch1_signal, digital_on))
            d_ch2_signal = np.concatenate((d_ch2_signal, digital_off))
            d_ch3_signal = np.concatenate((d_ch3_signal, trigger_signal))

        # combine analogue and digital sample dictionaries into single dictionary
        analog_samples = {'a_ch1': a_ch1_signal,
                          'a_ch2': a_ch2_signal}
        digital_samples = {'d_ch1': d_ch1_signal,
                           'd_ch2': d_ch2_signal,
                           'd_ch3': d_ch3_signal}

        # save waveform in awg.waveform_dict[name]:
        self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

        return

    def write_triggered_cw_odmr_list_sequence(self, frequency_list, sequence_name = 'odmr_cw_list'):
        """
        a_ch1: microwave        (frequency scan, on cw)
        a_ch2: rf               (off)
        d_ch1: green laser      (on cw)
        d_ch2: red laser        (unused)
        d_ch3: counter trigger  (pulse at start of each frequency step)

         """

        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

        self._frequency_list = frequency_list
        n_freq_steps = len(self._frequency_list)

        #### AWG setup #############

        # desired minimum number of memory segments for sequence
        # the '+1' is for a blank segment at the end of the sequence, which we add to match the odmr logic architecture
        n_segments = n_freq_steps + 1
        # AWG memory must be divided into 2^n segments
        n_actual_segments = int(pow(2, np.ceil(np.log2(n_segments))))

        # desired number of samples per segment
        desired_n_samples = 400 * 32  # should be a multiple of 32
        memory_size_bytes_per_channel = self._spcm_dwGetParam_i64(SPC_PCIMEMSIZE) / self._spcm_dwGetParam_i64(SPC_MIINST_CHPERMODULE)
        memory_size_samples_per_channel = memory_size_bytes_per_channel / self._spcm_dwGetParam_i64(SPC_MIINST_BYTESPERSAMPLE)
        max_n_samples = memory_size_samples_per_channel / n_actual_segments
        if desired_n_samples > max_n_samples:
            print('Desired_n_samples > max_n_samples. \nSetting memory allocation for each segment to max_n_samples')
            n_samples = max_n_samples
        else:
            n_samples = desired_n_samples
        llMemSamples = int64(n_samples)

        ##### Generate Data and Write Waveforms ########################
        sample_rate = self.get_sample_rate()
        for i in range(n_freq_steps):
            # print('\ni={}'.format(i))

            mw_signal = np.sin(np.arange(n_samples) * (self._frequency_list[i] / sample_rate) * 2 * np.pi)

            a_ch1_signal = mw_signal
            a_ch2_signal = np.zeros(n_samples)
            d_ch1_signal = np.ones(n_samples).astype(dtype=np.bool)
            d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
            d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)

            # combine analogue and digital sample dictionaries into single dictionary
            analog_samples = {'a_ch1': a_ch1_signal,
                              'a_ch2': a_ch2_signal}
            digital_samples = {'d_ch1': d_ch1_signal,
                               'd_ch2': d_ch2_signal,
                               'd_ch3': d_ch3_signal}

            # save waveform in awg.waveform_dict[name]:
            name = 'odmr_freqstep{}'.format(i)
            self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

        # write blank waveform for end of sequence:
        a_ch1_signal = np.zeros(n_samples)
        a_ch2_signal = np.zeros(n_samples)
        # d_ch1_signal = np.ones(n_samples).astype(dtype=np.bool)  # leaving the laser on
        d_ch1_signal = np.zeros(n_samples).astype(dtype=np.bool)
        d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
        d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)

        # combine analogue and digital sample dictionaries into single dictionary
        analog_samples = {'a_ch1': a_ch1_signal,
                          'a_ch2': a_ch2_signal}
        digital_samples = {'d_ch1': d_ch1_signal,
                           'd_ch2': d_ch2_signal,
                           'd_ch3': d_ch3_signal}

        # save waveform in awg.waveform_dict[name]:
        name = 'odmr_freqstep{}'.format(n_freq_steps)
        self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

        ##### Write Sequence ########################

        sequence_parameter_list = list()
        for i in range(n_freq_steps):
            # print('\ni={}'.format(i))

            name = 'odmr_freqstep{}'.format(i)
            parameters_dict = dict()
            parameters_dict['repetitions'] = 1
            parameters_dict['go_to'] = 0
            parameters_dict['event_jump_to'] = i + 1
            parameters_dict['event_trigger'] = 'OFF'
            parameters_dict['wait_for'] = 'ON'
            parameters_dict['flag_trigger'] = 'OFF'
            parameters_dict['flag_high'] = 'OFF'

            segment = [(name,), parameters_dict]
            sequence_parameter_list.append(segment)

        # append blank segment at end of sequence, to match with odmr logic architecture (e.g. extra trigger at end of sequence)
        name = 'odmr_freqstep{}'.format(n_freq_steps)
        parameters_dict = dict()
        parameters_dict['repetitions'] = 1
        parameters_dict['go_to'] = 0
        parameters_dict['event_jump_to'] = 0
        parameters_dict['event_trigger'] = 'OFF'
        parameters_dict['wait_for'] = 'ON'
        parameters_dict['flag_trigger'] = 'OFF'
        parameters_dict['flag_high'] = 'OFF'

        segment = [(name,), parameters_dict]
        sequence_parameter_list.append(segment)

        self.write_sequence(sequence_name, sequence_parameter_list)

        self.set_trigger('ext0')

        return

    def write_triggered_cw_odmr_sweep_sequence(self, freq_start, freq_stop, freq_step, sequence_name='odmr_cw_sweep'):
        """
        a_ch1: microwave        (frequency scan, on cw)
        a_ch2: rf               (off)
        d_ch1: green laser      (on cw)
        d_ch2: red laser        (unused)
        d_ch3: counter trigger  (pulse at start of each frequency step)

         """
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

        self._frequency_list = list(np.arange(freq_start, freq_stop, freq_step))
        n_freq_steps = len(self._frequency_list)

        #### AWG setup #############

        # desired minimum number of memory segments for sequence
        n_segments = n_freq_steps
        # AWG memory must be divided into 2^n segments
        n_actual_segments = int(pow(2, np.ceil(np.log2(n_segments))))

        # desired number of samples per segment
        desired_n_samples = 400 * 32  # should be a multiple of 32
        memory_size_bytes_per_channel = self._spcm_dwGetParam_i64(SPC_PCIMEMSIZE) / self._spcm_dwGetParam_i64(SPC_MIINST_CHPERMODULE)
        memory_size_samples_per_channel = memory_size_bytes_per_channel / self._spcm_dwGetParam_i64(SPC_MIINST_BYTESPERSAMPLE)
        max_n_samples = memory_size_samples_per_channel / n_actual_segments
        if desired_n_samples > max_n_samples:
            print('Desired_n_samples > max_n_samples. \nSetting memory allocation for each segment to max_n_samples')
            n_samples = max_n_samples
        else:
            n_samples = desired_n_samples
        llMemSamples = int64(n_samples)

        ##### Generate Data and Write Waveforms ########################
        sample_rate = self.get_sample_rate()
        for i in range(n_freq_steps):
            # print('\ni={}'.format(i))

            mw_signal = np.sin(np.arange(n_samples) * (self._frequency_list[i] / sample_rate) * 2 * np.pi)

            a_ch1_signal = mw_signal
            a_ch2_signal = np.zeros(n_samples)
            d_ch1_signal = np.ones(n_samples).astype(dtype=np.bool)
            d_ch2_signal = np.zeros(n_samples).astype(dtype=np.bool)
            d_ch3_signal = np.zeros(n_samples).astype(dtype=np.bool)

            # combine analogue and digital sample dictionaries into single dictionary
            analog_samples = {'a_ch1': a_ch1_signal,
                              'a_ch2': a_ch2_signal}
            digital_samples = {'d_ch1': d_ch1_signal,
                               'd_ch2': d_ch2_signal,
                               'd_ch3': d_ch3_signal}

            # save waveform in awg.waveform_dict[name]:
            name = 'odmr_freqstep{}'.format(i)
            self.write_waveform(name, analog_samples, digital_samples, True, True, n_samples)

        ##### Write Sequence ########################

        sequence_parameter_list = list()
        for i in range(n_freq_steps):
            # print('\ni={}'.format(i))

            name = 'odmr_freqstep{}'.format(i)
            parameters_dict = dict()
            parameters_dict['repetitions'] = 1
            parameters_dict['go_to'] = 0
            parameters_dict['event_jump_to'] = i + 1
            parameters_dict['event_trigger'] = 'OFF'
            parameters_dict['wait_for'] = 'ON'
            parameters_dict['flag_trigger'] = 'OFF'
            parameters_dict['flag_high'] = 'OFF'

            segment = [(name,), parameters_dict]
            sequence_parameter_list.append(segment)

        self.write_sequence(sequence_name, sequence_parameter_list)

        self.set_trigger('ext0')

        return

    def read_out_error(self):
        """checks the error state and prints it out. Errors must be read out before the AWG
        can accept further commands

        Return value is the error code:
            0 -> no error
            >0 -> error, check AWG manual for code meaning

        """

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
        # TODO -> shift all usage to query??
        value = int32(0)
        spcm_dwGetParam_i32(self._hCard, lRegister, byref(value))
        #if self.read_out_error():
        #   return -1
        outcome = self.read_out_error()
        if outcome:
            self.log.error('_spcm_dwGetParam_i32 error = {}'.format(outcome))
            return -1
        else:
            return value.value

    def _spcm_dwGetParam_i64(self, lRegister):
        """
        Reads an internal register.

        @param lRegister (int): register, that should be read out

        @return int: (-1: error , else parameter value)
        """
        # TODO -> shift all usage to query??
        value = int64(0)
        spcm_dwGetParam_i64(self._hCard, lRegister, byref(value))
        if self.read_out_error():
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
        # TODO -> shift all usage to write??
        spcm_dwSetParam_i32(self._hCard, lRegister, int32(plValue))
        if self.read_out_error():
           return -1
        else:
            return 0

    def _spcm_dwSetParam_i64(self, lRegister, plValue):
        """
        Sets an internal register.

        @param lRegister (int): register, that should be set
               plValue   (int): value of the register

        @return int: (-1: error , else parameter value)
        """
        # TODO -> shift all usage to write??
        spcm_dwSetParam_i64(self._hCard, lRegister, int64(plValue))
        if self.read_out_error():
           return -1
        else:
            return 0

    def _cycle_sample_rate(self):
        sample_rate = self.get_sample_rate()
        self.set_sample_rate(600e6)
        self.set_sample_rate(sample_rate)
        self.log.info('AWG sample rate cycled {} MS/s -> {} MS/s -> {} MS/s to fix "corrupted-output-bug"'.format(
            int(sample_rate / 1e6), int(600), int(sample_rate / 1e6)))

    def _create_combined_buffer_data(self, number_of_samples, channel_data):
        """
        digital channels are encoded in analog samples for synchronous readout
        -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)

        @param number_of_samples (int64): number of samples (must be identical for all channels)
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

        # sanity check:
        for key in channel_data.keys():
            if channel_data[key].size != number_of_samples.value:
                self.log.error('Channel data for channel {} has the wrong size!'.format(key))
                return -1

        # create buffer of zeros for active channels that have not been provided with an output signal
        for key in active_channels.keys():
            if key not in channel_data.keys():
                #print('Creating zero buffer for {}').format(key)
                channel_data[key] = np.zeros(number_of_samples.value, dtype=np.dtype(np.int16))

        # buffer length must be an integer multiple of 32. Pad input data to round to next integer multiple:
        padding, padded_number_of_samples = self._padded_number_of_samples(number_of_samples)
        if padding != 0:
            self.log.info('padding samples with {} zeros to increase sample number to integer multiple of [32 samples]'.format(padding))
        #print('padded_number_of_samples = {}, type {}'.format(padded_number_of_samples, type(padded_number_of_samples)))
            for key in active_channels.keys():
                channel_data[key] = np.append(channel_data[key], np.zeros(padding))
            #print('{}, new size = {}, padding = {}, type = {}'.format(key, channel_data[key].size, padding, type(padding)))

        if active_channels['d_ch2']:
           #print('Merging d_ch1 + d_ch2 with a_ch1')
            awg_ch1_data = ((channel_data['a_ch1'].astype(np.uint16) >> 2) |
                            (channel_data['d_ch1'].astype(bool) << 15) |
                            (channel_data['d_ch2'].astype(bool) << 14))
        elif active_channels['d_ch1']:
            #print('Merging d_ch1 with a_ch1')
            awg_ch1_data = ((channel_data['a_ch1'].astype(np.uint16) >> 1) |
                            (channel_data['d_ch1'].astype(bool) << 15))
        else:
            awg_ch1_data = channel_data['a_ch1'].astype(np.uint16)

        if active_channels['d_ch3']:
            #print('Merging d_ch3 with a_ch2')
            awg_ch2_data = ((channel_data['a_ch2'].astype(np.uint16) >> 1) |
                            (channel_data['d_ch3'].astype(np.uint16) << 15))
        else:
            awg_ch2_data = channel_data['a_ch2'].astype(np.uint16)

        # Combine the two channels into a single data transfer buffer.
        # Samples for the two channels are ordered according to: A0 B0 A1 B1 A2 B2 .... An Bn,
        # where An is the nth sample for channel 1, and Bn the nth sample for channel 2
        if active_channels['a_ch1'] and active_channels['a_ch2']:
            combined_channel_data = np.vstack((awg_ch1_data, awg_ch2_data)).ravel('F').astype(np.int16)
        elif active_channels['a_ch1']:
            combined_channel_data = awg_ch1_data.astype(np.int16)
        elif active_channels['a_ch2']:
            combined_channel_data = awg_ch2_data.astype(np.int16)
        else:
            self.log.error('No channel activated.')
            return -1

        # set up software buffer
        chcount = self._spcm_dwGetParam_i32(SPC_CHCOUNT)
        lBytesPerSample = self._spcm_dwGetParam_i32(SPC_MIINST_BYTESPERSAMPLE)

        #print('padded_number_of_samples.value={}, chcount={}, lBytesPerSample={}'.format(padded_number_of_samples.value,chcount,lBytesPerSample))
        #print('qwBuffer product = {}'.format(padded_number_of_samples.value * lBytesPerSample * chcount))
        qwBufferSize = uint64(padded_number_of_samples.value * lBytesPerSample * chcount)
        #print('qwBufferSize = {}, type {}'.format(qwBufferSize.value, type(qwBufferSize)))
        pvBuffer = create_string_buffer(qwBufferSize.value)

        # calculate the data
        pnBuffer = cast(pvBuffer, ptr16)
        # TODO: replace for loop with single line, e.g.: pnBuffer = combined_channel_data
        # or does it even have to be a pointer at all??
        for i in range(0, padded_number_of_samples.value * chcount, 1):
            pnBuffer[i] = combined_channel_data[i]

        return pvBuffer, qwBufferSize

    def _padded_number_of_samples(self, number_of_samples):
        """
        Spectrum AWG can only accept sample uploads with a length that is an
        integer multiple of [32 samples] long. Need to pad sample buffer to fulfill this requirement

        For sequence mode, there is also a minimum sample length for each segment:
        Using only 1 channel: minimum segment size = 384 samples
        Using both channels: minimum segment size = 192 samples

        For simplicity, we apply this minimum segment size of 384 samples to all waveforms.

        @param number_of_samples: desired number of samples input to the AWG
        @return padded_number_of_samples: input number rounded up to nearest multiple of [32 samples]
        @return padding: number of samples (e.g. zeros) required to pad the input data
        """

        min_awg_memory_stepsize = 32
        min_sample_length = 384
        padding = int((np.ceil(number_of_samples.value / min_awg_memory_stepsize) - number_of_samples.value / min_awg_memory_stepsize) * min_awg_memory_stepsize)
        padded_number_of_samples = int64(number_of_samples.value + padding)

        if padded_number_of_samples.value < min_sample_length:
            padded_number_of_samples = min_sample_length
            padding = padded_number_of_samples.value - number_of_samples.value

        return padding, padded_number_of_samples

    def prepare_in_sequence_mode(self, n_segments):
        """ Prepares AWG in sequence mode and divides memory into segments

        @param int n_segments: number of desired memory segments (e.g. number of ODMR frequency steps)

        @return int n_actual_segments: number of actual memory segments, which must be a power of 2
        """

        # AWG memory must be divided into 2^n segments
        n_actual_segments = int(pow(2, np.ceil(np.log2(n_segments))))

        # todo: check against AWG limits

        # Set up the card in sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)

        # break the memory up into segments
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, n_actual_segments)

        if self.read_out_error():
            self.log.error('prepare_in_sequence_mode error')
            return -1
        else:
            return n_actual_segments

    def set_software_trigger(self, software_trigger = False):
        """
        Sets up the trigger setting for the awg.

        @param software_trigger bool: enables (True) or disables (False) the software trigger
        If disabled, the AWG is set to have no trigger (and can only be triggered by self.force_trigger)

        @return int: error code (0:OK, -1:error)
        """

        # todo: functions to change between software and external triggers
        # todo: check that regular operations aren't affected by the trigger changes

        # Send stop command: prevents crash when changing trigger whilst card running, ignored if already stopped
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, SPC_TMASK_NONE)



        # Setting up the trigger
        if software_trigger:
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)

        if self.read_out_error():
            return -1
        else:
            return 0

    def set_trigger(self, trigger_src, trigger_level=2000):
        """
        Sets up the trigger setting for the awg.

        @param trigger_src str: software, external triggers: ext0 or ext1
        @param trigger_level int (optional): trigger level for external triggers, in mV

        @return int: error code (0:OK, -1:error)
        """
        # todo: AWG has much richer trigger functionality available than is currently implemented here

        # Send stop command: prevents crash when changing trigger whilst card running, ignored if already stopped
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)

        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, SPC_TMASK_NONE)  # define the trigger AND mask to be empty
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_DELAY, 0)  # delay units are [sample clocks]

        if trigger_src == 'ext0':
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_EXT0_LEVEL0, trigger_level)  # set trigger level, in mV
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_EXT0_MODE, SPC_TM_POS)  # set to trigger off positive edge
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_EXT0)  # define the trigger OR mask to only include ext0

        if trigger_src == 'ext1':
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_EXT1_LEVEL0, trigger_level)  # set trigger level, in mV
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_EXT1_MODE, SPC_TM_POS)  # set to trigger off positive edge
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_EXT1)  # define the trigger OR mask to only include ext1
            
        if trigger_src == 'software':
            spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)  # define the trigger OR mask to only include the software trigger

        if self.read_out_error():
            return -1
        else:
            return 0

    def _set_digital_async_on(self,channel):
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


        if self.read_out_error():
            return -1
        else:
            return 0

    def _set_digital_async_off(self, channel):
        """
        Turns off one of the Multi Purpose I/O Lines as a digital output

        @param channel int: select channel 0,1,2

        @return int: error code (0:OK, -1:error)
        """



        if self.read_out_error():
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
        self.read_out_error()
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
        self.read_out_error()

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
        llNext = 0 # Next step is Step  # 0
        llCondition = SPCSEQ_ENDLOOPONTRIG # Repeat current step until a trigger has occurred
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        # Start the card
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)
        print('start')
        # ... wait here or do something else ...

        #Stop the card
        # spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP);

    def _load_sequence_hardcoded(self, channel_data):
        """ Loads a sequence to the awg

        @param channel_data: list of dictionaries with channel_data:
                             [{'a_ch1' : np.array([...]), 'd_ch1' : np.array([])},
                              {'a_ch1' : np.array([...]), 'd_ch1' : np.array([])}
                             ]

        @return int: error code (0:OK, -1:error)
        """

        if len(channel_data) <= 1:
            self.log.error('Length of frequency list must be larger than 1')
            return -1

        # just zeros
        channel_data.append({'d_ch1': np.full(384, 1)})
        channel_data.append({'d_ch1': np.full(384, 0)})
        max_segments = len(channel_data)
        max_steps = 2 * max_segments
        if max_steps > self.get_constraints().sequence_num.max:
            self.log.error('Number of sequence steps ({}) is larger than the hardware constraints ({})'.format(
                max_steps, self.get_constraints().sequence_num.max))
            return -1

        # set trigger settings
        self.set_software_trigger()

        # Setting up the card mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)  # enable sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, max_segments)  # Divide on - board mem in parts
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_STARTSTEP, 0)  # Step#0 is the first step after card start

        # Setting up the data memory and transfer data
        for i in range(max_segments):
            llMemSamples = int64(list(channel_data[i].values())[0].size)

            # set current configuration switch to segment i
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, i)
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE,
                                llMemSamples)  # define size of current segment 0

            pvBuffer, qwBufferSize = self._create_combined_buffer_data(llMemSamples, channel_data[i])

            # it is assumed, that the Buffer memory has been allocated and is already filled with valid data
            spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

            print('\rUploaded segment {:4d} of {:4d}. '.format(i + 1, max_segments) +
                  '|' * round(50 * (i + 1) / max_segments) + '-' * round(50 - 50 * (i + 1) / max_segments), end='')
        print('')

    def _test_waveform_upload(self, frequency_1=100e6, frequency_2=200e6, number_of_samples=int64(KILO_B(64))):


        sample_rate = self.get_sample_rate()

        a_ch1_signal = (np.sin(np.arange(number_of_samples.value) * (frequency_1 / sample_rate) * 2 * np.pi) *
                        (2 ** 15 - 1)).astype(dtype=np.int16)

        a_ch2_signal = (np.sin(np.arange(number_of_samples.value) * (frequency_2 / sample_rate) * 2 * np.pi) *
                        (2 ** 15 - 1)).astype(dtype=np.int16)

        d_ch1_signal = np.zeros(number_of_samples.value).astype(dtype=np.bool)
        d_ch2_signal = np.zeros(number_of_samples.value).astype(dtype=np.bool)
        d_ch3_signal = np.zeros(number_of_samples.value).astype(dtype=np.bool)

        d_ch1_signal[0:10] = True

        analog_samples = {'a_ch1': a_ch1_signal,
                          'a_ch2': a_ch2_signal}

        digital_samples = {'d_ch1': d_ch1_signal,
                           'd_ch2': d_ch2_signal,
                           'd_ch3': d_ch3_signal}

        self.write_waveform('name', analog_samples, digital_samples, True, True, number_of_samples.value)

        # cycle sample rate to avoid AWG bug that corrupts the output signals
        self._cycle_sample_rate()

        self.pulser_on()



class StoredWaveform:

    def __init__(self):
        self.name = ''
        self.waveform_buffer = ()
        self.n_samples = ()
        self.buffersize = ()

class StoredSequence:

    def __init__(self):
        self.name = ''
        self.waveform_buffer = ()
        self.n_samples = ()
        self.buffersize = ()
