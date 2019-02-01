# -*- coding: utf-8 -*-
"""
Use Swabian Instruments PulseStreamer8/2 as a pulse generator.

Protobuf (pb2) and grpc files generated from pulse_streamer.proto
file available at https://www.swabianinstruments.com/static/documentation/PulseStreamer/sections/interface.html#grpc-interface.

Regenerate files for an update proto file using the following:
python3 -m grpc_tools.protoc -I=./ --python_out=. --grpc_python_out=. ./pulse_streamer.proto

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

from core.module import Base, ConfigOption
from core.util.modules import get_home_dir
from interface.pulser_interface import PulserInterface, PulserConstraints
from collections import OrderedDict

import grpc
import os
import hardware.swabian_instruments.pulse_streamer_pb2 as pulse_streamer_pb2
import dill

class PulseStreamer(Base, PulserInterface):
    """Methods to control PulseStreamer.
    """
    _modclass = 'pulserinterface'
    _modtype = 'hardware'

    _pulsestreamer_ip = ConfigOption('pulsestreamer_ip', '192.168.1.100', missing='warn')
    _laser_channel = ConfigOption('laser_channel', 0, missing='warn')
    _uw_x_channel = ConfigOption('uw_x_channel', 2, missing='warn')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self.current_status = -1
        self.sample_rate = 1e9
        self.current_loaded_asset = None

        self._channel = grpc.insecure_channel(self._pulsestreamer_ip + ':50051')

    def on_activate(self):
        """ Establish connection to pulse streamer and tell it to cancel all operations """
        self.pulse_streamer = pulse_streamer_pb2.PulseStreamerStub(self._channel)
        self.pulser_off()
        self.current_status = 0

    def on_deactivate(self):
        del self.pulse_streamer


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

        constraints.sample_rate.min = 1e9
        constraints.sample_rate.max = 1e9
        constraints.sample_rate.step = 0
        constraints.sample_rate.default = 1e9

        constraints.d_ch_low.min = 0.0
        constraints.d_ch_low.max = 0.0
        constraints.d_ch_low.step = 0.0
        constraints.d_ch_low.default = 0.0

        constraints.d_ch_high.min = 3.3
        constraints.d_ch_high.max = 3.3
        constraints.d_ch_high.step = 0.0
        constraints.d_ch_high.default = 3.3

        constraints.waveform_length.min = 1
        constraints.waveform_length.max = 64800000
        constraints.waveform_length.step = 1
        constraints.waveform_length.default = 80

        # sample file length max is not well-defined for PulseStreamer, which collates sequential identical pulses into
        # one. Total number of not-sequentially-identical pulses which can be stored: 1 M.
        constraints.waveform_length.min = 1
        constraints.waveform_length.max = 134217728
        constraints.waveform_length.step = 1
        constraints.waveform_length.default = 1


        """
        # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # channels. Here all possible channel configurations are stated, where only the generic
        # names should be used. The names for the different configurations can be customary chosen.
        activation_conf = OrderedDict()
        activation_conf['yourconf'] = {'a_ch1', 'd_ch1', 'd_ch2', 'a_ch2', 'd_ch3', 'd_ch4'}
        activation_conf['different_conf'] = {'a_ch1', 'd_ch1', 'd_ch2'}
        activation_conf['something_else'] = {'a_ch2', 'd_ch3', 'd_ch4'}
        constraints.activation_config = activation_conf
        """
        activation_config['all'] = ['d_ch1', 'd_ch2', 'd_ch3', 'd_ch4', 'd_ch5', 'd_ch6', 'd_ch7',
                                    'd_ch8']
        constraints.activation_config = activation_config

        return constraints

    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """
        # start the pulse sequence
        self.pulse_streamer.stream(self._sequence)
        self.log.info('Asset uploaded to PulseStreamer')
        self.pulse_streamer.startNow(pulse_streamer_pb2.VoidMessage())
        self.current_status = 1
        return 0

    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        # stop the pulse sequence
        channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        self.pulse_streamer.constant(pulse_streamer_pb2.PulseMessage(ticks=0, digi=channels, ao0=0, ao1=0))
        self.current_status = 0
        return 0

    def load_sequence(self, sequence_name):
        """ Loads a sequence to the channels of the device in order to be ready for playback.
        For devices that have a workspace (i.e. AWG) this will load the sequence from the device
        workspace into the channels.
        For a device without mass memory this will make the waveform/pattern that has been
        previously written with self.write_waveform ready to play.

        @param sequence_name:  dict|list, a dictionary with keys being one of the available channel
                                      index and values being the name of the already written
                                      waveform to load into the channel.
                                      Examples:   {1: rabi_ch1, 2: rabi_ch2} or
                                                  {1: rabi_ch2, 2: rabi_ch1}
                                      If just a list of waveform names if given, the channel
                                      association will be invoked from the channel
                                      suffix '_ch1', '_ch2' etc.

        @return dict: Dictionary containing the actually loaded waveforms per channel.
        """
        self.log.warning('PulseStremer doesn\'t perform streaming!\nload_sequence call ignored.')
        return

    def load_waveform(self, load_dict):
        """ Loads a waveform to the specified channel of the pulsing device.

        @param dict|list load_dict: a dictionary with keys being one of the
                                    available channel index and values being the
                                    name of the already written waveform to load
                                    into the channel. Examples:

                                        {1: rabi_ch1, 2: rabi_ch2}
                                    or
                                        {1: rabi_ch2, 2: rabi_ch1}

                                    If just a list of waveform names if given,
                                    the channel association will be invoked from
                                    the channel suffix '_ch1', '_ch2' etc. A
                                    possible configuration can be e.g.

                                        ['rabi_ch1', 'rabi_ch2', 'rabi_ch3']

        @return dict: Dictionary containing the actually loaded waveforms per
                      channel.

        For devices that have a workspace (i.e. AWG) this will load the waveform
        from the device workspace into the channel. For a device without mass
        memory, this will make the waveform/pattern that has been previously
        written with self.write_waveform ready to play.

        Please note that the channel index used here is not to be confused with the number suffix
        in the generic channel descriptors (i.e. 'd_ch1', 'a_ch1'). The channel index used here is
        highly hardware specific and corresponds to a collection of digital and analog channels
        being associated to a SINGLE wavfeorm asset.
        """
        # Since only one waveform can be present at a time check if only a single name is given
        if isinstance(load_dict, list):
            waveforms = list(set(load_dict))
        elif isinstance(load_dict, dict):
            waveforms = list(set(load_dict.values()))
        else:
            self.log.error('Method load_waveform expects a list of waveform names or a dict.')
            return self.get_loaded_assets()

        if len(waveforms) != 1:
            self.log.error('PulseStreamer pulser expects exactly one waveform name for load_waveform.')
            return self.get_loaded_assets()

        waveform = waveforms[0]
        if waveform != self.__current_waveform_name:
            self.log.error('No waveform by the name "{0}" generated for PulseStreamer pulser.\n'
                           'Only one waveform at a time can be held.'.format(waveform))
            return self.get_loaded_assets()


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
        if analog_samples:
            self.log.error('PulseStreamer pulse generator is purely digital (for now) and does not support waveform '
                           'generation with analog samples.')
            return -1, list()
        if not digital_samples:
            if total_number_of_samples > 0:
                self.log.warning('No samples handed over for waveform generation.')
                return -1, list()
            else:
                self.__current_waveform = bytearray([0])
                self.__current_waveform_name = ''
                return 0, list()


    def reset(self):
        """ Reset the device.

        @return int: error code (0:OK, -1:error)
        """
        channels = self._convert_to_bitmask([self._laser_channel, self._uw_x_channel])
        self.pulse_streamer.constant(pulse_streamer_pb2.PulseMessage(ticks=0, digi=channels, ao0=0, ao1=0))
        self.pulse_streamer.constant(laser_on)
        return 0

    def has_sequence_mode(self):
        """ Asks the pulse generator whether sequence mode exists.

        @return: bool, True for yes, False for no.
        """
        return False


    def get_waveform_names(self):
        """ Retrieve the names of all uploaded waveforms on the device.

        @return list: List of all uploaded waveform name strings in the device workspace.
        """
        return

    def get_sequence_names(self):
        """ Retrieve the names of all uploaded sequence on the device.

        @return list: List of all uploaded sequence name strings in the device workspace.
        """
        return list()

    def delete_waveform(self, waveform_name):
        """ Delete the waveform with name "waveform_name" from the device memory.

        @param str waveform_name: The name of the waveform to be deleted
                                  Optionally a list of waveform names can be passed.

        @return list: a list of deleted waveform names.
        """
        return

    def delete_sequence(self, sequence_name):
        """ Delete the sequence with name "sequence_name" from the device memory.

        @param str sequence_name: The name of the sequence to be deleted
                                  Optionally a list of sequence names can be passed.

        @return list: a list of deleted sequence names.
        """
        return list()

    def get_interleave(self):
        """ Check whether Interleave is ON or OFF in AWG.

        @return bool: True: ON, False: OFF

        Will always return False for pulse generator hardware without interleave.
        """
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
            self.log.error('No interleave functionality available in FPGA pulser.\n'
                           'Interleave state is always False.')
        return False

    def _convert_to_bitmask(self, active_channels):
        """ Convert a list of channels into a bitmask.
        @param numpy.array active_channels: the list of active channels like
                            e.g. [0,4,7]. Note that the channels start from 0.
        @return int: The channel-list is converted into a bitmask (an sequence
                     of 1 and 0). The returned integer corresponds to such a
                     bitmask.
        Note that you can get a binary representation of an integer in python
        if you use the command bin(<integer-value>). All higher unneeded digits
        will be dropped, i.e. 0b00100 is turned into 0b100. Examples are
            bin(0) =    0b0
            bin(1) =    0b1
            bin(8) = 0b1000
        Each bit value (read from right to left) corresponds to the fact that a
        channel is on or off. I.e. if you have
            0b001011
        then it would mean that only channel 0, 1 and 3 are switched to on, the
        others are off.
        Helper method for write_pulse_form.
        """
        bits = 0  # that corresponds to: 0b0
        for channel in active_channels:
            # go through each list element and create the digital word out of
            # 0 and 1 that represents the channel configuration. In order to do
            # that a bitwise shift to the left (<< operator) is performed and
            # the current channel configuration is compared with a bitwise OR
            # to check whether the bit was already set. E.g.:
            #   0b1001 | 0b0110: compare elementwise:
            #           1 | 0 => 1
            #           0 | 1 => 1
            #           0 | 1 => 1
            #           1 | 1 => 1
            #                   => 0b1111
            bits = bits | (1 << channel)
        return bits
