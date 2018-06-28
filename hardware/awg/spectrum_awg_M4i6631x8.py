# -*- coding: utf-8 -*-

"""
This file contains the Qudi hardware module for AWG70000 Series.

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


from core.util.modules import get_home_dir
import os
import sys
import time
import re
import visa
import numpy as np
import ctypes
import matplotlib.pyplot as plt
from socket import socket, AF_INET, SOCK_STREAM
from ftplib import FTP
from collections import OrderedDict
from fnmatch import fnmatch
from thirdparty.spectrum_instruments.pyspcm import *

from core.module import Base, ConfigOption
from interface.pulser_interface import PulserInterface, PulserConstraints

class AWGSpectrumM4i6631x8(Base, PulserInterface):
    """

    """
    _modclass = 'awgSpectrumM4i6631x8'
    _modtype = 'hardware'

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """
        config = self._configuration

        # open card
        hCard = spcm_hOpen(create_string_buffer(b'/dev/spcm0'))
        self._hCard = hCard
        if hCard == None:
            self.log.warning('no card found...')
            return 1

        #prepare the card in a defined state
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_RESET)

        # read type, function and sn and check for A/D card
        lCardType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCITYP, byref(lCardType))
        lSerialNumber = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCISERIALNO, byref(lSerialNumber))
        lFncType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_FNCTYPE, byref(lFncType))

        sCardName = 'M4i.6631-x8'

        spcm_dwSetParam_i32(hCard, SPC_TIMEOUT, 15000)  # timeout 15 s

        self._active_channels = {'a_ch1' : True,
                                 'a_ch2' : True,
                                 'd_ch1' : False,
                                 'd_ch2' : False,
                                 'd_ch3' : False,
                                }
        self.set_active_channels(self._active_channels)

    def on_deactivate(self):
        """ Required tasks to be performed during deactivation of the module.
        """
        if self._hCard != None:
            spcm_vClose(self._hCard)
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
            {<descriptor_str>: <channel_list>,
             <descriptor_str>: <channel_list>,
             ...}

        If the constraints cannot be set in the pulsing hardware (e.g. because it might have no
        sequence mode) just leave it out so that the default is used (only zeros).

        # Example for configuration with default values:
        constraints = PulserConstraints()

        # The file formats are hardware specific.
        constraints.waveform_format = ['wfm', 'wfmx']
        constraints.sequence_format = ['seq', 'seqx']

        constraints.sample_rate.min = 10.0e6
        constraints.sample_rate.max = 12.0e9
        constraints.sample_rate.step = 10.0e6
        constraints.sample_rate.default = 12.0e9

        constraints.a_ch_amplitude.min = 0.02
        constraints.a_ch_amplitude.max = 2.0
        constraints.a_ch_amplitude.step = 0.001
        constraints.a_ch_amplitude.default = 2.0

        constraints.a_ch_offset.min = -1.0
        constraints.a_ch_offset.max = 1.0
        constraints.a_ch_offset.step = 0.001
        constraints.a_ch_offset.default = 0.0

        constraints.d_ch_low.min = -1.0
        constraints.d_ch_low.max = 4.0
        constraints.d_ch_low.step = 0.01
        constraints.d_ch_low.default = 0.0

        constraints.d_ch_high.min = 0.0
        constraints.d_ch_high.max = 5.0
        constraints.d_ch_high.step = 0.01
        constraints.d_ch_high.default = 5.0

        constraints.sampled_file_length.min = 80
        constraints.sampled_file_length.max = 64800000
        constraints.sampled_file_length.step = 1
        constraints.sampled_file_length.default = 80

        constraints.waveform_num.min = 1
        constraints.waveform_num.max = 32000
        constraints.waveform_num.step = 1
        constraints.waveform_num.default = 1

        constraints.sequence_num.min = 1
        constraints.sequence_num.max = 8000
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

        constraints.trigger_in.min = 0
        constraints.trigger_in.max = 2
        constraints.trigger_in.step = 1
        constraints.trigger_in.default = 0

        constraints.event_jump_to.min = 0
        constraints.event_jump_to.max = 8000
        constraints.event_jump_to.step = 1
        constraints.event_jump_to.default = 0

        constraints.go_to.min = 0
        constraints.go_to.max = 8000
        constraints.go_to.step = 1
        constraints.go_to.default = 0

        # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # channels. Here all possible channel configurations are stated, where only the generic
        # names should be used. The names for the different configurations can be customary chosen.
        activation_conf = OrderedDict()
        activation_conf['yourconf'] = ['a_ch1', 'd_ch1', 'd_ch2', 'a_ch2', 'd_ch3', 'd_ch4']
        activation_conf['different_conf'] = ['a_ch1', 'd_ch1', 'd_ch2']
        activation_conf['something_else'] = ['a_ch2', 'd_ch3', 'd_ch4']
        constraints.activation_config = activation_conf
        """
        # Example for configuration with default values:
        constraints = PulserConstraints()

        # The file formats are hardware specific.
        constraints.waveform_format = ['wfm']
        constraints.sequence_format = ['seq']

        constraints.sample_rate.min = 50e6
        constraints.sample_rate.max = 1.25e9
        constraints.sample_rate.step = 1
        constraints.sample_rate.default = 1.25e9

        constraints.a_ch_amplitude.min = 0.08*2
        constraints.a_ch_amplitude.max = 2.0*2
        constraints.a_ch_amplitude.step = 0.001*2
        constraints.a_ch_amplitude.default = 2.0*2

        #offset is fixed to 0.0V for Spectrum M4i series
        constraints.a_ch_offset.min = 0.0
        constraints.a_ch_offset.max = 0.0
        constraints.a_ch_offset.step = 0.000
        constraints.a_ch_offset.default = 0.0

        constraints.d_ch_low.min = 0.0
        constraints.d_ch_low.max = 0.0
        constraints.d_ch_low.step = 0.00
        constraints.d_ch_low.default = 0.0

        constraints.d_ch_high.min = 3.3
        constraints.d_ch_high.max = 3.3
        constraints.d_ch_high.step = 0.00
        constraints.d_ch_high.default = 3.3

        constraints.sampled_file_length.min = 80
        constraints.sampled_file_length.max = 64800000
        constraints.sampled_file_length.step = 1
        constraints.sampled_file_length.default = 80

        constraints.waveform_num.min = 1
        constraints.waveform_num.max = 32000
        constraints.waveform_num.step = 1
        constraints.waveform_num.default = 1

        constraints.sequence_num.min = 1
        constraints.sequence_num.max = 8000
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

        constraints.trigger_in.min = 0
        constraints.trigger_in.max = 2
        constraints.trigger_in.step = 1
        constraints.trigger_in.default = 0

        constraints.event_jump_to.min = 0
        constraints.event_jump_to.max = 8000
        constraints.event_jump_to.step = 1
        constraints.event_jump_to.default = 0

        constraints.go_to.min = 0
        constraints.go_to.max = 8000
        constraints.go_to.step = 1
        constraints.go_to.default = 0

        # the name a_ch<num> and d_ch<num> are generic names, which describe UNAMBIGUOUSLY the
        # channels. Here all possible channel configurations are stated, where only the generic
        # names should be used. The names for the different configurations can be customary chosen.
        activation_conf = OrderedDict()
        activation_conf['all'] = ['a_ch1', 'a_ch2', 'd_ch1', 'd_ch2', 'd_ch3']
        activation_conf['digital_1'] = ['a_ch1', 'a_ch2', 'd_ch1']
        activation_conf['digital_2'] = ['a_ch1', 'a_ch2', 'd_ch1', 'd_ch3']
        constraints.activation_config = activation_conf

        return constraints

    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """
        spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 1)  # enable analog output 1
        spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 1)  # enable analog output 1
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_FORCETRIGGER)
        err = self._readOutError()
        if err:
            return -1
        else:
            return 0


    def pulser_off(self):
        """ Switches the pulsing device off.

        @return int: error code (0:OK, -1:error)
        """
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP)
        spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 0)  # enable analog output 1
        spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 0)  # enable analog output 1
        err = self._readOutError()
        if err:
            return -1
        else:
            return 0


    def upload_asset(self, asset_name=None):
        """ Upload an already hardware conform file to the device mass memory.
            Also loads these files into the device workspace if present.
            Does NOT load waveforms/sequences/patterns into channels.

        @param asset_name: string, name of the ensemble/sequence to be uploaded

        @return int: error code (0:OK, -1:error)

        If nothing is passed, method will be skipped.

        This method has no effect when using pulser hardware without own mass memory
        (i.e. PulseBlaster, FPGA)
        """
        self.log.warning('Uploading assets is not available for the Spectrum M4i AWG series!\n'
                         'Method call will be ignored.')
        return 0


    def load_asset(self, asset_name, load_dict=None):
        """ Loads a sequence or waveform to the specified channel of the pulsing device.
        For devices that have a workspace (i.e. AWG) this will load the asset from the device
        workspace into the channel.
        For a device without mass memory this will transfer the waveform/sequence/pattern data
        directly to the device so that it is ready to play.

        @param str asset_name: The name of the asset to be loaded

        @param dict load_dict:  a dictionary with keys being one of the available channel numbers
                                and items being the name of the already sampled waveform/sequence
                                files.
                                Examples:   {1: rabi_Ch1, 2: rabi_Ch2}
                                            {1: rabi_Ch2, 2: rabi_Ch1}
                                This parameter is optional. If none is given then the channel
                                association is invoked from the file name, i.e. the appendix
                                (_ch1, _ch2 etc.)

        @return int: error code (0:OK, -1:error)
        """

        try:
            asset_name
        except NameError:
            pass
        else:
            print(asset_name)

        # setup the mode
        llLoops = int64(0)  # loop continuously
        llMemSamples = int64(KILO_B(64))
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SINGLE)  # continuous
        spcm_dwSetParam_i64(self._hCard, SPC_MEMSIZE, llMemSamples)  # replay lenth
        spcm_dwSetParam_i64(self._hCard, SPC_LOOPS, llLoops)  # number of repetitions

        # setup the trigger mode -> no trigger
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE);
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE) # (SW trigger, no output)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, 0)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, 0)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, 0)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, 0)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, 0)
        #spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, 0)

        a_ch1_signal = (np.sin(np.arange(llMemSamples.value)/llMemSamples.value * 2*np.pi) *(2**15-1)).astype(dtype=np.int16)
        a_ch2_signal = (np.cos(np.arange(llMemSamples.value)/llMemSamples.value * 2*np.pi) *(2**15-1)).astype(dtype=np.int16)
        d_ch1_signal = (np.sin(np.arange(llMemSamples.value)/llMemSamples.value*2 * 2*np.pi) < 0).astype(dtype=np.bool)
        d_ch2_signal = (np.sin(np.arange(llMemSamples.value)/llMemSamples.value*4 * 2*np.pi) < 0).astype(dtype=np.bool)
        d_ch3_signal = (np.sin(np.arange(llMemSamples.value)/llMemSamples.value*8 * 2*np.pi) < 0).astype(dtype=np.bool)

        pvBuffer, qwBufferSize = self._create_combined_buffer_data(llMemSamples, {
            'a_ch1': a_ch1_signal, 'a_ch2': a_ch2_signal,
            'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal, 'd_ch3': d_ch3_signal})

        # we define the buffer for transfer and start the DMA transfer
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0),
                               qwBufferSize)
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

        self._readOutError()
        return 0


    def get_loaded_asset(self):
        """ Retrieve the currently loaded asset name of the device.

        @return str: Name of the current asset ready to play. (no filename)
        """
        return ""


    def clear_all(self):
        """ Clears all loaded waveforms from the pulse generators RAM/workspace.

        @return int: error code (0:OK, -1:error)
        """
        return -1


    def get_status(self):
        """ Retrieves the status of the pulsing hardware

        @return (int, dict): tuple with an interger value of the current status and a corresponding
                             dictionary containing status description for all the possible status
                             variables of the pulse generator hardware.

                             1h Acquisition modes only: the pretrigger area has been filled.
                             2h The first trigger has been detected.
                             4h The card has finished its run and is ready.
                             8h Multi/ABA/Gated acquisition of M4i/M4x/M2p only:
                                the pretrigger area of one segment has been filled.)
        """             
        
        status = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_M2STATUS, byref(status))
        self._readOutError()
        self.log.info('M4i6631x8 status: {0:x} (see manual for more information)'.format(status.value))

        return (0, {0 : 'dummy status'})


    def get_sample_rate(self):
        """ Get the sample rate of the pulse generator hardware

        @return float: The current sample rate of the device (in Hz)

        Do not return a saved sample rate from an attribute, but instead retrieve the current
        sample rate directly from the device.
        """

        sample_rate = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_SAMPLERATE, byref(sample_rate))
        self._readOutError()
        return float(sample_rate.value)


    def set_sample_rate(self, sample_rate):
        """ Set the sample rate of the pulse generator hardware.

        @param float sample_rate: The sampling rate to be set (in Hz)

        @return float: the sample rate returned from the device (in Hz).

        Note: After setting the sampling rate of the device, use the actually set return value for
              further processing.
        """
        constraints = self.get_constraints()
        if (sample_rate > constraints.sample_rate.max) or (sample_rate < constraints.sample_rate.min):
            self.log.info('Sample rate is out of constraints for the Spectrum M4i AWG series!\n'
                          'Method call will be ignored.')
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_SAMPLERATE, int32(int(sample_rate)))
            self._readOutError()

        return self.get_sample_rate()



    def get_analog_level(self, amplitude=None, offset=None):
        """ Retrieve the analog amplitude and offset of the provided channels.

        @param list amplitude: optional, if the amplitude value (in Volt peak to peak, i.e. the
                               full amplitude) of a specific channel is desired.
        @param list offset: optional, if the offset value (in Volt) of a specific channel is
                            desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel descriptor string
                               (i.e. 'a_ch1') and items being the values for those channels.
                               Amplitude is always denoted in Volt-peak-to-peak and Offset in volts.

        Note: Do not return a saved amplitude and/or offset value but instead retrieve the current
              amplitude and/or offset directly from the device.

        If nothing (or None) is passed then the levels of all channels will be returned. If no
        analog channels are present in the device, return just empty dicts.

        Example of a possible input:
            amplitude = ['a_ch1', 'a_ch4'], offset = None
        to obtain the amplitude of channel 1 and 4 and the offset of all channels
            {'a_ch1': -0.5, 'a_ch4': 2.0} {'a_ch1': 0.0, 'a_ch2': 0.0, 'a_ch3': 1.0, 'a_ch4': 0.0}

        The major difference to digital signals is that analog signals are always oscillating or
        changing signals, otherwise you can use just digital output. In contrast to digital output
        levels, analog output levels are defined by an amplitude (here total signal span, denoted in
        Voltage peak to peak) and an offset (a value around which the signal oscillates, denoted by
        an (absolute) voltage).

        In general there is no bijective correspondence between (amplitude, offset) and
        (value high, value low)!
        """

        ampl_ch1 = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_AMP0, byref(ampl_ch1))  # channel0 amplitude in mV
        self._readOutError()

        ampl_ch2 = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_AMP1, byref(ampl_ch2))  # channel1 amplitude in mV
        self._readOutError()

        ampl_dict = {'a_ch1': float(ampl_ch1.value*2/1000), 'a_ch2': float(ampl_ch2.value*2/1000)}

        if amplitude is not None:
            ampl_dict = {key: ampl_dict[key] for key in amplitude}
        if offset is not None:
            offset_dict = {key: 0.0 for key in offset}
        else:
            offset_dict = {'a_ch1' : 0.0, 'a_ch2' : 0.0}

        return ampl_dict, offset_dict


    def set_analog_level(self, amplitude=None, offset=None):
        """ Set amplitude and/or offset value of the provided analog channel(s).

        @param dict amplitude: dictionary, with key being the channel descriptor string
                               (i.e. 'a_ch1', 'a_ch2') and items being the amplitude values
                               (in Volt peak to peak, i.e. the full amplitude) for the desired
                               channel.
        @param dict offset: dictionary, with key being the channel descriptor string
                            (i.e. 'a_ch1', 'a_ch2') and items being the offset values
                            (in absolute volt) for the desired channel.

        @return (dict, dict): tuple of two dicts with the actual set values for amplitude and
                              offset for ALL channels.

        If nothing is passed then the command will return the current amplitudes/offsets.

        Note: After setting the amplitude and/or offset values of the device, use the actual set
              return values for further processing.

        The major difference to digital signals is that analog signals are always oscillating or
        changing signals, otherwise you can use just digital output. In contrast to digital output
        levels, analog output levels are defined by an amplitude (here total signal span, denoted in
        Voltage peak to peak) and an offset (a value around which the signal oscillates, denoted by
        an (absolute) voltage).

        In general there is no bijective correspondence between (amplitude, offset) and
        (value high, value low)!
        """

        if offset is not None:
            self.log.warning('Setting analog offset values is not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')

        constraints = self.get_constraints()
        if 'a_ch1' in amplitude.keys():
            if (amplitude['a_ch1'] >= constraints.a_ch_amplitude.min) and\
                (amplitude['a_ch1'] <= constraints.a_ch_amplitude.max):

                spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(int(amplitude['a_ch1']*1000/2)))  # channel0 amplitude in mV
                self._readOutError()
            else:
                self.log.warning('Amplitude voltage level for the analog output channel 1 is out of constraints for the'
                                 'Spectrum M4i AWG series!\nMethod call will be ignored.')
        if 'a_ch2' in amplitude.keys():
            if (amplitude['a_ch2'] >= constraints.a_ch_amplitude.min) and\
                (amplitude['a_ch2'] <= constraints.a_ch_amplitude.max):

                spcm_dwSetParam_i32(self._hCard, SPC_AMP1, int32(int(amplitude['a_ch2']*1000/2)))  # channel0 amplitude in mV
                self._readOutError()
            else:
                self.log.warning('Amplitude voltage level for the analog output channel 2 is out of constraints for the'
                                 'Spectrum M4i AWG series!\nMethod call will be ignored.')

        return self.get_analog_level()


    def get_digital_level(self, low=None, high=None):
        """ Retrieve the digital low and high level of the provided/all channels.

        @param list low: optional, if the low value (in Volt) of a specific channel is desired.
        @param list high: optional, if the high value (in Volt) of a specific channel is desired.

        @return: (dict, dict): tuple of two dicts, with keys being the channel descriptor strings
                               (i.e. 'd_ch1', 'd_ch2') and items being the values for those
                               channels. Both low and high value of a channel is denoted in volts.

        Note: Do not return a saved low and/or high value but instead retrieve
              the current low and/or high value directly from the device.

        If nothing (or None) is passed then the levels of all channels are being returned.
        If no digital channels are present, return just an empty dict.

        Example of a possible input:
            low = ['d_ch1', 'd_ch4']
        to obtain the low voltage values of digital channel 1 an 4. A possible answer might be
            {'d_ch1': -0.5, 'd_ch4': 2.0} {'d_ch1': 1.0, 'd_ch2': 1.0, 'd_ch3': 1.0, 'd_ch4': 4.0}
        Since no high request was performed, the high values for ALL channels are returned (here 4).

        The major difference to analog signals is that digital signals are either ON or OFF,
        whereas analog channels have a varying amplitude range. In contrast to analog output
        levels, digital output levels are defined by a voltage, which corresponds to the ON status
        and a voltage which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between (amplitude, offset) and
        (value high, value low)!
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

        @param dict low: dictionary, with key being the channel descriptor string
                         (i.e. 'd_ch1', 'd_ch2') and items being the low values (in volt) for the
                         desired channel.
        @param dict high: dictionary, with key being the channel descriptor string
                          (i.e. 'd_ch1', 'd_ch2') and items being the high values (in volt) for the
                          desired channel.

        @return (dict, dict): tuple of two dicts where first dict denotes the current low value and
                              the second dict the high value for ALL digital channels.
                              Keys are the channel descriptor strings (i.e. 'd_ch1', 'd_ch2')

        If nothing is passed then the command will return the current voltage levels.

        Note: After setting the high and/or low values of the device, use the actual set return
              values for further processing.

        The major difference to analog signals is that digital signals are either ON or OFF,
        whereas analog channels have a varying amplitude range. In contrast to analog output
        levels, digital output levels are defined by a voltage, which corresponds to the ON status
        and a voltage which corresponds to the OFF status (both denoted in (absolute) voltage)

        In general there is no bijective correspondence between (amplitude, offset) and
        (value high, value low)!
        """
        self.log.warning('Setting the digital voltage levels is not available for the Spectrum M4i AWG series!\n'
                         'Method call will be ignored.')
        return self.get_digital_level()


    def get_active_channels(self, ch=None):
        """ Get the active channels of the pulse generator hardware.

        @param list ch: optional, if specific analog or digital channels are needed to be asked
                        without obtaining all the channels.

        @return dict:  where keys denoting the channel string and items boolean expressions whether
                       channel are active or not.

        Example for an possible input (order is not important):
            ch = ['a_ch2', 'd_ch2', 'a_ch1', 'd_ch5', 'd_ch1']
        then the output might look like
            {'a_ch2': True, 'd_ch2': False, 'a_ch1': False, 'd_ch5': True, 'd_ch1': False}

        If no parameter (or None) is passed to this method all channel states will be returned.
        """
        #a_ch1 and a_ch2
        a_channel = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPC_CHENABLE, byref(a_channel))
        self._readOutError()

        self._active_channels['a_ch1'] = bool(a_channel.value & CHANNEL0)
        self._active_channels['a_ch2'] = bool(a_channel.value & CHANNEL1)

        #digital channels
        # FIXME -> configfile
        """
        digital channels are encoded in analog samples for synchonous readout
        -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)

        current implementation:
        analog channel 1: bit 0-13: a_ch1
                          bit   15: d_ch1
                          bit   14: d_ch2
        analog channel 2: bit 0-14: a_ch2
                          bit   15: d_ch3
        """

        #d_ch1
        x0_mode = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPCM_X0_MODE, byref(x0_mode))
        self._readOutError()
        self._active_channels['d_ch1']\
            = bool(x0_mode.value & SPCM_XMODE_DIGOUT) and bool(x0_mode.value & SPCM_XMODE_DIGOUTSRC_CH0)

        # d_ch2
        x1_mode = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPCM_X1_MODE, byref(x1_mode))
        self._readOutError()
        self._active_channels['d_ch2'] \
            = bool(x1_mode.value & SPCM_XMODE_DIGOUT) and bool(x1_mode.value & SPCM_XMODE_DIGOUTSRC_CH0)

        # d_ch3
        x2_mode = int32(0)
        spcm_dwGetParam_i32(self._hCard, SPCM_X2_MODE, byref(x2_mode))
        self._readOutError()
        self._active_channels['d_ch3'] \
            = bool(x2_mode.value & SPCM_XMODE_DIGOUT) and bool(x2_mode.value & SPCM_XMODE_DIGOUTSRC_CH1)

        if ch is None:
            return self._active_channels
        else:
            return {key: self._active_channels[key] for key in ch.keys()}


    def set_active_channels(self, ch=None):
        """ Set the active channels for the pulse generator hardware.

        @param dict ch: dictionary with keys being the analog or digital string generic names for
                        the channels (i.e. 'd_ch1', 'a_ch2') with items being a boolean value.
                        True: Activate channel, False: Deactivate channel

        @return dict: with the actual set values for ALL active analog and digital channels

        If nothing is passed then the command will simply return the unchanged current state.

        Note: After setting the active channels of the device,
              use the returned dict for further processing.

        Example for possible input:
            ch={'a_ch2': True, 'd_ch1': False, 'd_ch3': True, 'd_ch4': True}
        to activate analog channel 2 digital channel 3 and 4 and to deactivate
        digital channel 1.

        The hardware itself has to handle, whether separate channel activation is possible.
        """

        if ch is None:
            return self.get_active_channels()

        for channel in ch:
            if channel in self._active_channels:
                self._active_channels[channel] = ch[channel]
            else:
                self.log.error('Trying to (de)activate channel "{0}". This channel is not present '
                               'in AWG. Setting channels aborted.'.format(channel))

        # analog channels
        channel0 = CHANNEL0 * self._active_channels['a_ch1']
        channel1 = CHANNEL1 * self._active_channels['a_ch2']

        if (self._active_channels['d_ch1'] or self._active_channels['d_ch2']) \
            and (not self._active_channels['a_ch1']):
            self.log.warning('a_ch1 disable state is overwritten because d_ch1 or d_ch2 is enabled')
            channel0 = CHANNEL0

        if self._active_channels['d_ch3'] and (not self._active_channels['a_ch2']):
            self.log.warning('a_ch2 disable state is overwritten because d_ch3 is enabled')
            channel1 = CHANNEL1

        if channel0 | channel1:
            spcm_dwSetParam_i32(self._hCard, SPC_CHENABLE, int32(channel0 | channel1))
            self._readOutError()
        else:
            self.log.info('Disabling all channels is not implemented for the Spectrum M4i AWG series!')

        #digital channels
        #FIXME -> configfile
        """
        digital channels are encoded in analog samples for synchonous readout
        -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)
        
        current implementation:
        analog channel 1: bit 0-13: a_ch1
                          bit   15: d_ch1
                          bit   14: d_ch2
        analog channel 2: bit 0-14: a_ch2
                          bit   15: d_ch3
        """
        #d_ch1
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
        return self._active_channels


    def get_uploaded_asset_names(self):
        """ Retrieve the names of all uploaded assets on the device.

        @return list: List of all uploaded asset name strings in the current device directory.
                      This is no list of the file names.

        Unused for pulse generators without sequence storage capability (PulseBlaster, FPGA).
        """
        self.log.warning('Asset names are not available for the Spectrum M4i AWG series!\n'
                         'Method call will be ignored.')
        return []


    def get_saved_asset_names(self):
        """ Retrieve the names of all sampled and saved assets on the host PC. This is no list of
            the file names.

        @return list: List of all saved asset name strings in the current
                      directory of the host PC.
        """
        self.log.warning('Asset names are not available for the Spectrum M4i AWG series!\n'
                         'Method call will be ignored.')
        return []


    def delete_asset(self, asset_name):
        """ Delete all files associated with an asset with the passed asset_name from the device
            memory (mass storage as well as i.e. awg workspace/channels).

        @param str asset_name: The name of the asset to be deleted
                               Optionally a list of asset names can be passed.

        @return list: a list with strings of the files which were deleted.

        Unused for pulse generators without sequence storage capability (PulseBlaster, FPGA).
        """
        if asset_name:
            self.log.warning('Asset names are not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')
        return []


    def set_asset_dir_on_device(self, dir_path):
        """ Change the directory where the assets are stored on the device.

        @param str dir_path: The target directory

        @return int: error code (0:OK, -1:error)

        Unused for pulse generators without changeable file structure (PulseBlaster, FPGA).
        """
        if dir_path:
            self.log.warning('Asset dir paths are not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')
        return 0


    def get_asset_dir_on_device(self):
        """ Ask for the directory where the hardware conform files are stored on the device.

        @return str: The current file directory

        Unused for pulse generators without changeable file structure (i.e. PulseBlaster, FPGA).
        """
        if dir_path:
            self.log.warning('Asset dir paths are not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')
        return ''


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
            self.log.warning('Interleave mode not available for the Spectrum M4i AWG series!\n'
                             'Method call will be ignored.')
        return False


    def tell(self, command):
        """ Sends a command string to the device.

        @param string command: string containing the command

        @return int: error code (0:OK, -1:error)
        """
        pass


    def ask(self, question):
        """ Asks the device a 'question' and receive and return an answer from it.
a
        @param string question: string containing the command

        @return string: the answer of the device to the 'question' in a string
        """
        pass


    def reset(self):
        """ Reset the device.

        @return int: error code (0:OK, -1:error)
        """

        self._readOutError()
        if spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_RESET):
            self._readOutError()
            return -1
        else:
            return 0


    def has_sequence_mode(self):
        """ Asks the pulse generator whether sequence mode exists.

        @return: bool, True for yes, False for no.
        """

        availableModes = int32(0)
        spcm_dwSetParam_i32(self._hCard, SPC_AVAILCARDMODES, byref(availableModes))
        self._readOutError()
        if availableModes.value & SPC_REP_STD_SEQUENCE:
            return True
        else:
            return False

    def _readOutError(self):
        """checks the error state and prints it out. Errors must be read out before the AWG
        can accept further commands"""
        spcm_dwGetParam_i32
        errortext = (ctypes.c_char*ERRORTEXTLEN)()
        err = spcm_dwGetErrorInfo_i32(self._hCard, None, None, errortext)
        if err:
            self.log.warning(errortext.value)
        return err

    def _spcm_dwGetParam_i32(self, lRegister):

        value = int32(0)
        spcm_dwGetParam_i32(self._hCard, lRegister, byref(value))
        err = self._readOutError()
        if err:
           return -1
        else:
            return value.value

    def _spcm_dwSetParam_i32(self, lRegister, plValue):
        spcm_dwSetParam_i32(self._hCard, lRegister, int32(plValue))
        err = self._readOutError()
        if err:
           return -1
        else:
            return 0

    def _create_combined_buffer_data(self, llMemSamples, channel_data):
        """
        digital channels are encoded in analog samples for synchonous readout
        -> analog channel's resolution is reduced by 1 bit per digital channel (max 3 bits)

        #FIXME -> configfile
        current implementation:
        analog channel 1: bit 0-13: a_ch1
                          bit   15: d_ch1
                          bit   14: d_ch2
        analog channel 2: bit 0-14: a_ch2
                          bit   15: d_ch3
        """
        active_channels = self.get_active_channels()

        for key in channel_data.keys():
            if channel_data[key].size != llMemSamples.value:
                self.log.error('Channel data for channel {} has the wrong size!'.format(key))
                return -1


        if not 'a_ch1' in channel_data.keys():
            a_ch1_data = np.zeros(llMemSamples.value, dtype=np.dtype(np.int16))
        else:
            a_ch1_data = channel_data['a_ch1']

        if not 'a_ch2' in channel_data.keys():
            a_ch2_data = np.zeros(llMemSamples.value, dtype=np.dtype(np.int16))
        else:
            a_ch2_data = channel_data['a_ch2']

        if not 'd_ch1' in channel_data.keys():
            d_ch1_data = np.zeros(llMemSamples.value, dtype=np.dtype('?'))
        else:
            d_ch1_data = channel_data['d_ch1']

        if not 'd_ch2' in channel_data.keys():
            d_ch2_data = np.zeros(llMemSamples.value, dtype=np.dtype('?'))
        else:
            d_ch2_data = channel_data['d_ch2']

        if not 'd_ch3' in channel_data.keys():
            d_ch3_data = np.zeros(llMemSamples.value, dtype=np.dtype('?'))
        else:
            d_ch3_data = channel_data['d_ch3']

        if active_channels['d_ch2']:
            a_ch1_data = (a_ch1_data.astype(np.uint16) >> 2) | (d_ch1_data << 15) | (d_ch2_data << 14)
        elif active_channels['d_ch1']:
            a_ch1_data = (a_ch1_data.astype(np.uint16) >> 1) | (d_ch1_data << 15)

        if active_channels['d_ch3']:
            a_ch2_data = (a_ch2_data.astype(np.uint16) >> 1) | (d_ch3_data << 15)

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
        chcount = int32(0)  # number of activated channels
        spcm_dwGetParam_i32(self._hCard, SPC_CHCOUNT, byref(chcount))

        lBytesPerSample = int32(0)  # bytes per sample
        spcm_dwGetParam_i32(self._hCard, SPC_MIINST_BYTESPERSAMPLE, byref(lBytesPerSample))

        qwBufferSize = uint64(llMemSamples.value * lBytesPerSample.value * chcount.value)
        pvBuffer = create_string_buffer(qwBufferSize.value)

        # calculate the data
        pnBuffer = cast(pvBuffer, ptr16)
        for i in range(0, llMemSamples.value * chcount.value, 1):
            pnBuffer[i] = int16(combined_channel_data[i])

        return pvBuffer, qwBufferSize


    def _testfunction(self, channel, amplitude, offset, numberofsamples=int64(KILO_B(64))):
        """we try to produce a sine wave"""

        # setup the mode
        #qwChEnable = uint64(channel+1)
        llMemSamples = numberofsamples
        self.log.info(llMemSamples)
        llLoops = int64(3)  # loop continuously
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
                pnBuffer[i] = int((np.sin(i/int(llMemSamples.value) * 2*np.pi) + offset)*amplitude)
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
        self._readOutError()
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE);
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE) # (SW trigger, no output)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, 0)
        # spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, 0)

        # Setting up the card mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE); # enable sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, 2); # Divide on - board mem in two parts
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_STARTSTEP, 0); # Step#0 is the first step after card start

        # Setting up the data memory and transfer data
        llMemSamples = int64(KILO_B(64))
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 0); # set current configuration switch to segment 0
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, llMemSamples); # define size of current segment 0

        sample_rate = self.get_sample_rate()
        self._readOutError()

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
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize);
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        # Setting up the data memory and transfer data
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 1); # set current configuration switch to segment 1
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, llMemSamples); # define size of current segment 1

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
        spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize);
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        # Setting up the sequence memory(Only two steps used here as an example)
        lStep = 0; # current step is Step  # 0
        llSegment = 0; # associated with memory
        llLoop = 1; # Pattern will be repeated 10 times
        llNext = 1; # Next step is Step  # 1
        llCondition = SPCSEQ_ENDLOOPONTRIG; # Unconditionally leave current step

        # combine all the parameters to one int64 bit value
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment));
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue);

        lStep = 1; # current step is Step  # 1
        llSegment = 1; # associated with memory segment 1
        llLoop = 1; # Pattern will be repeated once before condition is checked
        llNext = 0; # Next step is Step  # 0
        llCondition = SPCSEQ_ENDLOOPONTRIG; # Repeat current step until a trigger has occurred
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment));
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue);

        # Start the card
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);

        # ... wait here or do something else ...

        #Stop the card
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_STOP);
