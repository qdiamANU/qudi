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

import os
import sys
import time
import re
import numpy as np
import ctypes
import matplotlib.pyplot as plt

from thirdparty.spectrum_instruments.pyspcm import *

class AWGDebugClass():
    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        # open card
        hCard = spcm_hOpen(create_string_buffer(b'/dev/spcm0'))
        self._hCard = hCard
        if hCard == None:
            print('no card found...')
            return 1

        # prepare the card in a defined state
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_RESET)

        # read type, function and sn and check for A/D card
        lCardType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCITYP, byref(lCardType))
        lSerialNumber = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_PCISERIALNO, byref(lSerialNumber))
        lFncType = int32(0)
        spcm_dwGetParam_i32(hCard, SPC_FNCTYPE, byref(lFncType))

        # sCardName = 'M4i.6631-x8'
        spcm_dwSetParam_i32(hCard, SPC_TIMEOUT, 30000)  # timeout 15 s

        self._active_channels = {'a_ch1': True,
                                 'a_ch2': True,
                                 'd_ch1': True,
                                 'd_ch2': True,
                                 'd_ch3': True,}

        self.set_active_channels(self._active_channels)
        self.pulser_off()

    def on_deactivate(self):
        """ Required tasks to be performed during deactivation of the module.
        """
        if self._hCard != None:
            spcm_vClose(self._hCard)
        return

    def pulser_on(self):
        """ Switches the pulsing device on.

        @return int: error code (0:OK, -1:error)
        """

        if self._spcm_dwGetParam_i32(SPC_ENABLEOUT0) == 0:
            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT0, 1)  # enable analog output 0
        if self._spcm_dwGetParam_i32(SPC_ENABLEOUT1) == 0:
            spcm_dwSetParam_i64(self._hCard, SPC_ENABLEOUT1, 1)  # enable analog output 1

        if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
            # spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_FORCETRIGGER)
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

        if self.read_out_error():
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

        if self.read_out_error():
            return -1
        else:
            return 0


    def set_analog_level(self, amplitude=None):
        spcm_dwSetParam_i32(self._hCard, SPC_AMP0, int32(int(amplitude['a_ch1'] * 1000 / 2)))
        spcm_dwSetParam_i32(self._hCard, SPC_AMP1, int32(int(amplitude['a_ch2'] * 1000 / 2)))
        self.read_out_error()
        return

    def set_active_channels(self, ch=None):
        if ch is None:
            return self.get_active_channels()

        for channel in ch:
            if channel in self._active_channels:
                self._active_channels[channel] = ch[channel]
            else:
                print('Trying to (de)activate channel "{0}". This channel is not present '
                               'in AWG. Setting channels aborted.'.format(channel))

        # analog channels
        channel0 = CHANNEL0 * (self._active_channels['a_ch1'] or self._active_channels['d_ch1'] or
                               self._active_channels['d_ch2'])
        channel1 = CHANNEL1 * (self._active_channels['a_ch2'] or self._active_channels['d_ch3'])
        spcm_dwSetParam_i32(self._hCard, SPC_CHENABLE, channel0 | channel1)

        # digital channels
        """
        analog channel 1: bit 0-13: a_ch1
                          bit   15: d_ch1
                          bit   14: d_ch2
        analog channel 2: bit 0-14: a_ch2
                          bit   15: d_ch3
        """
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

        self.read_out_error()
        return self._active_channels

    def load_sequence(self, channel_data, loop):
        """ Loads a sequence to the awg

        @param channel_data: list of dictionaries with channel_data:
                             [{'a_ch1' : np.array([...]), 'd_ch1' : np.array([])},
                              {'a_ch1' : np.array([...]), 'd_ch1' : np.array([])}
                             ]

        @return int: error code (0:OK, -1:error)
        """

        self.pulser_off()

        max_segments = len(channel_data)

        # Setting up the trigger
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_ANDMASK, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ORMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK0, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIG_CH_ANDMASK1, SPC_TMASK_NONE)
        spcm_dwSetParam_i32(self._hCard, SPC_TRIGGEROUT, SPC_TMASK_NONE)

        # Setting up the card mode
        spcm_dwSetParam_i32(self._hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)  # enable sequence mode
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_MAXSEGMENTS, max_segments * 2)  # Divide on - board mem in parts
        spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_STARTSTEP, 0)  # Step#0 is the first step after card start

        # Setting up the data memory and transfer data
        for i in range(max_segments):
            llMemSamples = int64(list(channel_data[i].values())[0].size)

            # set current configuration switch to segment 2*i
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 2 * i)
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, llMemSamples)  # define size of current segment 0

            pvBuffer, qwBufferSize = self._create_combined_buffer_data(llMemSamples, channel_data[i])

            # it is assumed, that the Buffer memory has been allocated and is already filled with valid data
            spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

            # set current configuration switch to segment 2*i+1
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, 2 * i + 1)
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE,
                                int64(384))  # define size of current segment 2*i+1

            # just 384 zeros for the timetagger
            pvBuffer, qwBufferSize = self._create_combined_buffer_data(int64(384), {})

            # it is assumed, that the Buffer memory has been allocated and is already filled with valid data
            spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

            print('Uploaded segment {} of {}.'.format(i + 1, max_segments))

        # Setting up the sequence memory
        for i in range(max_segments * 2):
            lStep = i  # current step is Step  #i
            llSegment = i  # associated with memory

            if i % 2 == 0:
                llLoop = loop  # Pattern will be repeated 1 times
            else:
                llLoop = 1  # Pattern will be repeated 1 time

            llNext = (i + 1) % (max_segments * 2)  # Next step is Step  # i+1

            if i != (max_segments * 2) - 1:
                llCondition = SPCSEQ_ENDLOOPALWAYS  # Unconditionally leave current step
            else:
                llCondition = SPCSEQ_END  # End of sequence

            # print(lStep, llSegment, llLoop, llNext, llCondition)
            # combine all the parameters to one int64 bit value
            llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
            spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

            llValue_test = int64(0)
            spcm_dwGetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, byref(llValue_test))

            if not llValue.value == llValue_test.value:
                print(
                    'Error! {:016x} != {:016x}'.format(np.uint64(llValue.value), np.uint64(llValue_test.value)))

        if self.read_out_error():
            return -1
        else:
            return 0

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

        status = self._spcm_dwGetParam_i32(SPC_M2STATUS)
        self.log.info('M4i6631x8 status: {0:x} (see manual for more information)'.format(status))
        self.read_out_error()
        return (status.value, {status.value: ''})

    def force_trigger(self, wait=False):
        if self._spcm_dwGetParam_i32(SPC_M2STATUS) & M2STAT_CARD_READY:
            return -1
        if wait:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER | M2CMD_CARD_WAITREADY)
        else:
            spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER)

        if self.read_out_error():
            return -1
        else:
            return 0

    def read_out_error(self):
        """checks the error state and prints it out. Errors must be read out before the AWG
        can accept further commands"""
        errortext = (ctypes.c_char * ERRORTEXTLEN)()
        err = spcm_dwGetErrorInfo_i32(self._hCard, None, None, errortext)
        if err:
            print(err, errortext.value)
        return err

    def _spcm_dwGetParam_i32(self, lRegister):
        value = int32(0)
        spcm_dwGetParam_i32(self._hCard, lRegister, byref(value))
        if self.read_out_error():
            return -1
        else:
            return value.value

    def _spcm_dwSetParam_i32(self, lRegister, plValue):
        spcm_dwSetParam_i32(self._hCard, lRegister, int32(plValue))
        if self.read_out_error():
            return -1
        else:
            return 0

    def _create_combined_buffer_data(self, llMemSamples, channel_data):
        """
        analog channel 1: bit 0-13: a_ch1
                          bit   15: d_ch1
                          bit   14: d_ch2
        analog channel 2: bit 0-14: a_ch2
                          bit   15: d_ch3
        """
        llMemSamples
        active_channels = self._active_channels

        for key in channel_data.keys():
            if channel_data[key].size != llMemSamples.value:
                print('Channel data for channel {} has the wrong size!'.format(key))
                return -1

        for key in active_channels.keys():
            if key not in channel_data.keys():
                channel_data[key] = np.zeros(llMemSamples.value, dtype=np.dtype(np.int16))

        if active_channels['d_ch2']:
            a_ch1_data = ((channel_data['a_ch1'].astype(np.uint16) >> 2)  |
                          (channel_data['d_ch1'].astype(bool)      << 15) |
                          (channel_data['d_ch2'].astype(bool)      << 14))
        elif active_channels['d_ch1']:
            a_ch1_data = ((channel_data['a_ch1'].astype(np.unint16) >> 1) |
                          (channel_data['d_ch1'].astype(bool)       << 15))

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
            print('No channel activated.')
            return -1

        # setup software buffer
        chcount = self._spcm_dwGetParam_i32(SPC_CHCOUNT)
        lBytesPerSample = self._spcm_dwGetParam_i32(SPC_MIINST_BYTESPERSAMPLE)
        qwBufferSize = uint64(llMemSamples.value * lBytesPerSample * chcount)
        pvBuffer     = create_string_buffer(qwBufferSize.value)

        # calculate the data
        pnBuffer = cast(pvBuffer, ptr16)
        for i in range(0, llMemSamples.value * chcount, 1):
            pnBuffer[i] = int16(combined_channel_data[i])

        return pvBuffer, qwBufferSize

    def set_list(self, frequency, amplitude= 4.0, llMemSamples=1024*256, sample_rate=1.25e9, loop=1000):
        """

        @param frequency: list of frequencies in Hz
        @amplitude:       amplitude (peak-peak) in Volt

        @return int: error code (0:OK, -1:error)
        """
        self.pulser_off()
        self.set_active_channels(ch={'a_ch1': True, 'a_ch2': False,
                                                 'd_ch1': True, 'd_ch2': True, 'd_ch3': False})
        self.set_analog_level(amplitude={'a_ch1': amplitude, 'a_ch2': amplitude})
        spcm_dwSetParam_i32(self._hCard, SPC_SAMPLERATE, int32(int(sample_rate)))

        channel_data = []
        for i in range(len(frequency)):
            t = np.arange(llMemSamples)/sample_rate
            a_ch1_signal = (np.sin(t * (frequency[i]) * 2*np.pi) * (2**15-1)).astype(dtype=np.int16)
            d_ch1_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)

            if i == 0:
                d_ch2_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)
            else:
                d_ch2_signal = np.full(llMemSamples, 0).astype(dtype=np.bool)

            channel_data.append({'a_ch1': a_ch1_signal, 'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal})
        self.load_sequence(channel_data, loop)

        self.channel_data = channel_data

    def restart(self):
        """
        Reset of MW list mode position to start (first frequency step)

        @return int: error code (0:OK, -1:error)
        """

        self.read_out_error()
        self.pulser_on()
        self.force_trigger(wait=True)
        return 0
