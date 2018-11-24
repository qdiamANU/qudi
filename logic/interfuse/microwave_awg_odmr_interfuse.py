# -*- coding: utf-8 -*-

"""
This file contains the Qudi Interfuse between Microwave hardware and AWG hardware.
It is for use in CW ODMR only - it is not appropriate for pulsed operations
Todo: many parameters are currently hard-coded for the Spectrum M4i6631-x8 AWG

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
from core.module import Connector, StatusVar
from core.module import Base, ConfigOption
from logic.generic_logic import GenericLogic
from interface.microwave_interface import MicrowaveInterface
from interface.microwave_interface import MicrowaveLimits
from interface.microwave_interface import MicrowaveMode
from interface.microwave_interface import TriggerEdge

from thirdparty.spectrum_instruments.py_header.regs import *


class MicrowaveAwgInterfuse(GenericLogic, MicrowaveInterface):

    _modclass = 'MicrowaveAwgInterfuse'
    _modtype = 'interfuse'

    # declare connectors, here you can see the interfuse action: the in
    # connector will cope a motor hardware, that means a motor device can
    # connect to the in connector of the logic.
    microwave = Connector(interface='MicrowaveInterface')
    awg = Connector(interface='PulserInterface')

    # additional stuff for AWG-microwave-combo
    #  fixme: these should be soft-coded
    awg_amplitude = StatusVar('awg_amplitude', 4.0)
    awg_offset_frequency = StatusVar('awg_offset', 0.0e6)
    number_of_samples = StatusVar('number_of_samples', 256)
    number_of_loops = StatusVar('number_of_loops', 1000)
    awg_sample_rate = StatusVar('awg_sample_rate', 1.25e9)
    awg_range = StatusVar('awg_range', 400e6)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._mode = MicrowaveMode.CW
        self._freq_list = []

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        self._microwave_device = self.microwave()
        self._awg_device = self.awg()

        return 0

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        return 0

    def off(self):
        """
        Switches off any microwave output.
        Must return AFTER the device is actually stopped.

        @return int: error code (0:OK, -1:error)
        """
        return_val_1 = self._microwave_device.off()
        return_val_2 = self._awg_device.pulser_off()
        if (return_val_1 == 0) and (return_val_2 == 0):
            return 0
        else:
            return -1

    def get_status(self):
        """
        Gets the current status of the MW source, i.e. the mode (cw, list or sweep) and
        the output state (stopped, running)

        @return str, bool: mode ['cw', 'list', 'sweep'], is_running [True, False]
        """
        _, mw_is_running = self._microwave_device.get_status()

        if self._mode == MicrowaveMode.CW:
            mode = 'cw'
        elif (self._mode == MicrowaveMode.SWEEP) or (self._mode == MicrowaveMode.ASWEEP):
            mode = 'sweep'
        elif self._mode == MicrowaveMode.LIST:
            mode = 'list'
        else:
            self.log.error('Unsupported microwave mode: {}'.format(self._mode))
            mode = 'cw'

        return mode, mw_is_running

    def get_power(self):
        """
        Gets the microwave output power for the currently active mode.

        @return float: the output power in dBm
        """
        return self._microwave_device.get_power()

    def get_frequency(self):
        """
        Gets the frequency of the microwave output.
        Returns single float value if the device is in cw mode.
        Returns list like [start, stop, step] if the device is in sweep mode.
        Returns list of frequencies if the device is in list mode.

        @return [float, list]: frequency(s) currently set for this device in Hz
        """
        if self._mode == MicrowaveMode.CW:
            return self._microwave_device.get_frequency()
        elif (self._mode == MicrowaveMode.SWEEP) or (self._mode == MicrowaveMode.ASWEEP):
            return [self._freq_list[0], self._freq_list[-1], self._freq_list[1] - self._freq_list[0]]
        elif self._mode == MicrowaveMode.LIST:
            return self._freq_list
        else:
            self.log.error('Unsupported microwave mode: {}'.format(self._mode))
            return -1

    def cw_on(self):
        """
        Switches on cw microwave output.
        Must return AFTER the device is actually running.

        @return int: error code (0:OK, -1:error)
        """
        return self._microwave_device.cw_on()

    def set_cw(self, frequency=None, power=None):
        """
        Configures the device for cw-mode and optionally sets frequency and/or power

        @param float frequency: frequency to set in Hz
        @param float power: power to set in dBm

        @return tuple(float, float, str): with the relation
            current frequency in Hz,
            current power in dBm,
            current mode
        """
        self._mode = MicrowaveMode.CW
        return self._microwave_device.set_cw(frequency=frequency, power=power)

    def list_on(self):
        """
        Switches on the list mode microwave output.
        Must return AFTER the device is actually running.

        @return int: error code (0:OK, -1:error)
        """

        self._microwave_device.cw_on()
        self._awg_device.pulser_on()
        return 0

    def set_list(self, frequency=None, power=None):
        """
        Configures the device for list-mode and optionally sets frequencies and/or power

        @param list frequency: list of frequencies in Hz
        @param float power: MW power of the frequency list in dBm

        @return list, float, str: current frequencies in Hz, current power in dBm, current mode
        """

        self._freq_list = frequency

        local_oscillator_freq = frequency[0] - self.awg_offset_frequency
        awg_freq_array = np.array(frequency) - local_oscillator_freq
        # awg_freq_list = list(np.array(frequency) - local_oscillator_freq)

        # sanity checks:
        if np.min(awg_freq_array) - local_oscillator_freq < 0:
            self.log.error('LO freq ({} MHz) larger than desired minimum MW freq ({} MHz). Need to adjust LO freq.'.format(local_oscillator_freq, np.min(awg_freq_array)))
            return frequency, power, MicrowaveMode.LIST

        if (frequency[-1] - frequency[0]) > awg_range:
            self.log.error('Sweep range ({} MHz) is larger than the awg bandwidth ({} MHz).'.format(frequency[-1] - frequency[0], awg_range))
            return frequency, power, MicrowaveMode.LIST

        if (frequency[-1] - local_oscillator_freq) > awg_range:
            self.log.error('Combined sweep range ({} MHz) and LO offset ({} MHz) is larger than the awg bandwidth ({} MHz).'.format(frequency[-1] - frequency[0], local_oscillator_freq, awg_range))
            return frequency, power, MicrowaveMode.LIST

        # set microwave source in cw mode to act as local oscillator
        self.set_cw(local_oscillator_freq, power)

        # prepare AWG in triggered sequence mode, and upload frequency list
        err = self._set_awg_list(awg_freq_array)

        # cycle sample rate to avoid bug when using Spectrum AWG in sequence mode:
        self._awg_device._cycle_sample_rate()

        if err:
            self.log.error('Could not set AWG in list mode. Returning to CW microwave mode.')
            self._mode = MicrowaveMode.CW
            power = self.get_power()
            frequency = self.get_frequency()
            return frequency, power, 'cw'
        else:
            self._mode = MicrowaveMode.LIST
            power = self.get_power()
            return frequency, power, 'list'

    def reset_listpos(self):
        """
        Reset of MW list mode position to start (first frequency step)

        @return int: error code (0:OK, -1:error)
        """
        # todo: not 100% sure this works - AH 2/11/2018
        if self._mode == MicrowaveMode.LIST:
            self._awg_device.pulser_off()
            self._awg_device.pulser_on()
            self._awg_device.force_trigger(wait=True)
            return 0
        else:
            return -1

    def sweep_on(self):
        """ Switches on the sweep mode.

        @return int: error code (0:OK, -1:error)
        """
        return self._microwave_device.sweep_on()

    def set_sweep(self, start=None, stop=None, step=None, power=None):
        """
        Configures the device for sweep-mode and optionally sets frequency start/stop/step
        and/or power

        @return float, float, float, float, str: current start frequency in Hz,
                                                 current stop frequency in Hz,
                                                 current frequency step in Hz,
                                                 current power in dBm,
                                                 current mode
        """
        return self._microwave_device.set_sweep(start=start, stop=stop, step=step, power=power)

    def reset_sweeppos(self):
        """
        Reset of MW sweep mode position to start (start frequency)

        @return int: error code (0:OK, -1:error)
        """
        return self._microwave_device.reset_sweeppos

    def set_ext_trigger(self, pol=TriggerEdge.RISING):
        """ Set the external trigger for this device with proper polarization.

        @param TriggerEdge pol: polarisation of the trigger (basically rising edge or falling edge)

        @return object: current trigger polarity [TriggerEdge.RISING, TriggerEdge.FALLING]
        """
        # Trigger is controlled by AWG
        return pol

    def trigger(self):
        """ Trigger the next element in the list or sweep mode programmatically.

        @return int: error code (0:OK, -1:error)

        Ensure that the Frequency was set AFTER the function returns, or give
        the function at least a save waiting time corresponding to the
        frequency switching speed.
        """
        self.log.warning('Software trigger is not implemented yet.')
        return -1

    def get_limits(self):
        """ Return the device-specific limits in a nested dictionary.

          @return MicrowaveLimits: Microwave limits object
        """
        limits = self._microwave_device.get_limits()
        limits.list_maxentries = self._awg_device.get_constraints().sequence_num.max
        limits.supported_modes = (MicrowaveMode.CW, MicrowaveMode.SWEEP, MicrowaveMode.ASWEEP, MicrowaveMode.LIST)
        return limits

    def _set_awg_list(self, awg_freq_list):
        """ Prepare the AWG in triggered sequence mode, and upload frequency list

        AWG MW output frequency is on analogue ch 1
            amplitude hard-coded to max (4.0 V)
        CW laser channel is d_ch1

        @param list awg_freq_list: list of output frequencies for the AWG

        @return int: error code (0:OK, -1:error)
        """


        # fixme: currently too much hard-coding.
        # fixme: need to integrate with pulser 3, e.g. to create waveform?? Not if using for CW ODMR only
        # fixme: for CW ODMR, need to integrate with scan trigger (traditionally from NI card to SMIQ)
        # would then need to use sequence mode, looping over each frequency-step section until receiving an external trigger to go to next section
        # alternative is to write something new with the Time Tagger, where the AWG triggers the TT. This wouldn't have to use sequence mode
        # might need a new time tagger script??

        frequency = awg_freq_list

        # check number of frequency steps, divide AWG memory into corresponding number of segments
        n_freq_steps = len(frequency)

        # todo: need to initialise AWG
        # set active AWG channels: only need a_ch1, d_ch1
        self._awg_device._active_channels = {'a_ch1': True,
                                 'a_ch2': False,
                                 'd_ch1': True,
                                 'd_ch2': False,
                                 'd_ch3': False, }
        self._awg_device.set_active_channels(self._awg_device._active_channels)

        # prepare AWG in sequence mode, and divide AWG memory into [next-power-of-two(n_freq_steps)] number of segments
        n_actual_segments = self._awg_device.prepare_in_sequence_mode(n_freq_steps)

        # create and upload data for each frequency step
        #  min data length is 384/192 for 1/2 channels, with a stepsize of 32
        # max data length is (Mem/n_channels) / SPC_SEQMODE_MAXSEGMENTS)



        # desired number of samples per segment
        default_samples = int64(500 * 32) # todo: check if this is an appropriate value
        mem_size = 2e9  # AWG memory size # todo: possible to read off card?
        max_samples = mem_size/n_actual_segments
        n_samples = min(default_samples, max_samples)

        # create list of waveforms for frequency list
        sample_rate = self.get_sample_rate()
        awg_waveformbuffer_list = list()
        for i in range(n_freq_steps):
            a_ch1_signal = (np.sin(np.arange(n_samples) * (frequency[i] / sample_rate) * 2 * np.pi) * (
                    2 ** 15 - 1)).astype(dtype=np.int16)

            d_ch1_signal = np.ones(n_samples).astype(dtype=np.bool)

            analog_samples = {'a_ch1': a_ch1_signal}
            digital_samples = {'d_ch1': d_ch1_signal}

            # fixme: this won't work until write_waveform is rewritten
            waveformbuffer, qwBufferSize = self._awg_device.write_waveform('freqstep_{}'.format(i), analog_samples, digital_samples, True, True, n_samples)
            awg_waveformbuffer_list.append(waveformbuffer)


        # Transfer data to AWG, segment by segment
        for i in range(n_freq_steps):
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_WRITESEGMENT, i)  # set current configuration switch to segment i
            spcm_dwSetParam_i32(self._hCard, SPC_SEQMODE_SEGMENTSIZE, qwBufferSize)  # define size of current segment i
            pvBuffer = awg_waveformbuffer_list[i]
            spcm_dwDefTransfer_i64(self._hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwBufferSize)

            # Setting up the sequence memory(Only two steps used here as an example)
            lStep = i  # current step is Step  # 0
            llSegment = i  # associated with memory
            llLoop = 0  # Pattern will be repeated 10 times
            if i == n_freq_steps-1:
                llNext = 0
            else:
                llNext = i+1  # Next step is Step  # 1
            llCondition = SPCSEQ_ENDLOOPONTRIG  # Unconditionally leave current step

            # combine all the parameters to one int64 bit value
            llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
            spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

        self._read_out_error()






        # Setting up the sequence memory(Only two steps used here as an example)
        lStep = 0  # current step is Step  # 0
        llSegment = 0  # associated with memory
        llLoop = 1  # Pattern will be repeated 10 times
        llNext = 1  # Next step is Step  # 1
        llCondition = SPCSEQ_ENDLOOPONTRIG  # Unconditionally leave current step

        # combine all the parameters to one int64 bit value
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        lStep = 1  # current step is Step  # 1
        llSegment = 1  # associated with memory segment 1
        llLoop = 1  # Pattern will be repeated once before condition is checked
        llNext = 0  # Next step is Step  # 0
        llCondition = SPCSEQ_ENDLOOPONTRIG  # Repeat current step until a trigger has occurred
        llValue = int64((llCondition << 32) | (llLoop << 32) | (llNext << 16) | (llSegment))
        spcm_dwSetParam_i64(self._hCard, SPC_SEQMODE_STEPMEM0 + lStep, llValue)

        # Start the card
        spcm_dwSetParam_i32(self._hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)
        print('start')
        # ... wait here or do something else ...










        llMemSamples = 1024*self.number_of_samples # fixme: no connection to desired ODMR settings
        llLoops = self.number_of_loops
        amplitude = self.awg_amplitude
        sample_rate = self.awg_sample_rate

        self._awg_device.pulser_off()
        self._awg_device.set_active_channels(ch={'a_ch1': True, 'a_ch2': False,
                                                 'd_ch1': True, 'd_ch2': True, 'd_ch3': False})
        self._awg_device.set_analog_level(amplitude={'a_ch1': amplitude})
        self._awg_device.set_sample_rate(sample_rate)

        channel_data = []
        t = np.arange(llMemSamples) / sample_rate
        for i in range(len(frequency)):
            a_ch1_signal = (np.sin(t * (frequency[i] * 2*np.pi) * (2**15-1)
                            ).astype(dtype=np.int16)
            d_ch1_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)

            #Fixme
            if i == 0:
                d_ch2_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)
            else:
                d_ch2_signal = np.full(llMemSamples, 0).astype(dtype=np.bool)

            channel_data.append({'a_ch1': a_ch1_signal, 'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal})

        # Fixme: don't generate all the data first and then upload it all at once -> needs a lot of memory
        self._awg_device.load_sequence(channel_data, llLoops) # fixme: now need to use write_waveform()

        err = self._awg_device._read_out_error()
        if err:
            return -1
        else:
            return 0
