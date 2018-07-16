# -*- coding: utf-8 -*-

"""
This file contains the Qudi Interfuse between Microwave hardware and AWG hardware.

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
from core.util.modules import get_main_dir
import ctypes
import os
import matplotlib.pyplot as plt
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
    awg_amplitude = StatusVar('awg_amplitude', 4.0)
    awg_offset = StatusVar('awg_offset', 0.0e6)
    number_of_samples = StatusVar('number_of_samples', 256)
    number_of_loops = StatusVar('number_of_loops', 1000)
    awg_sample_rate = StatusVar('awg_sample_rate', 1.25e9)

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

        freq_list = frequency
        self._freq_list = freq_list

        local_oscillator_freq = freq_list[0] - self.awg_offset
        awg_range = 400e6

        if (frequency[-1] - frequency[0]) > awg_range:
            self.log.error('Sweep range is larger than the awg bandwidth.')
            return frequency, power, MicrowaveMode.LIST

        self._freq_list = frequency

        # set microwave source in list mode
        self.set_cw(frequency[0]-local_oscillator_freq, power)

        err = self._set_awg_list()
        sample_rate = self._awg_device.get_sample_rate()
        self._awg_device.set_sample_rate(600e6)
        self._awg_device.set_sample_rate(sample_rate)
        self.log.warning('AWG sample rate cycled {} MS/s -> {} MS/s -> {} MS/s to fix "corrupted-output-bug"'.format(
            int(sample_rate/1e6), int(600), int(sample_rate/1e6)))

        if err:
            self._mode = MicrowaveMode.CW
            power = self.get_power()
            freq = self.get_frequency()
            return freq, power, 'cw'
        else:
            self._mode = MicrowaveMode.LIST
            power = self.get_power()
            return frequency, power, 'list'

    def reset_listpos(self):
        """
        Reset of MW list mode position to start (first frequency step)

        @return int: error code (0:OK, -1:error)
        """

        if self._mode == MicrowaveMode.LIST:
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
        self.log.warning('Trigger is not implemented yet.')
        return -1

    def get_limits(self):
        """ Return the device-specific limits in a nested dictionary.

          @return MicrowaveLimits: Microwave limits object
        """
        limits = self._microwave_device.get_limits()
        limits.list_maxentries = self._awg_device.get_constraints().sequence_num.max
        limits.supported_modes = (MicrowaveMode.CW, MicrowaveMode.SWEEP, MicrowaveMode.ASWEEP, MicrowaveMode.LIST)
        return limits

    def _set_awg_list(self):
        """
        @return int: error code (0:OK, -1:error)
        """

        frequency = self._freq_list

        llMemSamples = 1024*self.number_of_samples
        llLoops = self.number_of_loops
        amplitude = self.awg_amplitude
        frequency_local_oscillator = frequency[0] - self.awg_offset
        sample_rate = self.awg_sample_rate

        self._awg_device.pulser_off()
        self._awg_device.set_active_channels(ch={'a_ch1': True, 'a_ch2': False,
                                                 'd_ch1': True, 'd_ch2': True, 'd_ch3': False})
        self._awg_device.set_analog_level(amplitude={'a_ch1': amplitude})
        self._awg_device.set_sample_rate(sample_rate)

        channel_data = []
        t = np.arange(llMemSamples) / sample_rate
        for i in range(len(frequency)):
            a_ch1_signal = (np.sin(t * (frequency[i] - frequency_local_oscillator) * 2*np.pi) * (2**15-1)
                            ).astype(dtype=np.int16)
            d_ch1_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)

            #Fixme
            if i == 0:
                d_ch2_signal = np.full(llMemSamples, 1).astype(dtype=np.bool)
            else:
                d_ch2_signal = np.full(llMemSamples, 0).astype(dtype=np.bool)

            channel_data.append({'a_ch1': a_ch1_signal, 'd_ch1': d_ch1_signal, 'd_ch2': d_ch2_signal})

        # Fixme: don't generate all the data first and then upload it all at once -> needs a lot of memory
        err = self._awg_device.load_sequence(channel_data, llLoops)
        if err:
            return -1
        else:
            return 0
