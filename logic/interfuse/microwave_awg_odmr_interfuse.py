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

from core.util.interfaces import InterfaceMetaclass
import time
from core.module import Connector, StatusVar
from core.module import Base, ConfigOption
from logic.generic_logic import GenericLogic
from interface.microwave_interface import MicrowaveInterface
from interface.microwave_interface import MicrowaveLimits
from interface.microwave_interface import MicrowaveMode
from interface.microwave_interface import TriggerEdge

class MicrowaveAwgInterfuseAwgTriggered(GenericLogic, MicrowaveInterface):
    """This is the Interface class to define the controls for the simple
    microwave hardware.
    """

    _modclass = 'MicrowaveAwgInterfuseAwgTriggered'
    _modtype = 'interfuse'

    # declare connectors, here you can see the interfuse action: the in
    # connector will cope a motor hardware, that means a motor device can
    # connect to the in connector of the logic.
    microwave = Connector(interface='MicrowaveInterface')
    awg = Connector(interface='PulserInterface')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # todo: check if necessary
        self._mode = MicrowaveMode.CW
        self._freq_list = []

        # self.use_ext_microwave = True
        self.ext_microwave_frequency = 2.6e9
        # self.ext_microwave_power = -30

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """

        self._microwave_device = self.microwave()
        self._awg_device = self.awg()

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        pass
    
    
    def off(self):
        """
        Switches off any microwave output.
        Must return AFTER the device is actually stopped.

        @return int: error code (0:OK, -1:error)
        """
        print('interfuse off')
        return_val_1 = self._microwave_device.off()
        return_val_2 = self._awg_device.pulser_off()

        self.output_active = False

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

        mode, mw_is_running = self._microwave_device.get_status()

        if mode == 'cw':
            if self._mode == MicrowaveMode.LIST:
                mode = 'list'
            if self._mode == MicrowaveMode.SWEEP:
                mode = 'sweep'

        # awg_status, status_dic = self._awg_device.get_status()
        # # print('awg_status = {}, status dict = {}'.format(awg_status, status_dic))
        #
        # if awg_status < 0:
        #     return mode, -1
        # else:
        #     return mode, mw_is_running
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
        print('list on')
        self._microwave_device.cw_on()
        mode, status = self.get_status()
        #print("Status", mode, status)
        #time.sleep(2)
        self._awg_device.pulser_on()
        if self._awg_device.read_out_error():
            return -1
        else:
            self._mode = MicrowaveMode.LIST
            self.output_active = True
            return 0


    def set_list(self, frequency=None, power=None):
        """
        Configures the device for list-mode and optionally sets frequencies and/or power

        @param list frequency: list of frequencies in Hz
        @param float power: MW power of the frequency list in dBm

        @return list, float, str: current frequencies in Hz, current power in dBm, current mode
        """

        print('set list')
        self._freq_list = frequency

        # local_oscillator_freq = frequency[0] - self.awg_offset_frequency
        # awg_freq_array = np.array(frequency) - local_oscillator_freq
        # # awg_freq_list = list(np.array(frequency) - local_oscillator_freq)
        #
        # # sanity checks:
        # if np.min(awg_freq_array) - local_oscillator_freq < 0:
        #     self.log.error(
        #         'LO freq ({} MHz) larger than desired minimum MW freq ({} MHz). Need to adjust LO freq.'.format(
        #             local_oscillator_freq, np.min(awg_freq_array)))
        #     return frequency, power, MicrowaveMode.LIST
        #
        # if (frequency[-1] - frequency[0]) > awg_range:
        #     self.log.error(
        #         'Sweep range ({} MHz) is larger than the awg bandwidth ({} MHz).'.format(frequency[-1] - frequency[0],
        #                                                                                  awg_range))
        #     return frequency, power, MicrowaveMode.LIST
        #
        # if (frequency[-1] - local_oscillator_freq) > awg_range:
        #     self.log.error(
        #         'Combined sweep range ({} MHz) and LO offset ({} MHz) is larger than the awg bandwidth ({} MHz).'.format(
        #             frequency[-1] - frequency[0], local_oscillator_freq, awg_range))
        #     return frequency, power, MicrowaveMode.LIST

        # set microwave source in cw mode to act as local oscillator
        self.set_cw(self.ext_microwave_frequency, power)

        name = 'odmr_cw_list'
        self._awg_device.write_triggered_cw_odmr_list_sequence(frequency, name)
        self._awg_device.load_sequence(name)

        err = self._awg_device.read_out_error()

        self._awg_device.pulser_on()

        if err:
            self.log.error('Could not set AWG in list mode. Returning to CW microwave mode.')
            self._mode = MicrowaveMode.CW
            power = self.get_power()
            frequency = self.ext_microwave_frequency
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

        # todo: check if this should restart the output or not
        if self._mode == MicrowaveMode.LIST:
            # self._awg_device.pulser_off()
            # self._awg_device.pulser_on()
            self.output_active = True
            return 0
        else:
            return -1

    def sweep_on(self):
        """ Switches on the sweep mode.

        @return int: error code (0:OK, -1:error)
        """
        self._microwave_device.cw_on()
        self._awg_device.pulser_on()

        if self._awg_device.read_out_error():
            return -1
        else:
            self._mode = MicrowaveMode.SWEEP
            self.output_active = True
            return 0

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

        # set microwave source in cw mode to act as local oscillator
        self.set_cw(self.ext_microwave_frequency, power)

        name = 'odmr_cw_sweep'
        self._awg_device.write_triggered_cw_odmr_sweep_sequence(start, stop, step, name)
        self._awg_device.load_sequence(name)

        err = self._awg_device.read_out_error()

        if err:
            self.log.error('Could not set AWG in sweep mode. Returning to CW microwave mode.')
            self._mode = MicrowaveMode.CW
            power = self.get_power()
            frequency = self.ext_microwave_frequency
            return start, stop, step, power, 'cw'
        else:
            self._mode = MicrowaveMode.SWEEP
            power = self.get_power()
            return start, stop, step, power, 'sweep'

    def reset_sweeppos(self):
        """
        Reset of MW sweep mode position to start (start frequency)

        @return int: error code (0:OK, -1:error)
        """
        # todo: check if this should restart the output or not
        if self._mode == MicrowaveMode.SWEEP:
            # self._awg_device.pulser_off()
            # self._awg_device.pulser_on()
            self.output_active = True
            return 0
        else:
            return -1

    def set_ext_trigger(self, pol, timing):
        """ Set the external trigger for the AWG with proper polarization.

        @param TriggerEdge pol: polarisation of the trigger (basically rising edge or falling edge)
        @param timing: estimated time between triggers

        @return object, float: current trigger polarity [TriggerEdge.RISING, TriggerEdge.FALLING],
            trigger timing as queried from device
        """
        return pol, timing

    def trigger(self):
        """ Trigger the next element in the list or sweep mode programmatically.

        @return int: error code (0:OK, -1:error)

        Ensure that the Frequency was set AFTER the function returns, or give
        the function at least a save waiting time corresponding to the
        frequency switching speed.
        """
        print('trigger interfuse')
        self._awg_device.force_trigger()
        if self._awg_device.read_out_error():
            return -1
        else:
            return 0

    def get_limits(self):
        """ Return the device-specific limits in a nested dictionary.

          @return MicrowaveLimits: Microwave limits object

        """
        limits = self._microwave_device.get_limits()
        limits.list_maxentries = self._awg_device.get_constraints().sequence_num.max
        limits.supported_modes = (MicrowaveMode.CW, MicrowaveMode.SWEEP, MicrowaveMode.LIST)
        return limits

