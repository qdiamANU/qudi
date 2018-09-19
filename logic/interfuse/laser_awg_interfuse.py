# -*- coding: utf-8 -*-

"""
This file contains the Qudi Interfuse between Magnet Logic and Motor Hardware.

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

from core.module import Connector
from logic.generic_logic import GenericLogic
from interface.simple_laser_interface import SimpleLaserInterface
from interface.simple_laser_interface import LaserState
from interface.simple_laser_interface import ShutterState
from interface.simple_laser_interface import ControlMode


class LaserAWGInterfuse(GenericLogic, SimpleLaserInterface):

    _modclass = 'LaserAWGInterfuse'
    _modtype = 'interfuse'

    # declare connectors, here you can see the interfuse action: the in
    # connector will cope a motor hardware, that means a motor device can
    # connect to the in connector of the logic.
    # laser = Connector(interface='SimpleLaserInterface')
    awg = Connector(interface='PulserInterface')


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.LaserState = LaserState.OFF
        self.shutter = ShutterState.NOSHUTTER
        self.mode = ControlMode.POWER
        self.current_setpoint = 0
        self.power_setpoint = 0

    def on_activate(self):
        """ Activate module.
        """
        self._awg_device = self.awg()


    def on_deactivate(self):
        """ Deactivate module.
        """
        pass

    def get_power_range(self):
        """ Return optical power range

            @return (float, float): power range
        """
        return (0, 0.250)

    def get_power(self):
        """ Return laser power

            @return float: Laser power in watts
        """
        return 0

    def get_power_setpoint(self):
        """ Return optical power setpoint.

            @return float: power setpoint in watts
        """
        return self.power_setpoint

    def set_power(self, power):
        """ Set power setpoint.

            @param float power: power setpoint

            @return float: actual new power setpoint
        """
        self.power_setpoint = power
        self.current_setpoint = math.sqrt(4*self.power_setpoint)*100
        return self.power_setpoint

    def get_current_unit(self):
        """ Get unit for laser current.

            @return str: unit
        """
        return '%'

    def get_current_range(self):
        """ Get laser current range.

            @return (float, float): laser current range
        """
        return (0, 100)

    def get_current(self):
        """ Get current laser current

            @return float: laser current in current curent units
        """
        return 0

    def get_current_setpoint(self):
        """ Get laser curent setpoint

            @return float: laser current setpoint
        """
        return self.current_setpoint

    def set_current(self, current):
        """ Set laser current setpoint

            @prarm float current: desired laser current setpoint

            @return float: actual laser current setpoint
        """
        self.current_setpoint = current
        self.power_setpoint = math.pow(self.current_setpoint/100, 2) / 4
        return self.current_setpoint

    def allowed_control_modes(self):
        """ Get supported control modes

            @return list(): list of supported ControlMode
        """
        return [ControlMode.POWER, ControlMode.CURRENT]

    def get_control_mode(self):
        """ Get the currently active control mode

            @return ControlMode: active control mode
        """
        return self.mode

    def set_control_mode(self, control_mode):
        """ Set the active control mode

            @param ControlMode control_mode: desired control mode

            @return ControlMode: actual active ControlMode
        """
        self.mode = control_mode
        return self.mode

    def get_laser_state(self):
        """ Get laser operation state

        @return LaserState: laser state
        """
        # if self.psu == PSUTypes.SMD6000:
        #     state = self.inst.query('STAT?')
        # else:
        #     state = self.inst.query('STATUS?')S
        # if 'ENABLED' in state:
        #     return LaserState.ON
        # elif 'DISABLED' in state:
        #     return LaserState.OFF
        # else:
        #     return LaserState.UNKNOWN

        # check AWG output state, convert into on/off

        return self.LaserState

    def set_laser_state(self, status):
        """ Set desired laser state.

        @param LaserState status: desired laser state
        @return LaserState: actual laser state
        """
        actstat = self.get_laser_state()
        if actstat != status:
            if status == LaserState.ON:
                # set AWG output high
                # self.awg.pulser_on
                self.LaserState = LaserState.ON
            elif status == LaserState.OFF:
                # self.inst.query('OFF')
                # set AWG output low
                # self.awg.pulser_off
                self.LaserState = LaserState.OFF
        return self.get_laser_state()

    def on(self):
        """ Turn laser on.

            @return LaserState: actual laser state
        """
        return self.set_laser_state(LaserState.ON)

    def off(self):
        """ Turn laser off.

            @return LaserState: actual laser state
        """
        return self.set_laser_state(LaserState.OFF)

    def get_shutter_state(self):
        """ Get laser shutter state

            @return ShutterState: actual laser shutter state
        """
        return self.shutter

    def set_shutter_state(self, state):
        """ Set laser shutter state.

            @param ShutterState state: desired laser shutter state

            @return ShutterState: actual laser shutter state
        """
        time.sleep(1)
        self.shutter = state
        return self.shutter

    def get_temperatures(self):
        """ Get all available temperatures.

            @return dict: dict of temperature namce and value in degrees Celsius
        """
        return {
            'psu': 0,
            'head': 0
            }

    def set_temperatures(self, temps):
        """ Set temperatures for lasers with tunable temperatures.

            @return {}: empty dict, dummy not a tunable laser
        """
        return {}

    def get_temperature_setpoints(self):
        """ Get temperature setpoints.

            @return dict: temperature setpoints for temperature tunable lasers
        """
        return {'psu': 32.2, 'head': 42.0}

    def get_extra_info(self):
        """ Multiple lines of dignostic information

            @return str: much laser, very useful
        """
        return "Laser with only digital on/off remote control, e.g. current-pulsed lasers from Stuttgart\ncontrolled by AWG digital channel"



