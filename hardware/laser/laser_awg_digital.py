# -*- coding: utf-8 -*-
"""
This module is for controlling a laser that only accepts a digital on/off signal, using an AWG digital output channel

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

from core.module import Base
from interface.simple_laser_interface import SimpleLaserInterface
from interface.simple_laser_interface import LaserState
from interface.simple_laser_interface import ShutterState
from interface.simple_laser_interface import ControlMode
import math
import random
import time

class Laser_AWG_Digital(Base, SimpleLaserInterface):
    """
    This module is for controlling a laser that only accepts a digital on/off signal,
    using an AWG digital output channel
    """
    _modclass = 'laserawgdigital'
    _modtype = 'hardware'

    def __init__(self, **kwargs):
        """ """
        super().__init__(**kwargs)
        self.lstate = LaserState.OFF
        self.shutter = ShutterState.NOSHUTTER
        self.mode = ControlMode.POWER
        self.current_setpoint = 0
        self.power_setpoint = 0

    def on_activate(self):
        """ Activate module.
        """
        pass

    def on_deactivate(self):
        """ Deactivate module.
        """
        pass

    def get_power_range(self):
        """ Return optical power range

            @return (float, float): power range
        """
        pass

    def get_power(self):
        """ Return laser power

            @return float: Laser power in watts
        """
        pass

    def get_power_setpoint(self):
        """ Return optical power setpoint.

            @return float: power setpoint in watts
        """
        pass

    def set_power(self, power):
        """ Set power setpoint.

            @param float power: power setpoint

            @return float: actual new power setpoint
        """
        pass

    def get_current_unit(self):
        """ Get unit for laser current.

            @return str: unit
        """
        pass

    def get_current_range(self):
        """ Get laser current range.

            @return (float, float): laser current range
        """
        pass

    def get_current(self):
        """ Get current laser current

            @return float: laser current in current curent units
        """
        pass

    def get_current_setpoint(self):
        """ Get laser curent setpoint

            @return float: laser current setpoint
        """
        pass

    def set_current(self, current):
        """ Set laser current setpoint

            @prarm float current: desired laser current setpoint

            @return float: actual laser current setpoint
        """
        pass

    def allowed_control_modes(self):
        """ Get supported control modes

            @return list(): list of supported ControlMode
        """
        return [ControlMode.POWER]

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

    def on(self):
        """ Turn on laser.

            @return LaserState: actual laser state
        """
        time.sleep(1)
        self.lstate = LaserState.ON
        return self.lstate

    def off(self):
        """ Turn off laser.

            @return LaserState: actual laser state
        """
        time.sleep(1)
        self.lstate = LaserState.OFF
        return self.lstate

    def get_laser_state(self):
        """ Get laser state

            @return LaserState: actual laser state
        """
        return self.lstate

    def set_laser_state(self, state):
        """ Set laser state.

            @param LaserState state: desired laser state

            @return LaserState: actual laser state
        """
        time.sleep(1)
        self.lstate = state
        return self.lstate

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
        pass

    def get_temperatures(self):
        """ Get all available temperatures.

            @return dict: dict of temperature namce and value in degrees Celsius
        """
        pass

    def set_temperatures(self, temps):
        """ Set temperatures for lasers with tunable temperatures.

            @return {}: empty dict, dummy not a tunable laser
        """
        pass

    def get_temperature_setpoints(self):
        """ Get temperature setpoints.

            @return dict: temperature setpoints for temperature tunable lasers
        """
        pass

    def get_extra_info(self):
        """ Multiple lines of dignostic information

            @return str: much laser, very useful
        """
        return "Laser with only digital on/off\n remote control\n e.g. current-pulsed lasers from Stuttgart\n controlled by AWG digital channel"

