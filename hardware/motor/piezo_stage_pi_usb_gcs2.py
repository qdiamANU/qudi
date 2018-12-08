# -*- coding: utf-8 -*-

"""
This file contains the hardware control for PI piezo stages running GCS2.
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

from collections import OrderedDict

from core.module import Base, ConfigOption
from interface.motor_interface import MotorInterface


class PiezoStagePI(Base, MotorInterface):
    """unstable: Lachlan Rogers, Matt van Breugel
    This is the hardware module for communicating with PI Piezo scanning stages
    over USB (via the PI dll). It uses the General Command Set 2 (GCS2) from
    the PI documentation.
    This module has been developed for the E-725 Digital Piezo Controller,
    but probably works with any PI controller that talks GCS2 over USB.
    Example configuration:
    ```
    # piezo_pi:
    #     module.Class: 'motor.piezo_stage_pi_usb_gcs2.PiezoStagePI'
    #     axis_labels:
    #         - x
    #         - y
    #         - z
    #     x:
    #         channel: 0
    #         constraints:
    #             pos_min: 0e-6
    #             pos_max: 300e-6
    #     y:
    #         channel: 1
    #         constraints:
    #             pos_min: 0e-6
    #             pos_max: 300e-6
    #     z:
    #         channel: 2
    #         constraints:
    #             pos_min: 0e-6
    #             pos_max: 300e-6
    ```
    """

    _modclass = 'PiezoStagePI'
    _modtype = 'hardware'

    _devID = ctypes.c_int()

    _double3d = ctypes.c_double * 3  # This is creating a 3D double array object
    _double1d = ctypes.c_double * 1  # This is creating a 1D double object
    _bool1d = ctypes.c_bool * 1  # This is creating a 1D bool object

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def on_activate(self):
        """ Initialise and activate the hardware module.
            @return int error code (0:OK, -1:error)
        """

        path_dll = os.path.join(get_main_dir(),
                                'thirdparty',
                                'physik_instrumente',
                                'PI_GCS2_DLL_x64.dll'
                                )
        self._pidll = ctypes.windll.LoadLibrary(path_dll)

        # Find out what devices are connected
        emptybufferpy = ' ' * 1000  # TODO find out how long this should be?
        charBuffer = ctypes.c_char_p(emptybufferpy.encode())
        bufSize = ctypes.c_int(1000)

        numofdevs = self._pidll.PI_EnumerateUSB(charBuffer, bufSize, ctypes.c_char_p(b''))

        # read the device list out of ctype charBuffer into a regular python list of strings
        device_list = charBuffer.value.decode().split('\n')

        # split list into elements, check for PI devices
        pi_devices = [device for device in device_list if 'PI' in device]
        if len(pi_devices) == 1:
            device_name = ctypes.c_char_p(pi_devices[0].encode())
            self._pidll.PI_ConnectUSB(device_name)
            self._devID = ctypes.c_int(0)

        elif len(pi_devices) > 1:
            self.log.warning('There is more than 1 PI device connected, I do not know which one to choose!')

        else:
            self.log.warning('I cannot find any connected devices with "PI" in their name.')

        if not self._pidll.PI_IsConnected(self._devID):
            return 1
        else:
            self._set_servo_state(True)
            self._configured_constraints = self.get_constraints()
            #PIdevice.SVO(axis, switchOn);
            #self._pidll.PI_MOV(self._devID, ctypes.c_char_p('0'.encode()), self._double1d(0 * 1e3))
            return 0

    def on_deactivate(self):
        """ Deinitialise and deactivate the hardware module.
            @return: error code (0:OK, -1:error)
        """

        self._set_servo_state(False)
        self._pidll.PI_RTO(self._devID, ctypes.c_char_p(''.encode()))
        return 0

    def get_constraints(self):
        """ Retrieve the hardware constrains from the motor device.
            Provides all the constraints for the xyz stage  and rot stage (like total
            movement, velocity, ...)
            Each constraint is a tuple of the form
                (min_value, max_value, stepsize)

            @return dict: dict with constraints for the sequence generation and GUI
        """

        constraints = OrderedDict()

        config = self.getConfiguration()

        axis0 = {}
        axis0['label'] = 'x'
        axis0['channel'] = config['x']['channel']
        axis0['pos_min'] = config['x']['constraints']['pos_min']
        axis0['pos_max'] = config['x']['constraints']['pos_max']
        axis0['unit'] = 'm'
        axis0['pos_step'] = 10e-6
        axis0['vel_min'] = 0.0
        axis0['vel_max'] = 100.0
        axis0['vel_step'] = 0.01
        axis0['acc_min'] = 0.001
        axis0['acc_max'] = 0.001
        axis0['acc_step'] = 0.001

        axis1 = {}
        axis1['label'] = 'y'
        axis1['channel'] = config['y']['channel']
        axis1['pos_min'] = config['y']['constraints']['pos_min']
        axis1['pos_max'] = config['y']['constraints']['pos_max']
        axis1['unit'] = 'm'
        axis1['pos_step'] = 10e-6
        axis1['vel_min'] = 1.0
        axis1['vel_max'] = 20.0
        axis1['vel_step'] = 0.1
        axis1['acc_min'] = 0.0001
        axis1['acc_max'] = 0.0001
        axis1['acc_step'] = 0.001

        axis2 = {}
        axis2['label'] = 'z'
        axis2['channel'] = config['z']['channel']
        axis2['pos_min'] = config['z']['constraints']['pos_min']
        axis2['pos_max'] = config['z']['constraints']['pos_max']
        axis2['unit'] = 'm'
        axis2['pos_step'] = 10e-6
        axis2['vel_min'] = 1.0
        axis2['vel_max'] = 20.0
        axis2['vel_step'] = 0.1
        axis2['acc_min'] = 0.0001
        axis2['acc_max'] = 0.0001
        axis2['acc_step'] = 0.001

        # assign the parameter container for x to a name which will identify it
        constraints[axis0['label']] = axis0
        constraints[axis1['label']] = axis1
        constraints[axis2['label']] = axis2

        return constraints

    def move_rel(self, param_dict):
        """ Move the stage relatively in given direction
            @param dict param_dict : dictionary, which passes all the relevant
                                     parameters, which should be changed. Usage:
                                     {'axis_label': <the-abs-pos-value>}.
                                     'axis_label' must correspond to a label given
                                     to one of the axis.
            @return dict param_dict : dictionary with the current stage positions
        """

        invalid_axis = set(param_dict) - set(['x', 'y', 'z'])
        if invalid_axis:
            for axis in invalid_axis:
                self.log.warning('Desired axis {axis} is undefined'
                                 .format(axis=axis))

        self._configured_constraints = self.get_constraints()

        for axis in ['x', 'y', 'z']:
            if axis in param_dict.keys():
                if axis == 'x':
                    channel = self._configured_constraints[axis]['channel']
                    rel_move = param_dict[axis]
                    current_position = self.get_pos()
                    print(current_position[axis])
                    to_position = rel_move + current_position[axis]
                    self._do_move_abs(axis, channel, to_position)
                elif axis == 'y':
                    channel = self._configured_constraints[axis]['channel']
                    rel_move = param_dict[axis]
                    current_position = self.get_pos()
                    print(current_position[axis])
                    to_position = rel_move + current_position[axis]
                    self._do_move_abs(axis, channel, to_position)
                elif axis == 'z':
                    channel = self._configured_constraints[axis]['channel']
                    rel_move = param_dict[axis]
                    current_position = self.get_pos()
                    print(current_position[axis])
                    to_position = rel_move + current_position[axis]
                    self._do_move_abs(axis, channel, to_position)

        param_dict = self.get_pos()
        return param_dict


    def move_abs(self, param_dict):
        """ Move the stage to an absolute position
        @param dict param_dict : dictionary, which passes all the relevant
                                 parameters, which should be changed. Usage:
                                 {'axis_label': <the-abs-pos-value>}.
                                 'axis_label' must correspond to a label given
                                 to one of the axis.
                                 The values for the axes are in meter,
                                 the value for the rotation is in degrees.
        @return dict param_dict : dictionary with the current axis positions
        """

        invalid_axis = set(param_dict) - set(['x', 'y', 'z'])
        if invalid_axis:
            for axis in invalid_axis:
                self.log.warning('Desired axis {axis} is undefined'
                                 .format(axis=axis))

        self._configured_constraints = self.get_constraints()

        for axis in ['x', 'y', 'z']:
            if axis in param_dict.keys():
                if axis == 'x':
                    channel = self._configured_constraints[axis]['channel']
                    to_position = param_dict['x']
                    self._do_move_abs(axis, channel, to_position)
                elif axis == 'y':
                    channel = self._configured_constraints[axis]['channel']
                    to_position = param_dict['y']
                    self._do_move_abs(axis, channel, to_position)
                elif axis == 'z':
                    channel = self._configured_constraints[axis]['channel']
                    to_position = param_dict['z']
                    self._do_move_abs(axis, channel, to_position)

        param_dict = self.get_pos()
        return param_dict

    def abort(self):
        """ Stop movement of the stage

        BOOL PI_StopAll (int ID)
        Corresponding command: #24
        Stops the motion of all axes instantaneously. Sets error code to 10.
        Arguments:
        ID ID of controller
        Returns:
        TRUE if successful, FALSE otherwise

        @return int: error code (0:OK, -1:error)
        """
        print('magnet abort')

        successful = self._pidll.PI_StopAll(self._devID)

        if successful:
            return 0
        else:
            return -1

    def get_pos(self, param_list=None):
        """ Get rhe current position of the stage axis
        @param list param_list : optional, if a specific position of an axis
                                 is desired, then the labels of the needed
                                 axis should be passed in the param_list.
                                 If nothing is passed, then the positions of
                                 all axes are returned.
        @return dict param_dict : with keys being the axis labels and item the current
                                  position.
        """
        # Now we create an instance of that object
        posBuffer = self._double3d()
        axesBuffer = ctypes.c_char_p(''.encode())

        err = self._pidll.PI_qPOS(ctypes.c_int(0), axesBuffer, posBuffer)
        param_dict = {}
        param_dict['x'] = posBuffer[0] / 1e3  # unit conversion from communication
        param_dict['y'] = posBuffer[1] / 1e3  # unit conversion from communication
        param_dict['z'] = posBuffer[2] / 1e3  # unit conversion from communication

        if param_list:
            param_list = [x.lower() for x in param_list]  # make all param_list elements lower case
            for axis in list(set(param_dict.keys()) - set(param_list)):  # axes not in param_list
                del param_dict[axis]
                # print(param_dict)
            return param_dict
        else:
            # print(param_dict)
            return param_dict

    def get_status(self, param_list=None):
        """ Get the status of the position
            @param list param_list : optional, if a specific status of an axis
                                     is desired, then the labels of the needed
                                     axis should be passed in the param_list.
                                     If nothing is passed, then from each axis the
                                     status is asked.
            @return dict : with the axis label as key and the status number as item.

            The meaning of the return value is:
            Bit 0: Ready Bit 1: On target Bit 2: Reference drive active Bit 3: Joystick ON
            Bit 4: Macro running Bit 5: Motor OFF Bit 6: Brake ON Bit 7: Drive current active
        """
        self.log.info('Not yet implemented for this hardware')

        # BOOL PI_IsMoving(int ID, const char * szAxes, BOOL * pbValueArray)
        # BOOL PI_IsControllerReady(int ID, int * piControllerReady)
        # BOOL PI_IsRunningMacro(int ID, BOOL * pbRunningMacro)

    def calibrate(self, param_list=None):
        """ Calibrate the stage.
            @param dict param_list : param_list: optional, if a specific calibration
                                     of an axis is desired, then the labels of the
                                     needed axis should be passed in the param_list.
                                     If nothing is passed, then all connected axis
                                     will be calibrated.
            After calibration the stage moves to home position which will be the
            zero point for the passed axis.
            @return dict pos : dictionary with the current position of the axis
        """
        self.log.info('Not yet implemented for this hardware')

        pos = {}

        return pos

    def get_velocity(self, param_list=None):
        """ Get the current velocity for all connected axes in m/s.
            @param list param_list : optional, if a specific velocity of an axis
                                     is desired, then the labels of the needed
                                     axis should be passed as the param_list.
                                     If nothing is passed, then from each axis the
                                     velocity is asked.
            @return dict : with the axis label as key and the velocity as item.
        """
        self.log.info('Function not yet implemented for this stage')

    def set_velocity(self, param_dict):
        """ Write new value for velocity in m/s.
        @param dict param_dict : dictionary, which passes all the relevant
                                 parameters, which should be changed. Usage:
                                 {'axis_label': <the-velocity-value>}.
                                 'axis_label' must correspond to a label given
                                 to one of the axis.
        @return dict param_dict2 : dictionary with the updated axis velocity
        """
        self.log.info('Not yet implemented for this hardware')

    ########################## internal methods ##################################

    def _do_move_abs(self, axis, channel, to_pos):
        """ Make absolute axis move in meters  (internal method)
            @param str axis     : name of the axis that should be moved
                   int channel  : channel of the axis to be moved
                   float to_pos : desired position in meters
        """

        if not (self._configured_constraints[axis]['pos_min'] <= to_pos <= self._configured_constraints[axis][
            'pos_max']):
            self.log.warning('Cannot make the movement of the axis "{axis}"'
                             'since the border [{min},{max}] would be crossed! Ignore command!'
                             ''.format(axis=axis, min=self._configured_constraints[axis]['pos_min'],
                                       max=self._configured_constraints[axis]['pos_max']))
        else:
            self._write_axis_move(axis, channel, to_pos)

    def _write_axis_move(self, axis, channel, to_pos):
        """ Move a specified axis (internal method)
        @param axis string: name of the axis that should be moved
        @param int channel: channel of the axis to be moved
        @param float to_pos: desired position in meters
        """

        newpos = self._double1d(to_pos * 1e3)  # unit conversion for communication
        ax = ctypes.c_char_p(channel.encode())
        #ax = ctypes.c_char_p(channel)

        # send move command:
        self._pidll.PI_MOV(self._devID, ax, newpos)
        # check if stage has reached its target (PI_qONT)
        # todo: this should be removed, and replaced by a checkifmoving call higher up
        # as-is, the wait-til-there command prevents abort commands
        onT = self._bool1d(0)
        while not onT[0]:
            self._pidll.PI_qONT(self._devID, ax, onT)

    def _set_servo_state(self, to_state):
        """ Set the servo state (internal method)
            @param bool to_state : desired state of the servos
        """

        servo_state = self._bool1d()

        axis_list = ['1', '2', '3']

        for axis in axis_list:
            axesBuffer = ctypes.c_char_p(str(axis).encode())

            self._pidll.PI_qSVO(self._devID, axesBuffer, servo_state)

            if (servo_state[0] is False) and (to_state is True):
                self._pidll.PI_SVO(self._devID, axis, self._bool1d(1))
            elif (servo_state[0] is True) and (to_state is False):
                self._pidll.PI_SVO(self._devID, axis, self._bool1d(0))