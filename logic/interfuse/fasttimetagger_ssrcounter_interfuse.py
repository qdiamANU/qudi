# -*- coding: utf-8 -*-

"""
This file contains the Qudi Interfuse between single shot logic and fastcounter hardware.

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


"""
An interfuse file is indented to fuse/combine a logic with a hardware, which
was not indented to be used with the logic. The interfuse file extend the
ability of a hardware file by converting the logic calls (from a different
interface) to the interface commands, which suits the hardware.
In order to be addressed by the (magnet) logic it should inherit the (magnet)
interface, and given the fact that it will convert a magnet logic call to a
motor hardware call, that 'interfuse' file has to stick to the interfaces
methods of the motor interface.

Reimplement each call from the magnet interface and use only the motor interface
command to talk to a xyz motor hardware and a rotational motor hardware.
"""

import numpy as np
from core.util.network import netobtain
from core.module import Connector
from logic.generic_logic import GenericLogic
from interface.single_shot_interface import SingleShotInterface


class FastTimeTaggerSSRCounterInterfuse(GenericLogic, SingleShotInterface):

    _modclass = 'FastTimeTaggerSSRCounterInterfuse'
    _modtype = 'interfuse'

    # declare connectors, here you can see the interfuse action: the in
    # connector will cope a motor hardware, that means a motor device can
    # connect to the in connector of the logic.
    #FIXME:
    #fastcounter = Connector(interface='MotorInterface')
    #pulsedmeasurementlogic = Connector(interface='MotorInterface')
    fastcounter = Connector(interface='FastCounterInterface')
    pulsedmeasurementlogic = Connector(interface='PulsedMeasurementLogic')



    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """
        self._fastcounter = self.fastcounter()
        self._pulsedmeasurementlogic = self.pulsedmeasurementlogic()

        # self.charge_state_selection = False
        # self.charge_threshold = 5
        # self.subtract_mean = True

        #todo: check if these two below can be commented out - AJH 6/2/19
        # self.counts_per_readout = 1
        # self.countlength = 1

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        pass

    def get_constraints(self):
        """ Retrieve the hardware constrains from the magnet driving device.

        @return dict: dict with constraints for the magnet hardware. These
                      constraints will be passed via the logic to the GUI so
                      that proper display elements with boundary conditions
                      could be made.
        """
        constraints = self._fastcounter.get_constraints()
        return constraints


    def configure_ssr_counter(self, counts_per_readout, countlength):
        """ Configuration of the fast counter.

        from ssr_counter interface, may be out of date:
        @param int counts_per_readout: optional, number of readouts for one measurement
        @param int cycles: countlength, number of measurements

        @return tuple(binwidth_s, record_length_s, preset, cycles, sequences ):
                    binwidth_s: float the actual set binwidth in seconds
                    gate_length_s: the actual record length in seconds
                    preset: number of readout per measurement
                    cycles: number of measurements
                    sequences: number of sequences
        """
        # self._fastcounter.configure_ssr_counter(counts_per_readout, countlength)
        # self._fastcounter.configure(bin_width_s, record_length_s, number_of_gates=0)
        self.counts_per_readout = counts_per_readout
        self.countlength = countlength

    def get_status(self):
        """ Get the status of the position

        @param list param_list: optional, if a specific status of an axis
                                is desired, then the labels of the needed
                                axis should be passed in the param_list.
                                If nothing is passed, then from each axis the
                                status is asked.

        @return dict: with the axis label as key and the status number as item.
        """
        status = self._fastcounter.get_status()
        return status

    def start_measure(self):
        """ Start the fast counter. """
        status = self._fastcounter.start_measure()
        return status


    def stop_measure(self):
        """ Stop the fast counter. """
        status = self._fastcounter.stop_measure()
        return status


    def pause_measure(self):
        """ Pauses the current measurement.

        Fast counter must be initially in the run state to make it pause.
        """
        status = self._fastcounter.pause_measure()
        return status


    def continue_measure(self):
        """ Continues the current measurement.

        If fast counter is in pause state, then fast counter will be continued.
        """
        status = self._fastcounter.continue_measure()
        return status

    def get_data_trace(self, normalized=False, charge_state_selection=False, subtract_mean=False):
        """ Polls the current timetrace data from the fastcounter.

        Return value is a numpy array (dtype = int64).
        Will return a 1D-numpy-array
        """

        raw_data = netobtain(self._fastcounter.get_data_trace())
        print('netobtain raw_data = {}, type = {}'.format(raw_data, type(raw_data)))
        # #LOCALFIX Andrew
        # raw_data = self._fastcounter.get_data_trace()
        # print('raw_data = {}, type = {}'.format(raw_data, type(raw_data)))

        # remove all zeros at the end
        # first find all rows with only zeros
        data_size = np.where(~raw_data.any(axis=1))[0]
        # if there are no empty rows at the end (meaning the last row is not empty), delete nothing
        if len(data_size)==0 or len(raw_data)-1 not in data_size:
            ssr_data = raw_data
        else:
            # find all entries which have no next lower integer (max element should be the first empty row)
            test = [i for i, x in enumerate(data_size) if x - 1 not in data_size]
            # get the first row which should be deleted
            first_row = max(test)
            # delete the empty rows
            ssr_data = np.delete(raw_data, data_size[first_row:], 0)
        #ssr_data = np.delete(raw_data, data_size, 0)
        #ssr_data = raw_data
        self.raw_data = ssr_data

        if not charge_state_selection:
            charge_signal = np.array([])
            if normalized:
                print('normalised, not self.charge_state_selection')
                return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(ssr_data)
                laser1 = return_dict['laser_counts_arr']
                laser2 = return_dict['laser_counts_arr2']
                self.laser_data = [laser1, laser2]
                if laser1.any() and laser2.any():
                    tmp_signal1, tmp_error1 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser1)
                    tmp_signal2, tmp_error2 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser2)
                    tmp_signal = (tmp_signal1 - tmp_signal2) / (tmp_signal1 + tmp_signal2)
                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])
                #
                # # FIXME: here it is important that one laser pulse is on one side respectively
                # data_width = self.raw_data.shape[1]
                # print('data_width = {}'.format(data_width))
                # raw_data1 = self.raw_data[:, :int(data_width / 2)]
                # raw_data2 = self.raw_data[:, int(data_width / 2):]
                # print('raw_data1 = {}'.format(raw_data1))
                # return_dict1 = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(raw_data1)
                # print('return_dict1 = {}'.format(return_dict1))
                # laser1 = return_dict1['laser_counts_arr']
                # print('laser1 = {}'.format(laser1))
                # return_dict2 = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(raw_data2)
                # laser2 = return_dict2['laser_counts_arr']
                # print('laser2 = {}'.format(laser1))
            else:
                print('not normalised')
                return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(self.raw_data)
                self.laser_data = return_dict['laser_counts_arr']
                # analyze pulses and get data points for signal array.
                if self.laser_data.any():
                    tmp_signal, tmp_error = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(
                        self.laser_data)
                    if subtract_mean:
                        tmp_signal = tmp_signal  - np.mean(tmp_signal)
                        charge_signal = charge_signal - np.mean(charge_signal)
                    else:
                        tmp_signal = tmp_signal #- np.mean(tmp_signal)
                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])
                    # get rid of the last point since it is measured with less redouts
        else:
            # seperate the charge state and the ssr data
            ssr_raw = ssr_data[1::2, :]
            charge_raw = ssr_data[::2, :]
            # extract the orange laser pulses
            # charge_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(charge_raw)
            # self.charge_dict=charge_dict
            # charge_signal, charge_error = self._pulsedmeasurementlogic._pulseanalyzer.\
            #    analyse_laser_pulses(charge_dict['laser_counts_arr'])
            # self.charge_signal=charge_signal
            charge_signal = np.sum(charge_raw, axis=1)
            ## find all the entries which are below the threshold
            # incorrect_charge = np.where(charge_signal < self.charge_threshold)[0]
            if normalized:
                return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(ssr_raw)
                laser1 = return_dict['laser_counts_arr']
                laser2 = return_dict['laser_counts_arr2']
                self.laser_data = [laser1, laser2]
                if laser1.any() and laser2.any():
                    tmp_signal1, tmp_error1 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser1)
                    tmp_signal2, tmp_error2 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser2)
                    tmp_signal = (tmp_signal1 - tmp_signal2) / (tmp_signal1 + tmp_signal2)
                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])
            else:
                return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(ssr_raw)
                self.laser_data = return_dict['laser_counts_arr']
                # analyze pulses and get data points for signal array.
                if self.laser_data.any():
                    tmp_signal, tmp_error = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(
                        self.laser_data)
                    if subtract_mean:
                        tmp_signal = tmp_signal  - np.mean(tmp_signal)
                        charge_signal = charge_signal - np.mean(charge_signal)
                    else:
                        tmp_signal = tmp_signal #- np.mean(tmp_signal)
                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])
            # remove all the data points with the wrong charge state
            # tmp_signal = np.delete(tmp_signal, incorrect_charge[:-1], 0)
            # get rid of the last point since it is measured with less redouts
        print('tmp_signal[:-1] = {}'.format(tmp_signal[:-1]))
        print('charge_signal = {}'.format(charge_signal))
        return tmp_signal[:-1], charge_signal






