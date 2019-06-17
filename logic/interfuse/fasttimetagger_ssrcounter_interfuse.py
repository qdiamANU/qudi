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
import time

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

        self.next_channel = 3  # todo: define externally

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

    def configure_ssr_counter(self, histogram_length_s, num_of_fastcounter_rows, bin_width_s, rollover=True):
        """ Configuration of the fast counter for SSR. Counter must be gated

        @param float histogram_length_s: length of counts histogram for one readout, in seconds.
            Should be slightly longer than readout laser pulse duration
        @param int num_of_fastcounter_rows: number of required histograms, equal to the number of
            steps in the controlled variable (e.g. MW frequency) sweep

        @return tuple(binwidth_s, gate_length_s, number_of_gates):
                    binwidth_s: float the actual set binwidth in seconds
                    gate_length_s: the actual set gate length in seconds
                    number_of_gates: the number of gated, which are accepted
                    rollover: True = keep counting until stopped. Once histogram is filled, rollover to start of
                        histogram and add to counts already there. Generally, want rollover=True.
                        Want rollover=False for e.g. just-SSR measurements
        """

        print('tt configure ssr')
        self._fastcounter.configure(bin_width_s=bin_width_s,
                                    record_length_s=histogram_length_s,
                                    number_of_gates=num_of_fastcounter_rows,
                                    next_channel=self.next_channel,
                                    rollover=rollover)

        return (bin_width_s, histogram_length_s, num_of_fastcounter_rows)

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
        self.raw_data = raw_data

        # remove all zeros at the end - useful for analysing data during acquisition, where some fastcounter
        # histogram rows may still be unfilled
        self.ssr_data = raw_data[~np.all(raw_data == 0, axis=1)]

        if not charge_state_selection:
            charge_signal = np.array([])
            if normalized: ### Normalised, no charge state selection #####

                ### Pulse extraction from raw counts trace #####
                t3 = time.time()
                # for some measurement modes, only want 2 bins (one for each laser pulse), to reduce data size
                # in such cases, regular pulse extraction and pulse analysis are unnecessary, and they also don't work
                if self.ssr_data.shape[1] <= 1:
                    self.log.error('Error in SSR counter, get data: expected ssr_data.shape[1] > 1 for normalised SSR (got {})'.format(self.ssr_data.shape[1]))
                    return -1
                elif self.ssr_data.shape[1] <= 3:
                    laser1 = self.ssr_data[:, 0]
                    laser2 = self.ssr_data[:, 1]
                    perform_pulse_analysis = False
                    return_dict = dict()
                # otherwise, do pulse extraction and analysis as per normal
                else:
                    return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(self.ssr_data)
                    laser1 = return_dict['laser_counts_arr0']
                    laser2 = return_dict['laser_counts_arr1']
                    perform_pulse_analysis = True

                self.laser_data = [laser1, laser2]

                t4 = time.time()
                ### Pulse analysis of extracted data #####
                if laser1.any() and laser2.any():
                    if perform_pulse_analysis:
                        tmp_signal1, tmp_error1 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser1)
                        tmp_signal2, tmp_error2 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser2)
                    else:
                        tmp_signal1, tmp_error1 = laser1, np.sqrt(laser1)
                        tmp_signal2, tmp_error2 = laser2, np.sqrt(laser2)

                    tmp_signal = (tmp_signal1 - tmp_signal2) / (tmp_signal1 + tmp_signal2)

                    # set nan values to zero (nan occurs when both laser1 and laser2 counts are zero)
                    nan_indices = np.isnan(tmp_signal)
                    tmp_signal[nan_indices] = 0

                else:
                    # fixme: can't call shape on list object
                    tmp_signal = np.zeros(self.laser_data.shape[0])
                print('ssr_counter_interfuse pulse analysis: time taken = {} ms'.format((time.time() - t4) / 1e-3))

            else: ### Non-normalised, no charge state selection #####

                ### Pulse extraction from raw counts trace #####
                # for some measurement modes, only want 2 bins (one for each laser pulse), to reduce data size
                # in such cases, regular pulse extraction and pulse analysis are unnecessary, and they also don't work
                if self.ssr_data.shape[1] == 0:
                    self.log.error('Error in SSR counter, get data: ssr_data.shape[1] = 0!')
                    return -1
                elif self.ssr_data.shape[1] <= 2:
                    self.laser_data = self.ssr_data[:, 0]
                    perform_pulse_analysis = False
                    return_dict = dict()
                # otherwise, do pulse extraction and analysis as per normal
                else:
                    return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(self.ssr_data)
                    self.laser_data = return_dict['laser_counts_arr0']
                    perform_pulse_analysis = True

                t4 = time.time()
                ### Pulse analysis of extracted data #####
                # analyze pulses and get data points for signal array.
                if self.laser_data.any():
                    if perform_pulse_analysis:
                        tmp_signal, tmp_error = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(
                            self.laser_data)
                    else:
                        tmp_signal, tmp_error = self.laser_data, np.sqrt(self.laser_data)

                    if subtract_mean:
                        tmp_signal = tmp_signal - np.mean(tmp_signal)
                        charge_signal = charge_signal - np.mean(charge_signal)
                    else:
                        tmp_signal = tmp_signal

                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])

        ### Including charge state normalisation #####
        #todo: not currently working at ANU
        else:
            # separate the charge state and the ssr data
            ssr_raw = self.ssr_data[1::2, :]
            charge_raw = self.ssr_data[::2, :]
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
                laser1 = return_dict['laser_counts_arr0']
                laser2 = return_dict['laser_counts_arr1']
                # print('laser1.shape = {}'.format(self.laser1.shape))
                # print('laser2.shape = {}'.format(self.laser2.shape))
                self.laser_data = [laser1, laser2]
                if laser1.any() and laser2.any():
                    tmp_signal1, tmp_error1 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser1)
                    tmp_signal2, tmp_error2 = self._pulsedmeasurementlogic._pulseanalyzer.analyse_laser_pulses(laser2)
                    tmp_signal = (tmp_signal1 - tmp_signal2) / (tmp_signal1 + tmp_signal2)
                    nan_indices = np.isnan(tmp_signal)
                    tmp_signal[nan_indices] = 0
                else:
                    tmp_signal = np.zeros(self.laser_data.shape[0])
            else:
                return_dict = self._pulsedmeasurementlogic._pulseextractor.extract_laser_pulses(ssr_raw)
                self.laser_data = return_dict['laser_counts_arr0']
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
            # get rid of the last point since it is measured with less readouts (Ulm - not checked at ANU)
            # fixme Andrew 18/3/2019: pretty sure the last data point is fine... if anything, the first
            # data point should be removed, as we don't pre-polarise the nuclear spin

        # LOCALFIX Andrew: for debugging purposes:
        self.tmp_signal = tmp_signal # LOCALFIX Andrew 12/2/19: included to enable plotting of tmp_signal from jupyter script
        self.return_dict = return_dict

        return tmp_signal[:-1], charge_signal






