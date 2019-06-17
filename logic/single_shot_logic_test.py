# -*- coding: utf-8 -*-
"""
This file contains the Qudi counter logic class.

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

import datetime
import numpy as np
from core.util.mutex import Mutex
import time

from core.module import Connector, ConfigOption, StatusVar
from collections import OrderedDict
from core.module import Connector
from core.util import units
from core.util.network import netobtain
from logic.generic_logic import GenericLogic
from qtpy import QtCore
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import copy



class SingleShotLogic(GenericLogic):
    """ This class brings raw data coming from ssrcounter measurements (gated or ungated)
        into trace form processable by the trace_analysis_logic.
    """

    _modclass = 'SingleShotLogic'
    _modtype = 'logic'

    # declare connectors
    savelogic = Connector(interface='SaveLogic')
    fitlogic = Connector(interface='FitLogic')
    singleshotreadoutcounter = Connector(interface='SingleShotInterface')
    pulsedmeasurementlogic = Connector(interface='PulsedMeasurementLogic')
    pulsedmasterlogic = Connector(interface='PulsedMasterLogic')

    # ssr counter settings
    countlength = StatusVar(default=100)
    counts_per_readout = StatusVar(default=150)
    num_bins  = StatusVar(default=20)
    init_threshold0 = StatusVar(default=0)
    init_threshold1 = StatusVar(default=0)
    ana_threshold0 = StatusVar(default=0)
    ana_threshold1 = StatusVar(default=0)
    analyze_mode = StatusVar(default='full')
    ssr_mode = StatusVar(default='flip_probability')
    sequence_length = StatusVar(default=1)
    analysis_period = StatusVar(default=5)
    normalized = StatusVar(default=False)
    charge_state_selection = StatusVar(default=False)
    subtract_mean = StatusVar(default=True)


    # measurement timer settings
    timer_interval = StatusVar(default=5)


    # ssr measurement settings
    _number_of_ssr_readouts = StatusVar(default=1000)
    _controlled_variable = StatusVar(default=list(range(50)))
    # _normalized = StatusVar(default=False)
    # _charge_state_selection  = StatusVar(default=False)

    # Container to store measurement information about the currently loaded sequence
    _ssr_measurement_information = StatusVar(default=dict())

    # notification signals for master module (i.e. GUI)
    sigssrMeasurementDataUpdated = QtCore.Signal()
    sigTimerUpdated = QtCore.Signal(float, int, float)
    sigFitUpdated = QtCore.Signal(str, np.ndarray, object)
    sigMeasurementStatusUpdated = QtCore.Signal(bool, bool)
    sigPulserRunningUpdated = QtCore.Signal(bool)
    sigExtMicrowaveRunningUpdated = QtCore.Signal(bool)
    sigExtMicrowaveSettingsUpdated = QtCore.Signal(dict)
    sigFastCounterSettingsUpdated = QtCore.Signal(dict)
    sigssrMeasurementSettingsUpdated = QtCore.Signal(dict)
    sigAnalysisSettingsUpdated = QtCore.Signal(dict)
    sigExtractionSettingsUpdated = QtCore.Signal(dict)
    # Internal signals
    sigStartTimer = QtCore.Signal()
    sigStopTimer = QtCore.Signal()

    sigStatusSSRUpdated = QtCore.Signal(bool)
    sigSSRCounterSettingsUpdated = QtCore.Signal(dict)
    sigNumBinsUpdated = QtCore.Signal(int)
    sigSequenceLengthUpdated = QtCore.Signal(float)
    sigAnalysisPeriodUpdated = QtCore.Signal(float)
    sigNormalizedUpdated = QtCore.Signal(bool)
    sigChargeStateSelectionUpdated = QtCore.Signal(bool)
    sigSubtractMeanUpdated = QtCore.Signal(bool)
    sigAnalyzeModeUpdated = QtCore.Signal(str)
    sigSSRModeUpdated = QtCore.Signal(str)
    sigThresholdUpdated = QtCore.Signal(dict)
    sigTraceUpdated = QtCore.Signal(np.ndarray, np.ndarray, float, float, float, int, float)
    sigHistogramUpdated = QtCore.Signal(np.ndarray)
    sigFitUpdated = QtCore.Signal(dict)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self.log.debug('The following configuration was found.')
        # checking for the right configuration
        for key in config.keys():
            self.log.debug('{0}: {1}'.format(key, config[key]))

        # timer for measurement
        self.__analysis_timer = None
        self.__start_time = 0
        self.__elapsed_time = 0
        self.__elapsed_sweeps = 0 # not used yet, could be used for length of timetrace

        # threading
        self._threadlock = Mutex()

        # measurement data
        self.time_axis = np.zeros(25)
        self.trace = np.zeros(25)
        self.charge_threshold = 0
        self.hist_data = list([np.linspace(1,25,25),np.ones(24)])
        self.laser_data = np.zeros((10, 20), dtype='int64')
        self.raw_data = np.zeros((10, 20), dtype='int64')
        self.spin_flip_prob = 0
        self.spin_flip_error = 0
        self.mapped_state = 0
        self.lost_events = 0
        self.lost_events_percent = 0
        self.threshold_optimal = 0
        self.fidelity_left = 0
        self.fidelity_right = 0
        self.fidelity_total = 0
        self.fidelity_left_optimal = 0
        self.fidelity_right_optimal = 0
        self.fidelity_total_optimal = 0

        self.ssr_counter_rollover = False


        self._saved_raw_data = OrderedDict()  # temporary saved raw data
        self._recalled_raw_data_tag = None  # the currently recalled raw data dict key
        self._recalled_charge_data_tag = None  # the currently recalled raw data dict key

        # Paused measurement flag
        self.__is_paused = False


        # for fit:
        self.fc = None  # Fit container
        self.signal_fit_data = np.empty((2, 0), dtype=float)  # The x,y data of the fit result


        # min trace length for fitting histogram
        self.min_trace_length = 1

        return

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """


        # QTimer must be created here instead of __init__ because otherwise the timer will not run
        # in this logic's thread but in the manager instead.
        self.__analysis_timer = QtCore.QTimer()
        self.__analysis_timer.setSingleShot(False)
        self.__analysis_timer.setInterval(round(1000. * self.timer_interval))
        self.__analysis_timer.timeout.connect(self._ssr_analysis_loop,
                                              QtCore.Qt.QueuedConnection)

        # Fitting
        self.fc = self.fitlogic().make_fit_container('pulsed', '1d')
        self.fc.set_units(['s', 'arb.u.'])

        # Recall saved status variables
        if 'fits' in self._statusVariables and isinstance(self._statusVariables.get('fits'), dict):
            self.fc.load_from_dict(self._statusVariables['fits'])

        # Turn off pulse generator
        self.pulsedmeasurementlogic().pulse_generator_off()

        # Convert controlled variable list into numpy.ndarray
        self._controlled_variable = np.array(self._controlled_variable, dtype=float)


        # update gui
        # For as long as there are no status variables define here:
        self.set_number_of_histogram_bins(self.num_bins)
        self.set_threshold({'init_threshold0': self.init_threshold0,
                            'init_threshold1': self.init_threshold1,
                            'ana_threshold0': self.ana_threshold0,
                            'ana_threshold1': self.ana_threshold1})
        # self.set_ssr_counter_settings({'counts_per_readout': self.counts_per_readout,
        #                                'countlength': self.countlength})
        # self.set_ssr_counter_settings({'num_of_points': qm_dict['num_of_points'],
        #                                'countlength': qm_dict['countlength']})



        # recalled saved raw data dict key
        self._recalled_raw_data_tag = None
        self._recalled_charge_data_tag = None

        # Connect internal signals
        self.sigStartTimer.connect(self.__analysis_timer.start, QtCore.Qt.QueuedConnection)
        self.sigStopTimer.connect(self.__analysis_timer.stop, QtCore.Qt.QueuedConnection)
        return


    def on_deactivate(self):
        """ Deactivate the module properly.
        """
        if self.module_state() == 'locked':
            self.stop_ssr_measurement()

        self._statusVariables['_controlled_variable'] = list(self._controlled_variable)
        if len(self.fc.fit_list) > 0:
            self._statusVariables['fits'] = self.fc.save_to_dict()


        self.__analysis_timer.timeout.disconnect()
        self.sigStartTimer.disconnect()
        self.sigStopTimer.disconnect()
        return

        ############################################################################
        # ssr counter control methods and properties
        ############################################################################
    @property
    def ssr_counter_constraints(self):
        return self.singleshotreadoutcounter().get_constraints()


    def ssr_counter_on(self):
        """Switching on the ssr counter

        @return int: error code (0:OK, -1:error)
        """
        return self.singleshotreadoutcounter().start_measure()

    def ssr_counter_off(self):
        """Switching off the ssr counter

        @return int: error code (0:OK, -1:error)
        """
        return self.singleshotreadoutcounter().stop_measure()

    @QtCore.Slot(bool)
    def toggle_ssr_counter(self, switch_on):
        """
        """
        if not isinstance(switch_on, bool):
            return -1

        if switch_on:
            err = self.ssr_counter_on()
        else:
            err = self.ssr_counter_off()
        return err





    ############################################################################
    # Measurement control methods and properties
    ############################################################################


    def toggle_ssr_measurement(self, start, stash_raw_data_tag=''):
        """
        Convenience method to start/stop measurement

        @param bool start: Start the measurement (True) or stop the measurement (False)
        """
        if start:
            self.start_ssr_measurement(stash_raw_data_tag)
        else:
            self.stop_ssr_measurement(stash_raw_data_tag)
        return


    def start_ssr_measurement(self, stashed_raw_data_tag=''):
        """Start the ssr measurement."""
        with self._threadlock:
            if self.module_state() == 'idle':
                # Lock module state
                self.module_state.lock()

                self.pulsedmeasurementlogic().start_simple_pulsed_measurement(stashed_raw_data_tag)
                # initialize analysis_timer
                self.__elapsed_time = 0.0
                self.sigTimerUpdated.emit(self.__elapsed_time,
                                          self.__elapsed_sweeps,
                                          self.timer_interval)

                # Set starting time and start timer
                self.__start_time = time.time()
                self.sigStartTimer.emit()
                self.sigStatusSSRUpdated.emit(True)
        return


    def stop_ssr_measurement(self, stash_raw_data_tag=''):
        """
        Stop the measurement
        """
        try:
            # do one last analysis
            self._ssr_analysis_loop()
        except:
            pass
        with self._threadlock:
            if self.module_state() == 'locked':
                self.pulsedmeasurementlogic().stop_simple_pulsed_measurement(stash_raw_data_tag)
                self.sigStopTimer.emit()
                self.module_state.unlock()
                # update status
                self.sigStatusSSRUpdated.emit(False)
        return


    ############################### set & get ########################

    def set_parameters(self, settings_dict):
        # Set parameters if present
        if self.module_state() == 'idle':
            if 'num_of_points' in settings_dict:
                self.num_of_points = settings_dict['num_of_points']
            if 'analyze_mode' in settings_dict:
                self.set_analyze_mode(settings_dict['analyze_mode'])
            if 'ssr_mode' in settings_dict:
                self.set_ssr_mode(settings_dict['ssr_mode'])
            # if 'num_bins' in settings_dict:
            #     self.set_number_of_histogram_bins(settings_dict['num_bins'])
            if 'sequence_length' in settings_dict:
                self.set_sequence_length(settings_dict['sequence_length'])
            if 'analysis_period' in settings_dict:
                self.set_analysis_period(settings_dict['analysis_period'])
            if 'ssr_normalise' in settings_dict:
                self.set_normalized(bool(settings_dict['ssr_normalise']))
            if 'threshold_dict' in settings_dict:
                self.set_threshold(settings_dict['threshold_dict'])
            if 'charge_state_selection' in settings_dict:
                self.set_charge_state_selection(settings_dict['charge_state_selection'])
            if 'charge_threshold' in settings_dict:
                self.set_charge_threshold(settings_dict['charge_threshold'])
            if 'subtract_mean' in settings_dict:
                self.set_subtract_mean(bool(settings_dict['subtract_mean']))
            if 'x_data' in settings_dict:
                self.x_data = settings_dict['x_data']
            else:
                self.x_data = None

        else:
            self.log.warning('SSR measurement is running. CAnnot change parameters')
        return



    def set_ssr_counter_settings(self, settings_dict=None, **kwargs):
        """
        Either accept a settings dictionary as positional argument or keyword arguments.
        If both are present both are being used by updating the settings_dict with kwargs.
        The keyword arguments take precedence over the items in settings_dict if there are
        conflicting names.

        @param settings_dict:
        @param kwargs:
        @return:
        """
        # Check if ssr counter is running and do nothing if that is the case
        counter_status = self.singleshotreadoutcounter().get_status()
        if not counter_status >= 2 and not counter_status < 0:
            # Determine complete settings dictionary
            if not isinstance(settings_dict, dict):
                settings_dict = kwargs
            else:
                settings_dict.update(kwargs)

            # Set parameters if present
            if 'counts_per_readout' in settings_dict:
                self.counts_per_readout = int(settings_dict['counts_per_readout'])

            if 'num_of_fastcounter_rows' in settings_dict:
                # number of fastcounter histogram rows, corresponding to e.g. number of steps of controlled variable,
                self.num_of_fastcounter_rows = settings_dict['num_of_fastcounter_rows']
            if 'countlength' in settings_dict:
                # fast counter histogram length, roughly laser pulse duration
                self.histogram_length_s = settings_dict['countlength']
            if 'bin_width_s' in settings_dict:
                # fast counter histogram length, roughly laser pulse duration
                self.bin_width_s = settings_dict['bin_width_s']
            if 'ssr_counter_rollover' in settings_dict:
                # fast counter histogram length, roughly laser pulse duration
                self.ssr_counter_rollover = settings_dict['ssr_counter_rollover']

            # Apply the settings to hardware
            self.singleshotreadoutcounter().configure_ssr_counter(histogram_length_s=self.histogram_length_s,
                                                                  num_of_fastcounter_rows=self.num_of_fastcounter_rows,
                                                                  bin_width_s=self.bin_width_s,
                                                                  rollover=self.ssr_counter_rollover)

        else:
            self.log.warning('Gated counter is not idle (status: {0}).\n'
                             'Unable to apply new settings.'.format(counter_status))

        # emit update signal for master (GUI or other logic module)
        self.sigSSRCounterSettingsUpdated.emit({'countlength': self.countlength,
                                                'counts_per_readout': self.counts_per_readout})
        return self.countlength, self.counts_per_readout


    def get_ssr_counter_settings(self):
        counter_dict=dict()
        counter_dict['countlength'] = self.countlength
        counter_dict['counts_per_readout'] = self.counts_per_readout
        return counter_dict

    def set_analyze_mode(self, mode):
        self.analyze_mode = mode
        self.sigAnalyzeModeUpdated.emit(mode)
        return

    def get_analyze_mode(self):
        return self.analyze_mode

    def set_ssr_mode(self, mode):
        self.ssr_mode = mode
        self.sigSSRModeUpdated.emit(mode)
        return

    def get_ssr_mode(self):
        return self.ssr_mode


    def set_number_of_histogram_bins(self, num_bins):
        self.num_bins = num_bins
        self.sigNumBinsUpdated.emit(num_bins)
        return self.num_bins

    def get_number_of_bins(self):
        return self.num_bins

    def set_sequence_length(self, sequence_length):
        self.sequence_length = sequence_length
        self.sigSequenceLengthUpdated.emit(sequence_length)
        return self.sequence_length

    def get_sequence_length(self):
        return self.sequence_length

    def set_analysis_period(self, analysis_period):
        with self._threadlock:
            self.timer_interval = analysis_period
            if self.timer_interval > 0:
                self.__analysis_timer.setInterval(int(1000. * self.timer_interval))
                if self.module_state() == 'locked':
                    self.sigStartTimer.emit()
            else:
                self.sigStopTimer.emit()

            self.sigTimerUpdated.emit(self.__elapsed_time, self.__elapsed_sweeps,
                                      self.timer_interval)
        return

        self.timer_interval = analysis_period
        self.__analysis_timer.setInterval(round(1000. * self.timer_interval))
        self.sigAnalysisPeriodUpdated.emit(analysis_period)
        return self.timer_interval

    def get_analysis_period(self):
        return self.timer_interval

    def set_normalized(self, norm):
        self.normalized = norm
        self.sigNormalizedUpdated.emit(norm)
        return

    def get_normalized(self):
        return self.normalized

    def set_charge_state_selection(self, mode):
        self.charge_state_selection = mode
        self.sigChargeStateSelectionUpdated.emit(mode)
        return

    def get_charge_state_selection(self):
        return self.charge_state_selection

    def set_charge_threshold(self, threshold):
        self.charge_threshold = threshold
        #self.sigChargeStateSelectionUpdated.emit(mode)
        return

    def get_charge_threshold(self):
        return self.charge_threshold

    def set_subtract_mean(self, sub):
        self.subtract_mean = sub
        self.singleshotreadoutcounter().subtract_mean = sub
        self.sigSubtractMeanUpdated.emit(sub)
        return

    def get_subtract_mean(self):
        return self.subtract_mean

    def set_threshold(self, threshold_dict):
        if 'init_threshold0' in threshold_dict:
            self.init_threshold0=threshold_dict['init_threshold0']
        if 'init_threshold1' in threshold_dict:
            self.init_threshold1 = threshold_dict['init_threshold1']
        if 'ana_threshold0' in threshold_dict:
            self.ana_threshold0 = threshold_dict['ana_threshold0']
        if 'ana_threshold1' in threshold_dict:
            self.ana_threshold1 = threshold_dict['ana_threshold1']
        self.sigThresholdUpdated.emit(threshold_dict)
        #self._ssr_analysis_loop()
        return

    def get_threshold(self):
        threshold_dict = dict()
        threshold_dict['init_threshold'] = list()
        threshold_dict['ana_threshold'] = list()
        threshold_dict['init_threshold0'] = self.init_threshold0
        threshold_dict['init_threshold1'] = self.init_threshold1
        threshold_dict['ana_threshold0'] = self.ana_threshold0
        threshold_dict['ana_threshold1'] = self.ana_threshold1
        return threshold_dict

#################################################### Analysis ##################################################

    def manually_pull_data(self):
        """ Analyse and display the data
        """
        self._ssr_analysis_loop()
        return

    def _ssr_analysis_loop(self):
        """ Acquires laser pulses from ssr counter,
            calculates fluorescence signal and creates plots.
        """
        start_time = time.time()
        with self._threadlock:
            print('singleshot enter threadlock: time taken = {} ms'.format((time.time() - start_time) / 1e-3))
            # Update elapsed time
            self.__elapsed_time = time.time() - self.__start_time

            # self.log.error('_ssr_analysis_loop: get_raw_data')
            # Get counter raw data (including recalled raw data from previous measurement)
            t1 = time.time()
            self.trace, self.charge_state = self._get_raw_data()
            print('singleshot get data: time taken = {} ms'.format((time.time()-t1)/1e-3))

            self.time_axis = np.arange(1, len(self.trace) + 1) * self.sequence_length

            t2 = time.time()
            # compute spin flip probabilities
            self.analyze_ssr()
            print('singleshot analyze_ssr: time taken = {} ms'.format((time.time() - t2) / 1e-3))

            t5 = time.time()
            # update the trace in GUI
            self.sigTraceUpdated.emit(self.time_axis, self.trace, self.spin_flip_prob,
                                      self.spin_flip_error, self.mapped_state, self.lost_events,
                                      self.lost_events_percent)
            print('singleshot sigTraceUpdated: time taken = {} ms'.format((time.time() - t5) / 1e-3))

            # compute and fit histogram
            t3 = time.time()
            self.calculate_histogram()
            print('singleshot calculate_histogram: time taken = {} ms'.format((time.time() - t3)/1e-3))
            # t4 = time.time()
            self.double_gaussian_fit_histogram()
            # print('singleshot double_gaussian_fit_histogram: time taken = {} ms'.format((time.time() - t4)/1e-3))
            # update trace in pulsed gui
            if self.x_data is not None:
                self.update_spin_flip_result()

        time_taken = time.time() - start_time
        print('_ssr_analysis_loop: time = {}\n'.format(time_taken))

        return


    def _get_raw_data(self):
        """
        Get the raw count data from the ssr counting hardware and perform sanity checks.
        Also add recalled raw data to the newly received data.
        :return numpy.ndarray: The count data (1D for ungated, 2D for gated counter)
        """

        start_time = time.time()
        # get raw data from ssr counter
        # netobtain: Check if something is a rpyc remote object. If so, transfer it (i.e. convert to float??)
        ssr_data, charge_data = netobtain(self.singleshotreadoutcounter().get_data_trace(self.normalized,
                                                                                         self.charge_state_selection,
                                                                                         self.subtract_mean))
        time_taken = time.time() - start_time
        print('_get_raw_data.get_data_trace: time = {} ms, ssr_data.shape = {}'.format(time_taken/1e-3, ssr_data.shape))
        # add old raw data from previous measurements if necessary
        if self._saved_raw_data.get(self._recalled_raw_data_tag) is not None:
            if not ssr_data.any():
                ssr_data = self._saved_raw_data[self._recalled_raw_data_tag]
                charge_data = self._saved_raw_data[self._recalled_charge_data_tag]
            elif self._saved_raw_data[self._recalled_raw_data_tag].shape == ssr_data.shape:
                ssr_data = self._saved_raw_data[self._recalled_raw_data_tag] + ssr_data
                charge_data = self._saved_raw_data[self._recalled_charge_data_tag] + charge_data
            else:
                pass
                #self.log.warning('Recalled raw data has not the same shape as current data.'
                 #                '\nDid NOT add recalled raw data to current time trace.')
        if not ssr_data.any():
            self.log.warning('Only zeros received from ssr counter!')
            ssr_data = np.zeros(ssr_data.shape, dtype='int64')
            charge_data = []
        return ssr_data, charge_data


    def analyze_ssr(self):
        """
        self.trace is the ssr_data obtained from the SSR-fastcounter, a 1D np.array of total counts / laser pulse
        in a SSR cycle

        init_threshold = The threshold of counts to accept that the nuclear spin has been initialised into a particular
        state
        ana_threshould = The threshold of counts to accept that the nuclear spin has been measured in a particular state
        :return:
        """

        # start_time = time.time()

        incorrect_charge = np.where(self.charge_state < self.charge_threshold)[0]

        # create binary array where 1 => initialisation in low-counts state
        init_low_array = copy.deepcopy(self.trace)
        init_low_array[init_low_array < self.init_threshold0] = 1
        init_low_array[init_low_array != 1] = 0

        # create binary array where 1 => initialisation in high-counts state
        init_high_array = copy.deepcopy(self.trace)
        init_high_array[init_high_array >= self.init_threshold1] = 1
        init_high_array[init_high_array != 1] = 0

        # create binary array where 1 => measurement in low-counts state
        ana_low_array = copy.deepcopy(self.trace)
        ana_low_array[ana_low_array < self.ana_threshold0] = 1
        ana_low_array[ana_low_array != 1] = 0

        # create binary array where 1 => measurement in high-counts state
        ana_high_array = copy.deepcopy(self.trace)
        ana_high_array[ana_high_array >= self.ana_threshold1] = 1
        ana_high_array[ana_high_array != 1] = 0

        # time_taken = time.time() - start_time
        # print('analyze ssr: init/ana high/low time = {} ms'.format(time_taken/1e-3))

        if self.ssr_mode == 'flip_probability':
            amount_analyzed_data = self.analyze_ssr_flip_probability(init_low_array,
                                                                     init_high_array,
                                                                     ana_low_array,
                                                                     ana_high_array,
                                                                     incorrect_charge)

        elif self.ssr_mode == 'mapping':
            amount_analyzed_data = self.analyze_ssr_mapping(init_low_array,
                                                            init_high_array,
                                                            ana_low_array,
                                                            ana_high_array,
                                                            incorrect_charge)

        else:
            self.log.error('Unknown analysis mode')
            return

        # the number of lost events is given by the length of the time_trace minus the number of analyzed data points
        self.lost_events = len(self.trace) - 1 - (amount_analyzed_data)
        self.lost_events_percent = self.lost_events / len(self.trace) * 100
        return
        #return self.spin_flip_prob, self.spin_flip_error, self.lost_events, self.lost_events_percent

    
    def analyze_ssr_flip_probability(self, init_low_array, init_high_array,
                                     ana_low_array, ana_high_array, incorrect_charge):
        """
        Method which calculates the histogram, the fidelity and the flip probability of a time trace.

        todo: pass data as a variable to function, rather than through self.trace
        :return:
        """

        # spin flip type 1: was initialised in low-counts state, then measured in high-counts state
        self.spin_flip_array_low_raw = init_low_array[:-1] * ana_high_array[1:]
        # spin flip type 2: was initialised in high-counts state, then measured in low-counts state
        self.spin_flip_array_high_raw = init_high_array[:-1] * ana_low_array[1:]

        # todo: need to remove data from measurements with the incorrect NV charge state
        # we have the self.trace indices where the charge state was incorrect
        # for calculating the mean spin flip probability, we can simply remove the
        # relevant indicies from self.spin_flip_array_raw
        # I'm not sure right now how to do this for parameter sweeps without messing up the reshaping / averaging

        # If performing a parameter sweep (e.g. nuclear rabi), then need to reshape array of individual spin flip
        # measurements into an array corresponding to the parameter sweep. (Note that for just-SSR measurements,
        # we set self.num_of_points = 1)
        # First, we pad the data with zeros so that we end up with length that is an integer multiple of the number
        # of steps in our parameter sweep
        padded_len = int(np.ceil(len(self.spin_flip_array_low_raw) / self.num_of_points) * self.num_of_points)
        spin_flip_array_low_pad = np.pad(self.spin_flip_array_low_raw, (0, padded_len - len(self.spin_flip_array_low_raw)), 'constant')
        spin_flip_array_high_pad = np.pad(self.spin_flip_array_high_raw, (0, padded_len - len(self.spin_flip_array_high_raw)), 'constant')
        # We can then reshape the data to match the parameter sweep dimensions (e.g. a Rabi measurement with 50
        # timesteps will result in a 1D self.spin_flip_array with a length of 50)
        self.spin_flip_array_low = np.mean(np.reshape(spin_flip_array_low_pad, (-1, self.num_of_points)), 0)
        self.spin_flip_array_high = np.mean(np.reshape(spin_flip_array_high_pad, (-1, self.num_of_points)), 0)

        # We can also calculate the mean spin flip probability over the entire trace
        n_datapoints = len(self.spin_flip_array_low_raw)
        n_flipped = np.sum(self.spin_flip_array_low_raw)
        n_unflipped = np.sum(1-self.spin_flip_array_low_raw)
        self.spin_flip_prob = np.mean(self.spin_flip_array_low_raw)
        self.spin_flip_error = np.sqrt((n_flipped * (1 - self.spin_flip_prob) ** 2 + n_unflipped * self.spin_flip_prob ** 2) \
                                       / (n_datapoints * (n_datapoints - 1)))

        self.spin_flip_array = self.spin_flip_array_low

        return n_datapoints


    def analyze_ssr_mapping(self, init_low_array, init_high_array, ana_low_array, ana_high_array, incorrect_charge):
        """

        fixme: not currently working
        need to perform analysis on full vectors, not just averaged values

        :param init_low_array:
        :param init_high_array:
        :param ana_low_array:
        :param ana_high_array:
        :param incorrect_charge:
        :return:
        """

        # # remove all the data with the incorrect charge state
        # self.dark_states = [x for x in low if x not in incorrect_charge]
        # self.bright_states = [x for x in high if x not in incorrect_charge]
        # # compute the probability to be in the bright state
        # if (len(self.bright_states) + len(self.dark_states)) != 0:
        #     self.mapped_state = len(self.bright_states) / (len(self.bright_states) + len(self.dark_states))
        # else:
        #     self.mapped_state = 0
        # # return the amount of analyzed data
        # amount_analyzed_data = len(self.bright_states) + len(self.dark_states)
        #
        # print('self.mapped_state = {}'.format(self.mapped_state))
        # return amount_analyzed_data

        self.log.error('analyze_ssr_mapping not currently implemented')
        return 0


    def calculate_histogram(self, custom_bin_arr=None):
        """ Calculate the histogram of a given trace.
        @param np.array trace: a 1D trace
        @param int num_bins: number of bins between the minimal and maximal
                             value of the trace. That must be an integer greater
                             than or equal to 1.
        @param np.array custom_bin_arr: optional, 1D array. If a specific,
                                        non-uniform binning array is desired
                                        then it can be passed to the numpy
                                        routine. Then the parameter num_bins is
                                        ignored. Otherwise a uniform binning is
                                        applied by default.
        @return: np.array: a 2D array, where first entry are the x_values and
                           second entry are the count values. The length of the
                           array is normally determined by the num_bins
                           parameter.
        Usually the bins for the histogram are taken to be equally spaced,
        ranging from the minimal to the maximal value of the input trace array.
        """

        if custom_bin_arr is not None:
            hist_y_val, hist_x_val = np.histogram(self.trace, custom_bin_arr,
                                                  density=False)
        else:

            # analyze the trace, and check whether all values are the same
            difference = self.trace.max() - self.trace.min()

            # if all values are the same, run at least the method with an zero
            # array. That will ensure at least an output:
            if np.isclose(0, difference) and self.num_bins is None:
                # numpy can handle an array of zeros
                num_bins = 50
                hist_y_val, hist_x_val = np.histogram(self.trace, self.num_bins)

            # if no number of bins are passed, then take the integer difference
            # between the counts, that will prevent strange histogram artifacts:
            elif not np.isclose(0, difference) and self.num_bins is None:
                hist_y_val, hist_x_val = np.histogram(self.trace, int(difference))

            # a histogram with self defined number of bins
            else:
                hist_y_val, hist_x_val = np.histogram(self.trace, self.num_bins)

        self.hist_data = np.array([hist_x_val, hist_y_val])

        # update the histogram in GUI
        self.sigHistogramUpdated.emit(self.hist_data)
        return self.hist_data


    def double_gaussian_fit_histogram(self):

        t1 = time.time()

        # fit histogram with a double Gaussian
        axis = self.hist_data[0][:-1] + (self.hist_data[0][1] - self.hist_data[0][0]) / 2.
        data = self.hist_data[1]

        add_params = dict()
        add_params['offset'] = {'min': 0, 'max': data.max(), 'value': 1e-15, 'vary': False}
        if axis.min() < self.init_threshold0:
            add_params['g0_center'] = {'min': axis.min(), 'max': self.init_threshold0,
                                       'value': self.init_threshold0 - (self.init_threshold0 - axis.min()) / 10}
        else: # todo: +/-1 in else cases only makes sense for normalised measurements
            add_params['g0_center'] = {'min': axis.min() - 1, 'max': self.init_threshold0,
                                       'value': self.init_threshold0 - (self.init_threshold0 - axis.min()) / 10}
        if axis.max() > self.init_threshold1:
            add_params['g1_center'] = {'min': self.init_threshold1, 'max': axis.max(),
                                       'value': self.init_threshold1 + (axis.max() - self.init_threshold1) / 10}
        else:
            add_params['g1_center'] = {'min': self.init_threshold1, 'max': axis.max() + 1,
                                       'value': self.init_threshold1 + (axis.max() - self.init_threshold1) / 10}
        add_params['g0_amplitude'] = {'min': data.max() / 10, 'max': 1.3 * data.max(), 'value': data.max()}
        add_params['g1_amplitude'] = {'min': data.max() / 10, 'max': 1.3 * data.max(), 'value': data.max()}
        add_params['g0_sigma'] = {'min': (self.init_threshold0 - axis.min()) / 10,
                                  'max': (axis.max() - axis.min()) / 2,
                                  'value': (axis.max() - axis.min()) / 4}
        add_params['g1_sigma'] = {'min': (axis.max() - self.init_threshold1) / 10,
                                  'max': (axis.max() - axis.min()) / 2,
                                  'value': (axis.max() - axis.min()) / 4}

        t2 = time.time()

        # Fit histogram of SSR trace with double Gaussian
        # Only attempt fit if SSR trace is long enough to create a reasonably Gaussian histogram
        if len(self.trace) > self.min_trace_length:
            # try fitting the data histogram with a double Gaussian:
            try:
                hist_fit_x, hist_fit_y, param_dict, fit_result = self.do_doublegaussian_fit(axis, data,
                                                                                            add_params=add_params)
                fit_params = fit_result.best_values
                fit_success = True
            except:
                fit_result = None
                fit_success = False
                self.log.error('SSR error in double gaussian fit')
        else:
            fit_result = None
            fit_success = False
            self.log.warning('SSR trace too short: no fitting performed')

        print('double_gaussian_fit_histogram: fit time = {} ms'.format((time.time() - t2) / 1e-3))

        if fit_success:
            # create two Gaussian functions from the two Gaussians in the double-Gaussian fit
            # These will be used for fidelity analysis
            center1 = fit_params['g0_center']
            center2 = fit_params['g1_center']
            std1 = fit_params['g0_sigma']
            std2 = fit_params['g1_sigma']
            self.gaussian1 = lambda x: fit_params['g0_amplitude'] * np.exp(-(x - center1) ** 2 / (2 * std1 ** 2))
            self.gaussian2 = lambda x: fit_params['g1_amplitude'] * np.exp(-(x - center2) ** 2 / (2 * std2 ** 2))
            if center1 > center2:
                gaussian = self.gaussian1
                self.gaussian1 = self.gaussian2
                self.gaussian2 = gaussian

            t4 = time.time()
            self.area_left = integrate.quad(self.gaussian1, -np.inf, np.inf)[0]
            self.area_right = integrate.quad(self.gaussian2, -np.inf, np.inf)[0]

            #### calculate the fidelity for the left and right part from the input threshold ###################

            # LOCALFIX Andrew 1/4/2019: want to scan ana_threshold without perturbing fit
            # area_left1 = integrate.quad(gaussian1, -np.inf, self.init_threshold0)
            # area_left2 = integrate.quad(gaussian2, -np.inf, self.init_threshold0)
            # area_right1 = integrate.quad(gaussian1, self.init_threshold1, np.inf)
            # area_right2 = integrate.quad(gaussian2, self.init_threshold1, np.inf)
            area_left1 = integrate.quad(self.gaussian1, -np.inf, self.ana_threshold0)
            area_left2 = integrate.quad(self.gaussian2, -np.inf, self.ana_threshold0)
            area_right1 = integrate.quad(self.gaussian1, self.ana_threshold1, np.inf)
            area_right2 = integrate.quad(self.gaussian2, self.ana_threshold1, np.inf)
            self.fidelity_left = area_left1[0] / (area_left1[0] + area_left2[0])
            self.fidelity_right = area_right2[0] / (area_right1[0] + area_right2[0])
            self.fidelity_total = (area_left1[0] + area_right2[0]) / (
                area_left1[0] + area_left2[0] + area_right1[0] + area_right2[0])

            #### calculate the optimal readout fidelity #######################
            # individual Gaussians
            hist_fit_x12 = np.linspace(min(hist_fit_x), max(hist_fit_x), len(hist_fit_x) * 10)
            hist_fit_y1 = self.gaussian1(hist_fit_x12)
            hist_fit_y2 = self.gaussian2(hist_fit_x12)

            print('double_gaussian_fit_histogram: integrate time = {} ms'.format((time.time() - t4) / 1e-3))

            # find crossover point between the two Gaussians, and so the optimal analysis threshold
            self.index_crossover_raw = np.argwhere(np.diff(np.sign(hist_fit_y1 - hist_fit_y2))).flatten()
            # Sometimes there are crossing points at the outer tails of the two Gaussians
            # we only want the central crossing point. We can identify this central crossover from
            # its y-value, which should be the largest of the possible crossovers
            if len(self.index_crossover_raw) != 0:
                # find index with highest y-value
                idx = np.argmax(hist_fit_y1[self.index_crossover_raw])
                index_crossover = self.index_crossover_raw[idx]
            else:
                index_crossover = None

            # Calculate optimal threshold and fidelities
            if index_crossover is None:
                self.threshold_optimal = 0
                self.fidelity_left_optimal = 0
                self.fidelity_right_optimal = 0
                self.fidelity_total_optimal = 0
            else:
                self.threshold_optimal = hist_fit_x12[index_crossover]
                area_left1 = integrate.quad(self.gaussian1, -np.inf, self.threshold_optimal)
                area_left2 = integrate.quad(self.gaussian2, -np.inf, self.threshold_optimal)
                area_right1 = integrate.quad(self.gaussian1, self.threshold_optimal, np.inf)
                area_right2 = integrate.quad(self.gaussian2, self.threshold_optimal, np.inf)
                self.fidelity_left_optimal = area_left1[0] / (area_left1[0] + area_left2[0])
                self.fidelity_right_optimal = area_right2[0] / (area_right1[0] + area_right2[0])
                self.fidelity_total_optimal = (area_left1[0] + area_right2[0]) / (
                        area_left1[0] + area_left2[0] + area_right1[0] + area_right2[0])

            fit_dict = dict()
            fit_dict['fit_result'] = fit_result
            fit_dict['fit_x'] = hist_fit_x
            fit_dict['fit_y'] = hist_fit_y
            fit_dict['fit_y1'] = self.gaussian1(hist_fit_x)
            fit_dict['fit_y2'] = self.gaussian2(hist_fit_x)
            fit_dict['fidelity_left'] = self.fidelity_left
            fit_dict['fidelity_right'] = self.fidelity_right
            fit_dict['fidelity_total'] = self.fidelity_total
            fit_dict['threshold_optimal'] = self.threshold_optimal
            fit_dict['init_threshold0'] = self.init_threshold0
            fit_dict['init_threshold1'] = self.init_threshold1
            fit_dict['ana_threshold0'] = self.ana_threshold0
            fit_dict['ana_threshold1'] = self.ana_threshold1
            fit_dict['fidelity_left_optimal'] = self.fidelity_left_optimal
            fit_dict['fidelity_right_optimal'] = self.fidelity_right_optimal
            fit_dict['fidelity_total_optimal'] = self.fidelity_total_optimal
            fit_dict['relative_area_left'] = self.area_left/(self.area_left+self.area_right)
            fit_dict['relative_area_right'] = self.area_right/(self.area_left+self.area_right)
            fit_dict['gaussian1_func'] = self.gaussian1
            fit_dict['gaussian2_func'] = self.gaussian2

            # # update the histogram in GUI with fit traces
            # self.sigFitUpdated.emit(fit_dict)

        # If fit unsuccessful:
        else:

            fit_dict = dict()
            fit_dict['fit_result'] = None
            fit_dict['fit_x'] = axis
            fit_dict['fit_y'] = data
            fit_dict['fit_y1'] = data*0
            fit_dict['fit_y2'] = data*0
            fit_dict['fidelity_left'] = 0
            fit_dict['fidelity_right'] = 0
            fit_dict['fidelity_total'] = 0
            fit_dict['threshold_optimal'] = self.threshold_optimal
            fit_dict['init_threshold0'] = self.init_threshold0
            fit_dict['init_threshold1'] = self.init_threshold1
            fit_dict['ana_threshold0'] = self.ana_threshold0
            fit_dict['ana_threshold1'] = self.ana_threshold1
            fit_dict['fidelity_left_optimal'] = 0
            fit_dict['fidelity_right_optimal'] = 0
            fit_dict['fidelity_total_optimal'] = 0
            fit_dict['relative_area_left'] = 0
            fit_dict['relative_area_right'] = 0
            fit_dict['gaussian1_func'] = 0
            fit_dict['gaussian2_func'] = 0

        # update the histogram in GUI with fit traces
        self.sigFitUpdated.emit(fit_dict)

        self.fit_dict = fit_dict
        print('double_gaussian_fit_histogram: total time = {} ms'.format((time.time() - t1) / 1e-3))

        return fit_dict

    def do_doublegaussian_fit(self, axis, data, add_params=None):
        model, params = self.fitlogic().make_gaussiandouble_model()

        if len(axis) < len(params):
            self.log.warning('Fit could not be performed because number of '
                            'parameters is larger than data points')
            return self.do_no_fit()

        else:
            result = self.fitlogic().make_gaussiandouble_fit(axis, data, self.fitlogic().estimate_gaussiandouble_peak,
                                                             add_params=add_params)
            self.result = result # LOCALFIX Andrew 2/4/2019: included for debugging

            # 1000 points in x axis for smooth fit data
            hist_fit_x = np.linspace(axis[0], axis[-1], 1000)
            hist_fit_y = model.eval(x=hist_fit_x, params=result.params)


            # this dict will be passed to the formatting method
            param_dict = OrderedDict()

            # create the proper param_dict with the values:
            param_dict['sigma_0'] = {'value': result.params['g0_sigma'].value,
                                     'error': result.params['g0_sigma'].stderr,
                                     'unit': 'Counts/s'}

            param_dict['FWHM_0'] = {'value': result.params['g0_fwhm'].value,
                                    'error': result.params['g0_fwhm'].stderr,
                                    'unit': 'Counts/s'}

            param_dict['Center_0'] = {'value': result.params['g0_center'].value,
                                      'error': result.params['g0_center'].stderr,
                                      'unit': 'Counts/s'}

            param_dict['Amplitude_0'] = {'value': result.params['g0_amplitude'].value,
                                         'error': result.params['g0_amplitude'].stderr,
                                         'unit': 'Occurrences'}

            param_dict['sigma_1'] = {'value': result.params['g1_sigma'].value,
                                     'error': result.params['g1_sigma'].stderr,
                                     'unit': 'Counts/s'}

            param_dict['FWHM_1'] = {'value': result.params['g1_fwhm'].value,
                                    'error': result.params['g1_fwhm'].stderr,
                                    'unit': 'Counts/s'}

            param_dict['Center_1'] = {'value': result.params['g1_center'].value,
                                      'error': result.params['g1_center'].stderr,
                                      'unit': 'Counts/s'}

            param_dict['Amplitude_1'] = {'value': result.params['g1_amplitude'].value,
                                         'error': result.params['g1_amplitude'].stderr,
                                         'unit': 'Occurrences'}

            param_dict['chi_sqr'] = {'value': result.chisqr, 'unit': ''}

            return hist_fit_x, hist_fit_y, param_dict, result



#FIXME: Manual fitting is not working yet
    def do_fit(self, fit_function):
        return
    #     """ Makes the a fit of the current fit function.
    #     @param str fit_function: name of the chosen fit function.
    #     @return tuple(x_val, y_val, fit_results):
    #                 x_val: a 1D numpy array containing the x values
    #                 y_val: a 1D numpy array containing the y values
    #                 fit_results: a string containing the information of the fit
    #                              results.
    #     You can obtain with get_fit_methods all implemented fit methods.
    #     """
    #
    #     if self.hist_data is None:
    #         hist_fit_x = []
    #         hist_fit_y = []
    #         param_dict = OrderedDict()
    #         fit_result = None
    #         return hist_fit_x, hist_fit_y, param_dict, fit_result
    #     else:
    #
    #         # self.log.debug((self.calculate_threshold(self.hist_data)))
    #
    #         # shift x axis to middle of bin
    #         axis = self.hist_data[0][:-1] + (self.hist_data[0][1] - self.hist_data[0][0]) / 2.
    #         data = self.hist_data[1]
    #
    #         if fit_function == 'No Fit':
    #             hist_fit_x, hist_fit_y, fit_param_dict, fit_result = self.do_no_fit()
    #             return hist_fit_x, hist_fit_y, fit_param_dict, fit_result
    #         elif fit_function == 'Gaussian':
    #             hist_fit_x, hist_fit_y, fit_param_dict, fit_result = self.do_gaussian_fit(axis, data)
    #             return hist_fit_x, hist_fit_y, fit_param_dict, fit_result
    #         elif fit_function == 'Double Gaussian':
    #             hist_fit_x, hist_fit_y, fit_param_dict, fit_result = self.do_doublegaussian_fit(axis, data)
    #             return hist_fit_x, hist_fit_y, fit_param_dict, fit_result
    #         elif fit_function == 'Poisson':
    #             hist_fit_x, hist_fit_y, fit_param_dict, fit_result = self.do_possonian_fit(axis, data)
    #             return hist_fit_x, hist_fit_y, fit_param_dict, fit_result
    #         elif fit_function == 'Double Poisson':
    #             hist_fit_x, hist_fit_y, fit_param_dict, fit_result = self.do_doublepossonian_fit(axis, data)
    #             return hist_fit_x, hist_fit_y, fit_param_dict, fit_result

    def update_spin_flip_result(self):
        """
        send data to pulsed measurement gui for analysis

        pulsedmeasurementlogic.signal_data: 2D array, 1st dim = x data, 2nd dim = y data
        same for measurement_error

        :param spin_flip_tmp:
        :param error_spin_flip_tmp:
        :param rep:
        :param index:
        :return:
        """

        temp = np.zeros((2, len(self.spin_flip_array)))
        temp[0, :] = self.x_data
        temp[1, :] = self.spin_flip_array

        temp2 = np.zeros((2, len(self.spin_flip_array)))
        temp2[0, :] = self.x_data
        temp2[1, :] = np.ones(len(self.spin_flip_array)) * self.spin_flip_error

        self.pulsedmeasurementlogic().signal_data = temp
        self.pulsedmeasurementlogic().measurement_error = temp2

        self.pulsedmeasurementlogic().sigMeasurementDataUpdated.emit()
        return

    ############################################################################
    def save_measurement(self, save_tag):

        t1 = time.time()

        # filepath = self.savelogic().get_path_for_module('PulsedMeasurement')
        timestamp = datetime.datetime.now()

        data = OrderedDict()
        data['time'] = np.array(self.time_axis)
        if self.normalized:
            data['Norm. Counts'] = np.array(self.trace)
        else:
            data['Counts'] = np.array(self.trace)

        # create a parameters dictionary to be saved to file:
        parameters = OrderedDict()
        parameters['number bins'] = self.num_bins
        parameters['init_threshold0'] = self.init_threshold0
        parameters['init_threshold1'] = self.init_threshold1
        parameters['ana_threshold0'] = self.ana_threshold0
        parameters['ana_threshold1'] = self.ana_threshold1
        parameters['counts_per_readout'] = self.counts_per_readout
        parameters['countlength'] = self.countlength
        parameters['threshold_optimal'] = self.fit_dict['threshold_optimal']
        parameters['fidelity_total_optimal'] = self.fit_dict['fidelity_total_optimal']
        parameters['fit_dict'] = self.fit_dict

        ###### plot timetrace ###########
        plt.style.use(self.savelogic().mpl_qd_style)
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.plot(self.time_axis, self.trace, '-o', color='blue',
                 linestyle=':', linewidth=0.5, label='count trace')
        # add the spin flip probability to first plot
        # The position of the text annotation is controlled with the
        # relative offset in x direction
        rel_offset = 0.02
        if self.analyze_mode == 'mapping':
            ax1.text(1.00 + rel_offset, 0.99, 'mapping prob.: ' + str(round(self.mapped_state * 100, 2)),
                     verticalalignment='top', horizontalalignment='left', transform=ax1.transAxes, fontsize=18)
        else:
            ax1.text(1.00 + rel_offset, 0.99, 'spin flip prob.: ' + str(round(self.spin_flip_prob * 100, 2)),
                     verticalalignment='top', horizontalalignment='left', transform=ax1.transAxes, fontsize=18)
        ax1.set_xlabel('time [s]')
        ax1.set_ylabel('norm. intensity')

        print('save ssr measurement: plot1 time = {} ms'.format((time.time() - t1) / 1e-3))

        ###### plot histogram ###########
        t2 = time.time()
        ax2.hist(self.trace, self.num_bins)

        # include fit curve and fit parameters.

        # double-gaussian fit to total histogram
        fit_result = self.fit_dict['fit_result']
        hist_fit_x = self.fit_dict['fit_x']
        hist_fit_y = self.fit_dict['fit_y']

        # individual gaussians
        hist_fit_y1 = self.fit_dict['fit_y1']
        hist_fit_y2 = self.fit_dict['fit_y2']

        ax2.plot([self.threshold_optimal, self.threshold_optimal], [0, max(hist_fit_y)],
                 marker='None', linewidth=1.5, color='k')
        ax2.plot(hist_fit_x, hist_fit_y1,
                 marker='None', linewidth=1.5,  # color='o',
                 label='fit: Gaussian left')
        ax2.plot(hist_fit_x, hist_fit_y2,
                 marker='None', linewidth=1.5,  # color='o',
                 label='fit: Gaussian right')
        ax2.plot(hist_fit_x, hist_fit_y,
                 marker='None', linewidth=1.5,  # color='o',
                 label='fit: double Gaussian')

        # add then the fit result to the second plot:
        # Parameters for the text plot:
        # create the formatted fit text:
        if hasattr(fit_result, 'result_str_dict'):
            fit_res = units.create_formatted_output(fit_result.result_str_dict)
        else:
            self.savelogic().log.warning('The fit container does not contain any data '
                                  'from the fit! Apply the fit once again.')
            fit_res = ''
        # do reverse processing to get each entry in a list
        entry_list = fit_res.split('\n')
        entry_text = 'Fit results: \n'
        # for entry in entry_list:
        #     entry_text += entry + '\n'
        entry_text += 'set threshold: ' + str(round(np.mean((self.ana_threshold0, self.ana_threshold1)), 3)) + '\n'
        entry_text += 'fidelity left: ' + str(round(self.fidelity_left * 100, 2)) + '\n'
        entry_text += 'fidelity right: ' + str(round(self.fidelity_right * 100, 2)) + '\n'
        entry_text += 'fidelity total: ' + str(round(self.fidelity_total * 100, 2)) + '\n\n'
        entry_text += 'optimal threshold: ' + str(round(self.threshold_optimal, 3)) + '\n'
        entry_text += 'fidelity left (optimal): ' + str(round(self.fidelity_left_optimal * 100, 2)) + '\n'
        entry_text += 'fidelity right (optimal): ' + str(round(self.fidelity_right_optimal * 100, 2)) + '\n'
        entry_text += 'fidelity total (optimal): ' + str(round(self.fidelity_total_optimal * 100, 2))

        ax2.text(1.00 + rel_offset, 0.99, entry_text,
                 verticalalignment='top', horizontalalignment='left', transform=ax2.transAxes, fontsize=18)
        if self.normalized:
            ax2.set_xlabel('norm. intensity')
        else:
            ax2.set_xlabel('counts')
        ax2.set_ylabel('occurence')

        fig.tight_layout()

        print('save ssr measurement: plot2 time = {} ms'.format((time.time()-t2) / 1e-3))

        ###### save plots and parameter dictionary ###########

        t3 = time.time()
        self.savelogic().save_data(data, timestamp=timestamp,
                            parameters=parameters, fmt='%.15e',
                            filepath=self.savelogic().get_path_for_module('SSR'),
                            filelabel=save_tag,
                            delimiter='\t', plotfig=fig)
        print('save ssr measurement: savelogic().save_data time = {} ms'.format((time.time() - t3) / 1e-3))

        t4 = time.time()
        np.savez(self.savelogic().get_path_for_module('SSR') + '\\' + timestamp.strftime(
            '%Y%m%d-%H%M-%S') + '_' + save_tag + '_fastcounter', self.raw_data)
        # np.savez(save_directory + save_tag + '_summed_rows', tmp_signal)
        print('save ssr measurement: np.savez time = {} ms'.format((time.time() - t4) / 1e-3))

        plt.close()

        print('save ssr measurement: total time = {} ms'.format((time.time() - t1) / 1e-3))
        return



