# -*- coding: utf-8 -*-
"""
A hardware module for communicating with the fast counter FPGA.

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

from interface.fast_counter_interface import FastCounterInterface
import numpy as np
import TimeTagger as tt
from core.module import Base, ConfigOption
import os


class TimeTaggerFastCounter(Base, FastCounterInterface):
    """ Hardware class to controls a Time Tagger from Swabian Instruments.

    Example config for copy-paste:

    fastcounter_timetagger:
        module.Class: 'swabian_instruments.timetagger_fast_counter.TimeTaggerFastCounter'
        timetagger_channel_apd_0: 0
        timetagger_channel_apd_1: 1
        timetagger_channel_detect: 2
        timetagger_channel_sequence: 3
        timetagger_sum_channels: 4

    """
    _modclass = 'TimeTaggerFastCounter'
    _modtype = 'hardware'

    _channel_apd_0 = ConfigOption('timetagger_channel_apd_0', missing='error')
    _channel_apd_1 = ConfigOption('timetagger_channel_apd_1')
    _channel_detect = ConfigOption('timetagger_channel_detect', missing='error')
    _channel_sequence = ConfigOption('timetagger_channel_sequence')
    _sum_channels = ConfigOption('timetagger_sum_channels', False)

    def on_activate(self):
        """ Connect and configure the access to the FPGA.
        """
        self._tagger = tt.createTimeTagger()
        self._tagger.reset()

        self._number_of_gates = int(100)
        self._bin_width = 1
        self._record_length = int(4000)

        # self._tagger.setTestSignal(0, True)

        if self._sum_channels == True:
            self._channel_combined = tt.Combiner(self._tagger, channels=[self._channel_apd_0, self._channel_apd_1])
            self._channel_apd = self._channel_combined.getChannel()
        else:
            self._channel_apd = self._channel_apd_0

        self.log.info('TimeTagger (fast counter) configured to use  channel {0}'
                      .format(self._channel_apd))

        #self._tagger.setTestSignal(0, False)

        self.statusvar = 0

    def get_constraints(self):
        """ Retrieve the hardware constrains from the Fast counting device.

        @return dict: dict with keys being the constraint names as string and
                      items are the definition for the constaints.

         The keys of the returned dictionary are the str name for the constraints
        (which are set in this method).

                    NO OTHER KEYS SHOULD BE INVENTED!

        If you are not sure about the meaning, look in other hardware files to
        get an impression. If still additional constraints are needed, then they
        have to be added to all files containing this interface.

        The items of the keys are again dictionaries which have the generic
        dictionary form:
            {'min': <value>,
             'max': <value>,
             'step': <value>,
             'unit': '<value>'}

        Only the key 'hardware_binwidth_list' differs, since they
        contain the list of possible binwidths.

        If the constraints cannot be set in the fast counting hardware then
        write just zero to each key of the generic dicts.
        Note that there is a difference between float input (0.0) and
        integer input (0), because some logic modules might rely on that
        distinction.

        ALL THE PRESENT KEYS OF THE CONSTRAINTS DICT MUST BE ASSIGNED!
        """

        constraints = dict()

        # the unit of those entries are seconds per bin. In order to get the
        # current binwidth in seonds use the get_binwidth method.
        constraints['hardware_binwidth_list'] = [1 / 1000e6]

        # TODO: think maybe about a software_binwidth_list, which will
        #      postprocess the obtained counts. These bins must be integer
        #      multiples of the current hardware_binwidth

        return constraints

    def on_deactivate(self):
        """ Deactivate the FPGA.
        """
        if self.module_state() == 'locked':
            self.pulsed.stop()
        self.pulsed.clear()
        self.pulsed = None

    def configure(self, bin_width_s, record_length_s, number_of_gates=0):

        """ Configuration of the fast counter.

        @param float bin_width_s: Length of a single time bin in the time trace
                                  histogram in seconds.
        @param float record_length_s: Total length of the timetrace/each single
                                      gate in seconds.
        @param int number_of_gates: optional, number of gates in the pulse
                                    sequence. Ignore for not gated counter.

        @return tuple(binwidth_s, gate_length_s, number_of_gates):
                    binwidth_s: float the actual set binwidth in seconds
                    gate_length_s: the actual set gate length in seconds
                    number_of_gates: the number of gated, which are accepted
        """


        self._number_of_gates = number_of_gates
        self._bin_width = bin_width_s * 1e9
        self._record_length = 1 + int(record_length_s / bin_width_s)
        self.statusvar = 1
        # Detection channel input voltage is quite low, due to signal being split between laser and timetagger
        # Need to lower the trigger level on this channel, to e.g. 0.9 V
        self._tagger.setTriggerLevel(self._channel_detect, 0.9)
        #print('TriggerLevel ch{} = {} V'.format(self._channel_detect, self._tagger.getTriggerLevel(self._channel_detect)))
        #number_of_gates = 50

        print('tt configure')
        print('tagger={}, click_channel={}, start/next_channel={}'.format(
            self._tagger, self._channel_apd, self._channel_detect))
        print('bin_width={} ns, '.format(self._bin_width))
        binwidth_ps = int(np.round(self._bin_width * 1000))
        print('binwidth={} ps, n_bins={}, n_histograms={}'.format(binwidth_ps, int(self._record_length), number_of_gates))
        self.pulsed = tt.TimeDifferences(
            tagger=self._tagger,
            click_channel=self._channel_apd,
            start_channel=self._channel_detect,
            next_channel=self._channel_detect,
            sync_channel=tt.CHANNEL_UNUSED,
            binwidth=int(np.round(self._bin_width * 1000)),
            n_bins=int(self._record_length),
            n_histograms=number_of_gates)

        self.pulsed.stop()

        return (bin_width_s, record_length_s, number_of_gates)

    def start_measure(self):
        """ Start the fast counter. """
        print('tt start measure')
        self.module_state.lock()
        self.pulsed.clear()
        self.pulsed.start()
        self.statusvar = 2
        return 0

    def stop_measure(self):
        """ Stop the fast counter. """
        print('tt stop measure')
        if self.module_state() == 'locked':
            self.pulsed.stop()
            self.module_state.unlock()
        self.statusvar = 1
        return 0

    def pause_measure(self):
        """ Pauses the current measurement.

        Fast counter must be initially in the run state to make it pause.
        """
        if self.module_state() == 'locked':
            self.pulsed.stop()
            self.statusvar = 3
        return 0

    def continue_measure(self):
        """ Continues the current measurement.

        If fast counter is in pause state, then fast counter will be continued.
        """
        if self.module_state() == 'locked':
            self.pulsed.start()
            self.statusvar = 2
        return 0

    def is_gated(self):
        """ Check the gated counting possibility.

        Boolean return value indicates if the fast counter is a gated counter
        (TRUE) or not (FALSE).
        """
        return True

    def get_data_trace(self):
        """ Polls the current timetrace data from the fast counter.

        @return numpy.array: 2 dimensional array of dtype = int64. This counter
                             is gated the the return array has the following
                             shape:
                                returnarray[gate_index, timebin_index]

        The binning, specified by calling configure() in forehand, must be taken
        care of in this hardware class. A possible overflow of the histogram
        bins must be caught here and taken care of.
        """
        print('tt get_data_trace')
        return np.array(self.pulsed.getData(), dtype='int64')


    def get_status(self):
        """ Receives the current status of the Fast Counter and outputs it as
            return value.

        0 = unconfigured
        1 = idle
        2 = running
        3 = paused
        -1 = error state
        """
        return self.statusvar

    def get_binwidth(self):
        """ Returns the width of a single timebin in the timetrace in seconds. """
        width_in_seconds = self._bin_width * 1e-9
        return width_in_seconds

################################### Methods for SSR interface ####################################

    def configure_ssr_counter(self, counts_per_readout=None, countlength=None):
        """ Configuration of the fast counter.

        necessary for some fast counters. Does not appear necessary for TimeTagger

                @param int counts_per_readout: optional, number of readouts for one measurement
                @param int cycles: countlength, number of measurements


                @return tuple(binwidth_s, record_length_s, preset, cycles, sequences ):
                            binwidth_s: float the actual set binwidth in seconds
                            gate_length_s: the actual record length in seconds
                            preset: number of readout per measurement
                            cycles: number of measurements
                            sequences: number of sequences
        """

        pass

    def change_sweep_mode(self, gated=True):
        """
                necessary for some fast counters. Does not appear necessary for TimeTagger
        """
        pass
