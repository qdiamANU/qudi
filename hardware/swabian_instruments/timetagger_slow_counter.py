# -*- coding: utf-8 -*-

"""
This file contains the Qudi hardware module to use TimeTagger as a counter.

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

import TimeTagger as tt
import time
import numpy as np


from core.module import Base, ConfigOption
from interface.slow_counter_interface import SlowCounterInterface
from interface.slow_counter_interface import SlowCounterConstraints
from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface

class TimeTaggerSlowCounter(Base, SlowCounterInterface, ODMRCounterInterface):

    """ Using the TimeTagger as a counter."""

    _modtype = 'TTCounter'
    _modclass = 'hardware'

    _channel_apd_0 = ConfigOption('channel_apd_0', missing='error')
    _channel_apd_1 = ConfigOption('channel_apd_1', None)
    _channel_apd = _channel_apd_0
    _sum_channels = ConfigOption('timetagger_sum_channels', False)
    _odmr_trigger_channel = ConfigOption('channel_odmr_trigger', None, missing='warn')

    _nvalues = 4096
    def on_activate(self):
        """ Start up TimeTagger interface
        """
        self._tagger = tt.createTimeTagger()
        self._tagger.reset()
        self._count_frequency = 50  # Hz

        self._tagger.setTestSignal([0, 1], False)

        self.odmr_counter = None

        if self._sum_channels and self._channel_apd_1 is None:
            self.log.error('Cannot sum channels when only one apd channel given')

        # self._mode can take 3 values:
        # 0: single channel, no summing
        # 1: single channel, summed over apd_0 and apd_1
        # 2: dual channel for apd_0 and apd_1
        if self._sum_channels:
            self._mode = 1
        elif self._channel_apd_1 is None:
            self._mode = 0
        else:
            self._mode = 2


    def on_deactivate(self):
        """ Shut down the TimeTagger.
        """
        self.close_counter()
        self.close_odmr()
        return 0

    def set_up_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the TimeTagger for timing

        @param float clock_frequency: if defined, this sets the frequency of
                                      the clock
        @param string clock_channel: if defined, this is the physical channel
                                     of the clock

        @return int: error code (0:OK, -1:error)
        """

        self._count_frequency = clock_frequency
        return 0

    def set_up_counter(self,
                       counter_channels=None,
                       sources=None,
                       clock_channel=None,
                       counter_buffer=None):
        """ Configures the actual counter with a given clock.

        @param str counter_channels: optional, physical channel of the counter
        @param str sources: optional, physical channel where the photons
                                  are to count from
        @param str clock_channel: optional, specifies the clock channel for the
                                  counter
        @param int counter_buffer: optional, a buffer of specified integer
                                   length, where in each bin the count numbers
                                   are saved.

        @return int: error code (0:OK, -1:error)
        """

        # currently, parameters passed to this function are ignored -- the channels used and clock frequency are
        # set at startup
        #self._tagger.setTestSignal(0, True)
        #self._tagger.setTestSignal(1, True)
        if self._mode == 1:
            channel_combined = tt.Combiner(self._tagger, channels = [self._channel_apd_0, self._channel_apd_1])
            self._channel_apd = channel_combined.getChannel()

            self.counter = tt.Counter(
                self._tagger,
                channels=[self._channel_apd],
                binwidth=int((1 / self._count_frequency) * 1e12),
                n_values=1
            )
        elif self._mode == 2:
            self.counter0 = tt.Counter(
                self._tagger,
                channels=[self._channel_apd_0],
                binwidth=int((1 / self._count_frequency) * 1e12),
                n_values=1
            )

            self.counter1 = tt.Counter(
                self._tagger,
                channels=[self._channel_apd_1],
                binwidth=int((1 / self._count_frequency) * 1e12),
                n_values=1
            )
        else:
            self._channel_apd = self._channel_apd_0
            self.counter = tt.Counter(
                self._tagger,
                channels=[self._channel_apd],
                binwidth=int((1 / self._count_frequency) * 1e12),
                n_values=1
            )

        self.log.info('set up counter with {0}'.format(self._count_frequency))
        return 0

    def get_counter_channels(self):
        """ Returns the list of counter channel names.

                @return tuple(str): channel names

                Most methods calling this might just care about the number of channels, though.
                """
        if self._mode < 2:
            return [self._channel_apd, ]
        else:
            return [self._channel_apd_0, self._channel_apd_1]

    def get_constraints(self):
        """ Get hardware limits the device

        @return SlowCounterConstraints: constraints class for slow counter

        FIXME: ask hardware for limits when module is loaded
        """
        constraints = SlowCounterConstraints()
        constraints.max_detectors = 2
        constraints.min_count_frequency = 1e-3
        constraints.max_count_frequency = 10e9
        constraints.counting_mode = [CountingMode.CONTINUOUS]
        return constraints

    def get_counter(self, samples=None):
        """ Returns the current counts per second of the counter.

        @param int samples: if defined, number of samples to read in one go

        @return numpy.array(uint32): the photon counts per second
        """

        time.sleep(2 / self._count_frequency)
        if self._mode < 2:
            return self.counter.getData() * self._count_frequency
        else:
            return np.array([self.counter0.getData() * self._count_frequency,
                             self.counter1.getData() * self._count_frequency])

    def close_counter(self):
        """ Closes the counter and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        return 0

    def close_clock(self):
        """ Closes the clock and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        return 0


    #ODMR counter methods
    def set_up_odmr_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the timetagger to give the timing.

        @param float clock_frequency: if defined, this sets the frequency of the
                                      clock
        @param str clock_channel: if defined, this is the physical channel of
                                  the clock

        @return int: error code (0:OK, -1:error)
        """
        return 0

    def set_up_odmr(self, counter_channel=None, photon_source=None,
                    clock_channel=None, odmr_trigger_channel=None):
        """ Configures the actual counter with a given clock.

        @param str counter_channel: if defined, this is the physical channel of
                                    the counter
        @param str photon_source: if defined, this is the physical channel where
                                  the photons are to count from
        @param str clock_channel: if defined, this specifies the clock for the
                                  counter
        @param str odmr_trigger_channel: if defined, this specifies the trigger
                                         output for the microwave

        @return int: error code (0:OK, -1:error)
        """
        # currently, parameters passed to this function are ignored -- the channels used and clock frequency are
        # set at startup
        self._tagger.setTestSignal(0, False)
        self._tagger.setTestSignal(1, False)
        if self._mode == 0:
            self._channel_apd = self._channel_apd_0
            self.odmr_counter = tt.CountBetweenMarkers(
                self._tagger,
                click_channel=self._channel_apd,
                begin_channel=self._odmr_trigger_channel,
                end_channel=self._odmr_trigger_channel+8,
                n_values=self._nvalues
            )
        elif self._mode == 1:
            channel_combined = tt.Combiner(self._tagger, channels=[self._channel_apd_0, self._channel_apd_1])
            self._channel_apd = channel_combined.getChannel()
            self.odmr_counter = tt.CountBetweenMarkers(
                self._tagger,
                click_channel=channel_combined,
                begin_channel=self._channel_apd,
                end_channel=self._odmr_trigger_channel+8,
                n_values=self._nvalues)
        else:
            self.log.error('Cannot do the dual channel mode for ODMR')
            return -1

        self.log.info('set up odmr counter with')
        return 0

    def set_odmr_length(self, length=100):
        """Set up the trigger sequence for the ODMR and the triggered microwave.

        @param int length: length of microwave sweep in pixel

        @return int: error code (0:OK, -1:error)
        """
        if length > self._nvalues:
            self.log.error('ODMR length is higher than the hardcoded limit of {} pixels for the timetagger!'.format(
                self._nvalues))
            return -1
        else:
            return 0

    def count_odmr(self, length=100):
        """ Sweeps the microwave and returns the counts on that sweep.

        @param int length: length of microwave sweep in pixel

        @return float[]: the photon counts per second
        """

        self.set_odmr_length(length)

        #start_count = time.time()
        for i in range(50):
            bin_widths = self.odmr_counter.getBinWidths()
            if np.count_nonzero(bin_widths) < length + 1:
                time.sleep(0.001)
                continue
            else:
                break
        #end_count = time.time()

        data = self.odmr_counter.getData()[:length]
        self.odmr_counter.clear()
        #print('collect data: {:0.3f} s'.format(end_count - start_count))

        return np.reshape(data, (1, data.size))

    def close_odmr(self):
        """ Close the odmr and clean up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        if not self.odmr_counter == None:
            self.odmr_counter.stop()
        self._tagger.clearOverflows()
        return 0

    def close_odmr_clock(self):
        """ Close the odmr and clean up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        return 0

    def get_odmr_channels(self):
        """ Return a list of channel names.

        @return list(str): channels recorded during ODMR measurement
        """
        return ['', ]
