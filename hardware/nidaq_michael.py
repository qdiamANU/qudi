import numpy as np
import re

import PyDAQmx as daq

from core.module import Base, ConfigOption
from interface.slow_counter_interface import SlowCounterInterface
from interface.slow_counter_interface import SlowCounterConstraints
from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface
from interface.confocal_scanner_interface import ConfocalScannerInterface


class PulseTrainCounter:
    """Outputs pulsed train and performs gated count. Used for CW ODMR"""

    def __init__(self, CounterIn, CounterOut, TickSource):

        self._CounterIn = CounterIn
        self._CounterOut = CounterOut
        self._TickSource = TickSource

    def configure(self, SampleLength, SecondsPerPoint, DutyCycle=0.9, MaxCounts=1e7, RWTimeout=1.0):

        if hasattr(self, '_CITask') or hasattr(self, '_COTask'):
            self.clear()

        f = 1. / SecondsPerPoint

        # nidaq Tasks
        self._COTask = daq.TaskHandle()
        self._CITask = daq.TaskHandle()

        daq.DAQmxCreateTask('', daq.byref(self._COTask))
        daq.DAQmxCreateTask('', daq.byref(self._CITask))

        # ctr1 generates a continuous square wave with given duty cycle. This serves simultaneously
        # as sampling clock for AO (update DAC at falling edge), and as gate for counter (count between
        # rising and falling edge)
        daq.DAQmxCreateCOPulseChanFreq(self._COTask,
                                       self._CounterOut, '',
                                       daq.DAQmx_Val_Hz,
                                       daq.DAQmx_Val_Low,
                                       0,
                                       f,
                                       DutyCycle)

            # ctr0 is used to count photons. Used to count ticks in N+1 gates
        daq.DAQmxCreateCIPulseWidthChan(self._CITask,
                                        self._CounterIn,
                                        '',
                                        0.,
                                        MaxCounts * DutyCycle / f,
                                        daq.DAQmx_Val_Ticks,
                                        daq.DAQmx_Val_Rising, '')

        daq.DAQmxSetCIPulseWidthTerm(self._CITask, self._CounterIn, self._CounterOut + 'InternalOutput')
        daq.DAQmxSetCICtrTimebaseSrc(self._CITask, self._CounterIn, self._TickSource)

        daq.DAQmxCfgImplicitTiming(self._COTask, daq.DAQmx_Val_ContSamps, SampleLength)
        daq.DAQmxCfgImplicitTiming(self._CITask, daq.DAQmx_Val_FiniteSamps, SampleLength)

        # read samples from beginning of acquisition, do not overwrite
        daq.DAQmxSetReadRelativeTo(self._CITask, daq.DAQmx_Val_CurrReadPos)
        daq.DAQmxSetReadOffset(self._CITask, 0)
        daq.DAQmxSetReadOverWrite(self._CITask, daq.DAQmx_Val_DoNotOverwriteUnreadSamps)

        self._CIData = np.empty((SampleLength,), dtype=np.uint32)
        self._CINread = daq.int32()

        self._SampleLength = SampleLength
        self._TaskTimeout = 4 * SampleLength / f
        self._RWTimeout = RWTimeout

    def run(self):
        daq.DAQmxStartTask(self._CITask)
        daq.DAQmxStartTask(self._COTask)
        daq.DAQmxWaitUntilTaskDone(self._CITask, daq.c_double(self._TaskTimeout))
        daq.DAQmxReadCounterU32(self._CITask,
                                    self._SampleLength,
                                    self._RWTimeout,
                                    self._CIData,
                                    self._SampleLength,
                                    daq.byref(self._CINread), None)
        daq.DAQmxStopTask(self._COTask)
        daq.DAQmxStopTask(self._CITask)
        return self._CIData

    def clear(self):
        daq.DAQmxClearTask((self._CITask))
        daq.DAQmxClearTask((self._COTask))
        del self._CITask
        del self._COTask

    def __del__(self):
        try:
            self.clear()
        except Exception as e:
            print(str(e))


def test():
    counter = PulseTrainCounter('/Dev1/ctr0', '/Dev1/ctr1', '/Dev1/PFI8')
    counter.configure(100, 0.001)
    print(counter.run())
    print(counter._CITask)
    print(counter._COTask)
    counter.clear()
    #del (counter)

if __name__ == '__main__':
    test()