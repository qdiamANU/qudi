import numpy as np
import re

import PyDAQmx as daq

from core.module import Base, Connector, ConfigOption
#from interface.slow_counter_interface import SlowCounterInterface
#from interface.slow_counter_interface import SlowCounterConstraints
#from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface
#from interface.confocal_scanner_interface import ConfocalScannerInterface


class PulseTrainCounter(Base, ODMRCounterInterface):
    """Outputs pulsed train and performs gated count. Used for CW ODMR"""

    _modclass = 'ODMRCounterDummy'
    _modtype = 'hardware'

    # connectors
    fitlogic = Connector(interface='FitLogic')

    # config options
    _clock_frequency = ConfigOption('clock_frequency', 100, missing='warn')
    _counter_in = ConfigOption('counter_in',  missing='error')
    _counter_out = ConfigOption('counter_out', missing='error')
    _photon_source = ConfigOption('photon_source', missing='error')

    MaxCount = 1.e7
    DutyCycle = 0.8

    def on_activate(self):

        pass

    def set_up_odmr(self, counter_channel=None, photon_source=None,
                    clock_channel=None, odmr_trigger_channel=None):

        self._counter_channel = counter_channel
        self._photon_source = photon_source
        self._clock_channel = clock_channel
        self._odmr_trigger_channel = odmr_trigger_channel

        self._counter_task = daq.TaskHandle()
        daq.DAQmxCreateTask('', daq.byref(self._counter_task))

            # ctr0 is used to count photons. Used to count ticks in N+1 gates
        daq.DAQmxCreateCIPulseWidthChan(self._counter_task,
                                        self._counter_channel,
                                        '',
                                        0.,
                                        MaxCount * DutyCycle / self._clock_frequency,
                                        daq.DAQmx_Val_Ticks,
                                        daq.DAQmx_Val_Rising, '')

        daq.DAQmxSetCIPulseWidthTerm(self._counter_task, self._CounterIn, self._CounterOut + 'InternalOutput')
        daq.DAQmxSetCICtrTimebaseSrc(self._counter_task, self._CounterIn, self._TickSource)


    def set_up_odmr_clock(self, clock_frequency=None, clock_channel=None):

        self._clock_frequency = clock_frequency
        self._clock_channel = clock_channel

        # nidaq Tasks
        self._clock_Task = daq.TaskHandle()
        daq.DAQmxCreateTask('', daq.byref(self._clock_task))


        # ctr1 generates a continuous square wave with given duty cycle. This serves simultaneously
        # as sampling clock for AO (update DAC at falling edge), and as gate for counter (count between
        # rising and falling edge)
        daq.DAQmxCreateCOPulseChanFreq(self._clock_task,
                                       self._clock_channel, '',
                                       daq.DAQmx_Val_Hz,
                                       daq.DAQmx_Val_Low,
                                       0,
                                       self._clock_frequency,
                                       DutyCycle)

    def set_odmr_length(self, length=100):

        daq.DAQmxCfgImplicitTiming(self._clock_task, daq.DAQmx_Val_ContSamps, SampleLength)
        daq.DAQmxCfgImplicitTiming(self._counter_task, daq.DAQmx_Val_FiniteSamps, SampleLength)

        # read samples from beginning of acquisition, do not overwrite
        daq.DAQmxSetReadRelativeTo(self._counter_task, daq.DAQmx_Val_CurrReadPos)
        daq.DAQmxSetReadOffset(self._counter_task, 0)
        daq.DAQmxSetReadOverWrite(self._counter_task, daq.DAQmx_Val_DoNotOverwriteUnreadSamps)


    def run(self):

        self._CIData = np.empty((SampleLength,), dtype=np.uint32)
        self._CINread = daq.int32()

        self._SampleLength = SampleLength
        self._TaskTimeout = 4 * SampleLength / f
        self._RWTimeout = RWTimeout

        daq.DAQmxStartTask(self._counter_task)
        daq.DAQmxStartTask(self._clock_task)
        daq.DAQmxWaitUntilTaskDone(self._counter_task, daq.c_double(self._TaskTimeout))
        daq.DAQmxReadCounterU32(self._counter_task,
                                    self._SampleLength,
                                    self._RWTimeout,
                                    self._CIData,
                                    self._SampleLength,
                                    daq.byref(self._CINread), None)
        daq.DAQmxStopTask(self._clock_task)
        daq.DAQmxStopTask(self._counter_task)
        return self._CIData

    def clear(self):
        daq.DAQmxClearTask((self._counter_task))
        daq.DAQmxClearTask((self._clock_task))
        del self._counter_task
        del self._clock_task

    def __del__(self):
        try:
            self.clear()
        except Exception as e:
            print(str(e))


def test():
    counter = PulseTrainCounter('/Dev1/ctr0', '/Dev1/ctr1', '/Dev1/PFI8')
    counter.configure(100, 0.001)
    print(counter.run())
    counter.clear()
    #del (counter)

if __name__ == '__main__':
    test()