import PyDAQmx as daq
import numpy as np

from core.module import Base
from core.configoption import ConfigOption


# output = np.array([3.])
# num_written = daq.int32()
# xhandle = daq.TaskHandle()
# daq.DAQmxCreateTask('xDaqTask', daq.byref(xhandle))
# daq.DAQmxCreateAOVoltageChan(xhandle, '/Dev3/AO1', '', -10., 10., daq.DAQmx_Val_Volts, '')
# #daq.DAQmxStartTask(xhandle)
# daq.DAQmxWriteAnalogF64(xhandle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, output, daq.byref(num_written),None)


class NI_AOut(Base):
    '''
    Aout from an National Instruments card (static only).
    '''

    _modtype = 'NI_AOut'
    _modclass = 'hardware'
    _channels = ConfigOption('channels', missing='warn')

    def on_activate(self):

        self._num_channels = len(self._channels)
        self._channels = ','.join(self._channels)
        print(self._channels)

        self.handle = daq.TaskHandle()
        daq.DAQmxCreateTask('AOutDaq', daq.byref(self.handle))

        daq.DAQmxCreateAOVoltageChan(self.handle, self._channels, '', -10, 10, daq.DAQmx_Val_Volts, '')
        volts = np.zeros(self._num_channels)
        self.num_written = daq.int32()
        daq.DAQmxWriteAnalogF64(self.handle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, volts,
                                daq.byref(self.num_written), None)
        self._current_volts = volts

    def on_deactivate(self):
        self.set_all_zero()
        daq.DAQmxStopTask(self.handle)
        del self.handle

    def write_volts(self, volts):
        if not len(volts) == self._num_channels:
            self.log.error('nimagnet - Number of channels not equal to number of voltages {} {}'.format(len(volts),
                                                                                                        self._num_channels))
            return -1
        else:
            volts = np.array(volts).astype(np.float64)
            daq.DAQmxWriteAnalogF64(self.handle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, volts,
                                    daq.byref(self.num_written), None)
            self._current_volts = volts
            return self.num_written

    def set_all_zero(self, ):
        volts = np.zeros(self._num_channels)
        return self.write_volts(volts)