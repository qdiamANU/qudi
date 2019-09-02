import PyDAQmx as daq
import numpy as np

from core.module import Base, ConfigOption
from interface.magnet_interface import MagnetInterface

# output = np.array([3.])
# num_written = daq.int32()
# xhandle = daq.TaskHandle()
# daq.DAQmxCreateTask('xDaqTask', daq.byref(xhandle))
# daq.DAQmxCreateAOVoltageChan(xhandle, '/Dev3/AO1', '', -10., 10., daq.DAQmx_Val_Volts, '')
# #daq.DAQmxStartTask(xhandle)
# daq.DAQmxWriteAnalogF64(xhandle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, output, daq.byref(num_written),None)


class NIMagnet(Base):

    '''
    A 3 axis electromagnet controlled by analog outs from an National Instruments (static only).
    '''

    _modtype = 'NIMagnet'
    _modclass = 'hardware'
    _channels = ConfigOption('channels',missing='error')
    _xy_calib = 1/(150./10.)  #150 G is at 0.3 A and that corresponds to 10V
    _z_calib = 1./(130./10.)   #zfield is a bit weaker at the position of the diamond


    def on_activate(self):

        self._num_channels = len(self._channels)
        self._channels = ','.join(self._channels)
        print(self._channels)


        self.handle = daq.TaskHandle()
        daq.DAQmxCreateTask('MagnetDaq', daq.byref(self.handle))

        daq.DAQmxCreateAOVoltageChan(self.handle,self._channels,'',-10,10, daq.DAQmx_Val_Volts, '')
        volts = np.zeros(self._num_channels)
        self.num_written = daq.int32()
        daq.DAQmxWriteAnalogF64(self.handle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, volts, daq.byref(self.num_written),None)
        self._current_volts = volts
        self._current_fields = self._volts_to_fields(volts)

    def on_deactivate(self):
        self.set_all_zero()
        daq.DAQmxStopTask(self.handle)
        del self.handle

    def write_volts(self, volts):
        if not len(volts) == self._num_channels:
            self.log.error('nimagnet - Number of channels not equal to number of voltages {} {}'.format(len(volts),self._num_channels))
            return -1
        else:
            volts = np.array(volts).astype(np.float64)
            daq.DAQmxWriteAnalogF64(self.handle, 1, True, 5., daq.DAQmx_Val_GroupByChannel, volts, daq.byref(self.num_written), None)
            self._current_volts = volts
            self._current_fields = self._volts_to_fields(volts)
            return self.num_written

    def set_all_zero(self,):
        volts = np.zeros(self._num_channels)
        return self.write_volts(volts)

    def write_fields(self, fields):
        fields = np.array(fields)
        return self.write_volts(self._fields_to_volts(fields))

    def _fields_to_volts(self,fields):
        return fields*np.array([self._xy_calib,self._xy_calib,self._z_calib])

    def _volts_to_fields(self,volts):
        return volts/np.array([self._xy_calib,self._xy_calib,self._z_calib])

    def write_polar_fields(self,r,theta,phi):
        Bx = r*np.sin(theta)*np.cos(phi)
        By = r*np.sin(theta)*np.sin(phi)
        Bz = r*np.cos(theta)
        return self.write_fields(np.array([Bx,By,Bz]))





