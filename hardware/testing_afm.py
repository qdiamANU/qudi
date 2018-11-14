import numpy as np
import PyDAQmx as daq

from core.module import Base, ConfigOption

class NIAFM(): # Base):

    _modtype = 'hardware'

    analog_Z_in_channels = ['/Dev2/ai0','/Dev2/ai1','/Dev2/ai2','/Dev2/ai3']
    #analog_XY_in_channels = ['/Dev2/ai1','/Dev2/ai2']
    analog_XYZ_out_channels = ['/Dev2/ao0','/Dev2/ao1','/Dev2/a02']
    analog_XYZ_in_channels = ['/Dev2/ai4','/Dev2/ai5','/Dev2/ai6']
    clock_in_channel = '/Dev2/PFI10'
    max_clock_rate = 1000
    _RWTimeout = 10

    analog_XYZ_out_voltage_range = [[-5,5],[-5,5], [-5, 5]]
    analog_XYZ_in_voltage_range = [[-10, 10], [-10, 10], [-10, 10]]
    analog_Z_in_voltage_range = [[-5,5]]*4

    def set_up_scanner(self):

        self.analog_Z_in_task = daq.TaskHandle()
        self.analog_XYZ_out_task = daq.TaskHandle()
        self.analog_XYZ_in_task = daq.TaskHandle()

        daq.DAQmxCreateTask('', daq.byref(self.analog_Z_in_task))
        daq.DAQmxCreateTask('', daq.byref(self.analog_XYZ_out_task))
        daq.DAQmxCreateTask('', daq.byref(self.analog_XYZ_in_task))

        for n, chan in enumerate(self.analog_Z_in_channels):
            daq.DAQmxCreateAIVoltageChan(self.analog_Z_in_task,
                                         chan,
                                         '',
                                         daq.DAQmx_Val_RSE,
                                         min(self.analog_Z_in_voltage_range[n]),
                                         max(self.analog_Z_in_voltage_range[n]),
                                         daq.DAQmx_Val_Volts,
                                         '')

        for n,chan in enumerate(self.analog_XYZ_out_channels):
            daq.DAQmxCreateAOVoltageChan(self.analog_XYZ_out_task,
                                         chan,
                                         '',
                                         min(self.analog_XYZ_out_voltage_range[n]),
                                         max(self.analog_XYZ_out_voltage_range[n]),
                                         daq.DAQmx_Val_Volts,
                                         '')


        for n,chan in enumerate(self.analog_XYZ_in_channels):
            daq.DAQmxCreateAIVoltageChan(self.analog_XYZ_in_task,
                                         chan,
                                         '',
                                         daq.DAQmx_Val_RSE,
                                         min(self.analog_XYZ_in_voltage_range[n]),
                                         max(self.analog_XYZ_in_voltage_range[n]),
                                         daq.DAQmx_Val_Volts,
                                         '')


    def set_up_timing(self,NumberSamps):

        daq.DAQmxCfgSampClkTiming(self.analog_Z_in_task,
                                  self.clock_in_Z_channel,
                                  self.max_clock_rate,
                                  daq.DAQmx_Val_Rising,
                                  daq.DAQmx_Val_FiniteSamps,
                                  NumberSamps)

        daq.DAQmxCfgSampClkTiming(self.analog_XYZ_out_task,
                                  self.clock_in_XYZ_channel,
                                  self.max_clock_rate,
                                  daq.DAQmx_Val_Rising,
                                  daq.DAQmx_Val_ContSamps,
                                  NumberSamps)

        daq.DAQmxCfgSampClkTiming(self.analog_XYZ_in_task,
                                  self.clock_in_XYZ_channel,
                                  self.max_clock_rate,
                                  daq.DAQmx_Val_Rising,
                                  daq.DAQmx_Val_ContSamps,
                                  NumberSamps)

    def read_Z_in(self, numtoread):

        numread = daq.int32()

        analog_data = np.full((len(self.analog_Z_in_channels),numtoread),
                              0.,
                              dtype=np.float64)

        daq.DAQmxReadAnalogF64(self.analog_Z_in_task,
                               int(numtoread),
                               float(self._RWTimeout),
                               daq.DAQmx_Val_GroupByChannel,
                               analog_data,
                               int(len(self.analog_Z_in_channels)*numtoread),
                               daq.byref(numread),
                               None)

        return numread.value, analog_data

    def write_XY_out(self,voltages,length=1,start=False):

        numout = daq.int32()

        daq.DAQmxWriteAnalogF64(self.analog_XY_out_task,
                                int(length),
                                start,
                                float(self._RWTimeout),
                                daq.DAQmx_Val_GroupByChannel,
                                voltages,
                                daq.byref(numout),
                                None)
        return numout.value

    def _stop_analog_output(self):

        if self.analog_XY_out_task is None:
            return -1
        retval =0
        try:
            daq.DAQmxStopTask(self.analog_XY_out_task)
        except:
            self.log.exception('Error stopping analog input')
            retval = -1

        return retval




if __name__ == '__main__':
    afm = NIAFM()
    afm.set_up_scanner()
    afm.set_up_timing(100)
    import matplotlib.pyplot as plt
    dummydata = np.array([np.linspace(0,1,100)]*3)
    plt.plot(dummydata[0,:])
    nout = afm.write_XY_out(dummydata,100)
    n, data = afm.read_Z_in(101)
    daq.DAQmxStartTask(afm.analog_XY_out_task)
    print(n,'\n',data,'\n',data.shape)
    plt.plot(data[0,:])
    plt.show()