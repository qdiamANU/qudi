# -*- coding: utf-8 -*-
"""
This is just a dummy hardware class to be used with TemplateLogic.

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

from core.module import Base
from core.configoption import ConfigOption
import visa, time

class StanfordDS345(Base):

    _gpib_address = ConfigOption('gpib_address', missing='error')
    _gpib_timeout = ConfigOption('gpib_timeout', 10, missing='warn')
    _gpib_baud_rate = ConfigOption('gpib_baud_rate', None)

    def on_activate(self):

        #connect
        """ Initialisation performed during activation of the module. """
        self._gpib_timeout = self._gpib_timeout * 1000
        # trying to load the visa connection to the module
        self.rm = visa.ResourceManager()
        try:
            if self._gpib_baud_rate is None:
                self._gpib_connection = self.rm.open_resource(self._gpib_address,
                                                            timeout=self._gpib_timeout)
            else:
                self._gpib_connection = self.rm.open_resource(self._gpib_address,
                                                            timeout=self._gpib_timeout,
                                                            baud_rate=self._gpib_baud_rate)
        except:
            self.log.error('This is Stanford: could not connect to GPIB address >>{}<<.'
                           ''.format(self._gpib_address))
            raise

    def on_deactivate(self):

        # turn off
        self.output_off()
        return 0

    def _command_wait(self, command_str):
        """
        Writes the command in command_str via GPIB and waits until the device has finished
        processing it.

        @param command_str: The command to be written
        """
        self._gpib_connection.write(command_str)
        self._gpib_connection.write('*WAI')
        while int(float(self._gpib_connection.query('*OPC?'))) != 1:
            time.sleep(0.2)
        return

    def set_wave_type(self, wavetype='sin'):

        # sqaure wave etc
        if wavetype == 'sin':
            self._gpib_connection.write('FUNC 0')
        elif wavetype == 'square':
            self._gpib_connection.write('FUNC 1')

    def set_amplitude(self,amplitude=1,type='Vpp'):

        if type == 'Vpp':
            self._gpib_connection.write('AMPL {0}VP'.format(amplitude))
        else:
            print('Stanford: voltage units not understood')

    def set_offset(self,offset=0.):

        pass

    def set_frequency(self,frequency=1.e3):

        self._gpib_connection.write('FREQ {0}'.format(frequency))


    def get_amplitude(self):
        return self._gpib_connection.query('AMPL?')

    def get_frequency(self):
        return self._gpib_connection.query('FREQ?')