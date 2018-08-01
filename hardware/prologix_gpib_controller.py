import serial
import re

class PrologixController():
    
    def __init__(self, usb_com_address='COM4', visa_address=28, timeout=0.5):
             
        self.s = serial.Serial(usb_com_address, timeout = timeout)
        self.s.write('++mode 1\n')
        self.s.write('++addr '+str(visa_address) + '\n')
        self.s.write('++eoi 1\n')
        self.s.write('++auto 0\n')
        
        
    def change_instrument(self, visa_address=28):
        self.s.write('++mode 1\n')
        self.s.write('++addr '+str(visa_address) + '\n')
        self.s.write('++eoi 1\n')
        self.s.write('++auto 0\n')
    
    def query(self, cmd):
        
        out = None
        #self.s.write('++auto 1\n')
        if str(cmd):
            self.s.write(str(cmd)+'\n')
        else:
            print 'Command not string-able'
        try:
            self.s.write('++read eoi\n')
            out = self.s.read(256)
            out = re.search(r'(.*?)\n', out).group(1).rstrip()       
        except:
            pass
        #self.s.write('++auto 0\n')
        
        return out        
        
            
    def write(self, cmd):
        if str(cmd):
            try:
                self.s.write(str(cmd)+'\n')
            except:
                print 'Prologix write failed\tYou suck!'
        else:
            print 'Command not string-able'
            
    def __del__(self):
        del self.s
            
    
if __name__ == '__main__':
    
    instr = MyPrologixController()
    