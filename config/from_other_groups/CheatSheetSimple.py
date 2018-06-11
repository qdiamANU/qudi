#######################################################################
###############                                         ###############
###############                INIT                     ###############
###############                                         ###############
#######################################################################


cd 'c:/software/diamond'
import diamond
import numpy, threading, pickle, time, pylab
import TimeTagger
import Pulsed
%matplotlib qt


#######################################################################
###############                                         ###############
###############      Pulsed ODMR                         ###############
###############                                         ###############
#######################################################################

odseq = [(['mw', 'mw_gate'], 120), ([], 100), (['laser', 'aom'], 300), ([], 500)]
Pulsed.PulseGen.Sequence(odseq)


#######################################################################
###############                                         ###############
###############     fit ODMR                           ###############
###############                                         ###############
#######################################################################

import Analysis
import numpy
import pylab
%matplotlib
odmr = diamond.widget.ODMRWidget# .Measurement
fig = pylab.figure()
ax = fig.add_subplot(111)
lower_bound = 60
# lower_bound = 150
# upper_bound = 200
upper_bound = 90
# upper_bound = odmr.x.shape[0]
fitpeak = Analysis.ODMR(odmr.x[:-2], odmr.y[:-2])#[lower_bound:upper_bound],odmr.y[lower_bound:upper_bound])
a=fitpeak.FitLorentzian(odmr.y[:-2], odmr.x[:-2])#[lower_bound:upper_bound],odmr.x[lower_bound:upper_bound])
# fitpeak = Analysis.ODMR(odmr.x[lower_bound:upper_bound],odmr.y[lower_bound:upper_bound])
# a=fitpeak.FitLorentzian(odmr.y[lower_bound:upper_bound],odmr.x[lower_bound:upper_bound])
fres = a[0][0]
ax.plot(odmr.x,odmr.y, label = 'data')
# ax.plot(odmr.x,fitpeak.Lorentzian(*a[0])(odmr.x), label = 'fit, fres = %s'%str(a[0][0]))
ax.plot(odmr.x,fitpeak.Lorentzian(*a[0])(odmr.x), label = 'fit, fres = %s'%str(a[0][0]))
# ax.plot(odmr.x[lower_bound:upper_bound],fitpeak.Lorentzian(*a[0])(odmr.x[lower_bound:upper_bound]), label = 'fit, fres = %s'%str(a[0][0]))
ax.legend()
pylab.show()

print 'contrast: ' + str(numpy.abs((100/a[0][3])*a[0][2]/(numpy.pi*a[0][1])))
print 'B = ' + str((2.87e9+fres)/2.8/1e6)



#######################################################################
###############                                         ###############
###############                 RABI                    ###############
###############                                         ###############
#######################################################################

import numpy
import pickle
import Pulsed
power = -15.0
power_deer = -40
binwidth=10
AcqTime=1500
LaserOn=1000
LaserWait=1000

t0 = 2
t1 = 500
dt =  5

rabi_seq = "[([ch('mw'), ch('mwy'), ch('mw_gate')], tau)]"

try:
    rabi.Stop()
except:
    pass
rabi = Pulsed.Generic_double_resonance(fres, power, 0, power_deer, binwidth, 100, 0, 100, 0, AcqTime, LaserOn, LaserWait, numpy.arange(float(t0), float(t1), float(dt)), rabi_seq, rabi_seq, False)
a = rabi.Start()

# import numpy
# import pickle
# import Pulsed
# # power = -1.0
# power = 5
# #fres = 2.7078e9
# power_deer = 5
# binwidth=10
# AcqTime=2000
# LaserOn=2000
# LaserWait=1000
#
# t0 = 10
# t1 = 400
# dt =  4
# rabi_seq = "[([ch('mw'), ch('mw_gate')], tau)]"
#
# try:
#     rabi.Stop()
# except:
#     pass
# rabi = Pulsed.Generic_double_resonance(fres, power, fres, power, binwidth, 100, 0, 100, 0, AcqTime, LaserOn, LaserWait, numpy.arange(float(t0), float(t1), float(dt)), rabi_seq, '', False)
# rabi.Start()


#############
rabi.Stop()
############

#PLOT

rd = rabi.export_rabi()
rd.AnalyzeSequence(300, 250, 50, 680) #diode laser
# rd.AnalyzeSequence(200, 250, 50, 680) #dpss laser

# fit and plot the Rabi fringes
rd.FitwithDecay()
# rd.Fit()
rd.Plot()

diamond.RabiPeriod = rd.RabiPeriod
diamond.Rabix0 = rd.Rabix0
diamond.power = power
diamond.fres = fres
diamond.AcqTime = AcqTime
diamond.binwidth = binwidth
diamond.LaserOn = LaserOn
diamond.LaserWait = LaserWait

diamond.widget.ODMRWidget.RabiPeriod.setValue(rd.RabiPeriod)
diamond.widget.ODMRWidget.RabiOffset.setValue(rd.Rabix0)
diamond.widget.ODMRWidget.PulsedPower.setValue(power)

#SAVE
prefix = "D:/data/170510_tips/NV11/scratching/"
try:
    os.makedirs(prefix)
except:
    pass

fname = "rabi_after_1st_scratch_scan"

file = prefix + fname + ".pyd"
fil = open(file, 'wb')
pickle.dump(rd, fil,1 )
fil.close()
rd.Plot(file=prefix + fname + ".png")


#######################################################################
###############                                         ###############
###############                 HAHN                    ###############
###############                                         ###############
#######################################################################

#Wait = numpy.linspace(100, 200, 10)
Wait = numpy.asarray(numpy.linspace(5, 5000, 10), dtype=int)
# hahn = "[([ch('mw', 'mw_gate'], Tpi2), ([   ], tau), (['mw', 'mw_gate'], Tpi), ([   ], tau), (['mw', 'mw_gate'], Tpi2)]"
# hahni = "[(['mw', 'mw_gate'], Tpi2), ([   ], tau), (['mw', 'mw_gate'], Tpi), ([   ], tau), (['mw', 'mw_gate'], T3pi2)]"
hahn = "[([ch('mw'), ch('mw_gate')], Tpi2), ([   ], tau), ([ch('mw'), ch('mw_gate')], Tpi), ([   ], tau), ([ch('mw'), ch('mw_gate')], Tpi2)]"
hahni = "[([ch('mw'), ch('mw_gate')], Tpi2), ([   ], tau), ([ch('mw'), ch('mw_gate')], Tpi), ([   ], tau), ([ch('mw'), ch('mw_gate')], T3pi2)]"
# mes = Pulsed.Generic_double_resonance(fres, power, fres, power, binwidth, rd.RabiPeriod,rd.Rabix0, rd.RabiPeriod, rd.Rabix0, AcqTime, LaserOn, LaserWait, Wait, hahn, hahni)
mes = Pulsed.Generic_double_resonance(fres, power, 0, -40, binwidth, rd.RabiPeriod,rd.Rabix0, 50, 0, AcqTime, LaserOn, LaserWait, Wait, hahn, hahni)
mes.Start()


#############
mes.Stop()

md = mes.export()
md.AnalyzeSequence(300, 250, 50, 680)

prefix = "D:/data/170510_tips/NV11/scratching/"
try:
    os.mkdir(prefix)
except:
    pass

fname = "hahn_echo_after_1st_scratch"
file = prefix + fname + ".pyd"
fil = open(file, 'wb')
pickle.dump(md, fil,1 )
fil.close()

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [mus]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(md.t[0::2]/1.0e3, md.z[0::2]-md.z[1::2], label = '0-pi')
# ax.plot(d.t[1::2], , label = 'pi')
ax.legend()
fig.savefig(prefix + fname + "_contrast.png")
pylab.close(fig)
# pylab.show()

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [ns]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(md.t[0::2], md.z[0::2], label = 'pi/2')
ax.plot(md.t[1::2], md.z[1::2], label = '3pi/2')
ax.legend()
# pylab.show()
fig.savefig(prefix + fname + "_bgplot.png")
pylab.close(fig)



#######################################################################
###############                                         ###############
###############                 T1                      ###############
###############                                         ###############
#######################################################################

# LaserWait = 1000
T1_seq = "[([  ], tau), ([ ], Tpi), ([  ], 100)]"  # pulse at the end
T1_iseq = "[([  ], tau), ([ch('mw'), ch('mw_gate')], Tpi), ([  ], 100)]" # pulse at the end

Wait = numpy.linspace(20, 400000, 5)
# Wait = numpy.array([1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000])
T1 = Pulsed.Generic_double_resonance(fres, power, 0, -40, binwidth, rd.RabiPeriod, rd.Rabix0, rd.RabiPeriod, rd.Rabix0, AcqTime, LaserOn, LaserWait, Wait, T1_seq, T1_iseq)
T1.Start()


#############
T1.Stop()

#save the measurement
md = T1.export()
md.AnalyzeSequence(300, 250, 50, 680)


prefix = "D:/data/170510_tips/NV11/scratching/"
try:
    os.makedirs(prefix)
except:
    pass


fname = "T1_after_1st_scratch"
file = prefix + fname+ ".pyd"
fil = open(file, 'wb')
pickle.dump(md, fil,1 )
fil.close()
md.Plot(file=prefix + fname + ".png")

import pickle
fil = open(file, 'rb')
d = pickle.load(fil)
fil.close()

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [ns]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(d.t[0::2], d.z[0::2], label = '0')
ax.plot(d.t[1::2], d.z[1::2], label = 'pi')
ax.legend()
fig.savefig(prefix + fname + "_bgplot.png")
pylab.close(fig)
# pylab.show()

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [mus]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(md.t[0::2]/1.0e3, md.z[0::2]-md.z[1::2], label = '0-pi')
# ax.plot(d.t[1::2], , label = 'pi')
ax.legend()
fig.savefig(prefix + fname + "_contrast.png")
pylab.close(fig)
# pylab.show()

# fig = pylab.figure()
# ax = fig.add_subplot(111)
# #ax.set_xscale('log')
# ax.set_xlabel('Time [ns]')
# ax.set_ylabel('Fluorescence [1]')
# ax.plot(d.t[0::2], d.z[0::2], label = '0')
# ax.plot(d.t[1::2], d.z[1::2], label = 'pi')
# ax.legend()
# fig.savefig(prefix + fname + "_bgplot-log.png")
# pylab.close(fig)
# pylab.show()


#######################################################################
###############                                         ###############
###############                 XY8                     ###############
###############                                         ###############
#######################################################################

def XY8_seq( N , end32=False):
    '''generates a XY8_seq sequence of order N'''
    #abbrevation definitions
    Pi2i = "(['mw', 'mw_gate'], Tpi2), "
    Pi2f = "(['mw', 'mw_gate'], Tpi2), "
    Pi32f = "(['mw', 'mw_gate'], T3pi2), "
    Pi = "(['mw', 'mw_gate'], Tpi), "
    Piy = "(['mwy', 'mw_gate'], Tpi), "
    delay = "([  ], tau), "
    delay2 = "([  ], tau/2.), "
    end = "([  ], 300) "
    if end32:
        return "[" + Pi2i + (N) *((delay2 + Pi + delay + Piy + delay + Pi + delay + Piy + delay2) + (delay2 + Piy + delay + Pi + delay + Piy + delay + Pi + delay2)) + Pi32f + "]"

    else:
        return "[" + Pi2i + (N) *((delay2 + Pi + delay + Piy + delay + Pi + delay + Piy + delay2) + (delay2 + Piy + delay + Pi + delay + Piy + delay + Pi + delay2)) + Pi2f + "]"





Wait = numpy.linspace(10,4000, 40)
# Wait = numpy.linspace(10, 1000, 30)
# Wait = numpy.concatenate((numpy.linspace(1500, 2500, 20),numpy.linspace(7500, 8500, 20)))
XY8 = XY8_seq(2, False)
XY8i = XY8_seq(2, True)
mes = Pulsed.Generic_double_resonance(fres, power, 0, -40, binwidth, rd.RabiPeriod,rd.Rabix0, rd.RabiPeriod, rd.Rabix0, AcqTime, LaserOn, LaserWait, Wait, XY8, XY8i, False, False)
mes.Start()
#############


mes.Stop()
Pulsed.Night()
#save the measurement
md = mes.export()
md.AnalyzeSequence(300, 250, 50, 680)

prefix = "C:/data/151014_n15_mtssl_bsa/C1/nv3/realigned/"
try:
    os.mkdir(prefix)
except:
    pass



fname = "xy8_order_1"
file = prefix + fname + ".pyd"
fil = open(file, 'wb')
pickle.dump(md, fil,1 )
fil.close()
md.Plot(file=prefix + fname + ".png")

#convert hahn/fid echo data files into graphics
# splits into pi/2 and 3pi/2 pulses
# import autopilot
import pickle
fil = open(file, 'rb')
d = pickle.load(fil)
fil.close()

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [ns]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(d.t[0::2], d.z[0::2], label = 'pi/2')
ax.plot(d.t[1::2], d.z[1::2], label = '3pi/2')
ax.legend()
fig.savefig(prefix + fname + "_bgplot.png")
pylab.close(fig)
# pylab.show()


fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time [ns]')
ax.set_ylabel('Fluorescence [1]')
ax.plot(d.t[0::2], (d.z[1::2]-d.z[0::2]), label = 'pi/2-3pi/2')
# ax.plot(d.t[1::2], , label = 'pi')
ax.legend()
fig.savefig(prefix + fname + "_contrast.png")
# pylab.close(fig)
pylab.show()



#######################################################################
###############                                         ###############
###############         Home AFM z to surface           ###############
###############                                         ###############
#######################################################################

#1: engage AFM to surface
#2:
afm_handler.start_piezo_loops([2])
#3: "Read ALL" in AFM PIDS Loop Panel
#4:
afm_handler.home_piezo_loop(2)



#######################################################################
###############                                         ###############
###############         scanning AFM                    ###############
###############                                         ###############
#######################################################################

import afm_scanning
import Pulsed
import numpy as np
afm_scanner = afm_scanning.AFMScanner(None, Pulsed)

scan_height = 0.05e-6
center = [x00+150e-9, y00-450e-9]
size = [0.5e-6, 0.5e-6]
pts = [40, 40]
integration_time = 1

x0 = center[0]-size[0]/2.
y0 = center[1]-size[1]/2.

plane = afm_scanning.Plane()
plane.origin = np.array([x0,y0,scan_height])
plane.vectors[0] = np.array([size[0], 0, 0])
plane.vectors[1] = np.array([0, size[1], 0])
plane.pixels = pts

afm_scanner.StartScan(plane, integrationtime=integration_time, measurement_type='pulsed fluorescence')





