import time
import numpy as np
import pickle

import TimeTagger as tt
import pulsestreamer as ps


def calculate_rabi_power(current_power, current_rabi, desired_rabi):
    '''
    simple numerical calculation to determine a desired power for a given desired rabi frequency for a given power and rabi frequency
    @param float current_power: current power in dBm
    @param float current_rabi: current Rabi frequency
    @param float desired_rabi: desired Rabi frequency

    @return float: desired_power in dBm
    '''
    power_ratio = (current_rabi / desired_rabi) ** (-2)
    dBdifference = np.log10(power_ratio) * 10.
    return current_power + dBdifference


def set_scan_range_for_coscan(scan_size=1.e-6):
    '''
    sets a confocal scan range with the current confocal position in the lower left corner.

    @param float scan_size: desired scan size in units of metres
    '''
    current_pos = confocal._scanning_logic.get_position()
    confocal._mw.x_min_InputWidget.setValue(current_pos[0])
    confocal._mw.y_min_InputWidget.setValue(current_pos[1])
    xmax = confocal._mw.x_min_InputWidget.value() + scan_size
    ymax = confocal._mw.y_min_InputWidget.value() + scan_size
    confocal._mw.x_max_InputWidget.setValue(xmax)
    confocal._mw.y_max_InputWidget.setValue(ymax)
    confocal.change_x_image_range()
    confocal.change_y_image_range()
    return


def refocus(poi=None, wait_time=1):
    '''
    refocuses using the poimanagerlogic for a given poi and returns the location, waits for focus to finish.
    @param str poi: poimanager's string label, if None - refocuses at current location

    @return array: XYZ array of poi's location
    '''
    if poi is None:
        confocal.refocus_clicked()
    else:
        poimanagerlogic.optimise_poi_position(poi)

    time.sleep(1)
    while True:
        if not optimizerlogic.module_state() == 'locked':
            time_last = time.time()
            break
        else:
            time.sleep(wait_time)
    if poi is None:
        return confocal._scanning_logic.get_position()
    else:
        return poimanagerlogic.poi_positions[poi]


def measure_ODMR_freqs(poi=None, do_refocus=True, odmr_time=7, odmr_power=-5, Rmag=None, save_data=False):
    # work out magnetic field to guess sweep variables
    if Rmag is None:
        try:
            Rmag = np.sqrt(np.sum((nimagnet._current_fields) ** 2))
        except:
            Rmag = 0
    f1, f2 = 2870 - Rmag * 2.8, 2870. + Rmag * 2.8
    Nsteps = 100
    fstep = (f2 - f1) / Nsteps

    # ODMR settings
    odmrlogic.set_sweep_parameters(f1 * 1e6, f2 * 1e6, fstep * 1e6, odmr_power)
    odmrlogic.set_runtime(odmr_time)

    # recfocus if desired
    if do_refocus:
        refocus(poi)
    elif poi is not None:
        poimanagerlogic.go_to_poi(poi)
    else:
        print('no poi to go to')
        return -1

    # preapre loop breaker and start
    odmrlogic.loop_breaker = False
    odmrlogic.start_odmr_scan()
    while True:
        if odmrlogic.loop_breaker == True or odmrlogic.module_state() == 'idle':
            break
        else:
            odmrlogic.do_fit('Two Lorentzian dips')
            time.sleep(1)
    odmrlogic.do_fit('Two Lorentzian dips')

    if save_data:
        odmrlogic.save_odmr_data(tag=str(poi) + "_odmr_")

    # return fit frequencies
    return odmrlogic.fc.current_fit_param['l0_center'].value, odmrlogic.fc.current_fit_param['l1_center'].value


def measure_rabi_frequency(poi, freq=2.87e9, power=5., runtime=10, do_refocus=True, save_data=False):
    # do refocus if necessary
    if poi is None:
        refocus(poi)
    else:
        if do_refocus:
            poimanagerlogic.optimise_poi_position(poi)
        else:
            poimanagerlogic.go_to_poi(poi)

    pulsedmeasurement._pa.ana_param_invoke_settings_CheckBox.setChecked(True)
    pulsedmeasurementlogic.set_microwave_settings({'power': float(power), 'frequency': freq, 'use_ext_microwave': True})
    pulsedmeasurement._pm.samplo_buttons['rabi'].click()
    time.sleep(0.1)
    pulsedmeasurementlogic.toggle_pulsed_measurement(True)
    while True:
        if pulsedmeasurementlogic.elapsed_time > runtime or pulsedmeasurementlogic.module_state() == 'idle':
            break
        else:
            time.sleep(2)
            fit = pulsedmeasurementlogic.do_fit('sin')
    pulsedmeasurementlogic.toggle_pulsed_measurement(False)
    fit = pulsedmeasurementlogic.do_fit('sin')
    rabi_freq = fit[-1].values['frequency']

    if save_data:
        pulsedmeasurementlogic.save_measurement_data(tag=str(poi) + "_rabi_")

    return rabi_freq


def measure_spin_echo(poi, freq=2.87e9, power=5., rabi_period=140.e9, runtime=10, do_refocus=True, save_data=False):
    # do refocus if necessary
    if poi is None:
        refocus(poi)
    else:
        if do_refocus:
            poimanagerlogic.optimise_poi_position(poi)
        else:
            poimanagerlogic.go_to_poi(poi)

    pulsedmeasurement._pa.ana_param_invoke_settings_CheckBox.setChecked(True)
    pulsedmeasurementlogic.set_microwave_settings({'power': float(power), 'frequency': freq, 'use_ext_microwave': True})
    pulsedmasterlogic.set_generation_parameters({'rabi_period': float(rabi_period)})
    pulsedmeasurement._pm.samplo_buttons['hahnecho'].click()
    pulsedmeasurementlogic.set_alternative_data_type('Delta')
    pulsedmeasurementlogic.toggle_pulsed_measurement(True)

    while True:
        if pulsedmeasurementlogic.elapsed_time > runtime or pulsedmeasurementlogic.module_state() == 'idle':
            break
        else:
            time.sleep(2)
            fit = pulsedmeasurementlogic.do_fit('exp')
            fit2 = pulsedmeasurementlogic.do_fit('exp', use_alternative_data=True)
    pulsedmeasurementlogic.toggle_pulsed_measurement(False)
    fit = pulsedmeasurementlogic.do_fit('exp')
    fit2 = pulsedmeasurementlogic.do_fit('exp', use_alternative_data=True)
    t2 = fit2[-1].values['lifetime']

    if save_data:
        pulsedmeasurementlogic.save_measurement_data(tag=str(poi) + "_hahnecho_")

    return t2


def align_field(poi=None, sweep_direction='azimuthal', iterations=3, grid_points=5, Rmag=50,
                polar=np.arccos(1 / np.sqrt(3)), azimuthal=0, do_refocus=True, refocus_interval=100, odmr_power=0,
                odmr_time=5):
    # work out microwave limits
    f1, f2 = 2870 - Rmag * 2.8, 2870. + Rmag * 2.8
    Nsteps = 100
    fstep = (f2 - f1) / Nsteps
    odmrlogic.set_sweep_parameters(f1 * 1e6, f2 * 1e6, fstep * 1e6, odmr_power)
    odmrlogic.set_runtime(odmr_time)
    refocus(poi)
    refocus_last = time.time()

    # azimuthal sweep first
    start = 0
    stop = np.pi * 2
    breaknext = False
    odmrlogic.loop_breaker = False
    best_angles, allangles, allfrequencies, allsplittings, allmeans = [], [], [], [], []

    for run in range(iterations):
        data, fits = [], []
        sweep_variable = np.linspace(start, stop, grid_points, endpoint=False)

        for angle in sweep_variable:

            if sweep_direction == 'azimuthal':
                nimagnet.write_polar_fields(Rmag, polar, angle)
            elif sweep_direction == 'polar':
                nimagnet.write_polar_fields(Rmag, angle, azimuthal)
            else:
                print('sweep direction not understood!')
                return -1
            time.sleep(0.5)

            if time.time() - refocus_last > refocus_interval:
                refocus(poi)
                refocus_last = time.time()
            time.sleep(0.5)

            odmrlogic.start_odmr_scan()

            while True:
                if odmrlogic.module_state() == 'idle' or odmrlogic.loop_breaker == True:
                    break
                else:
                    odmrlogic.do_fit('Two Lorentzian dips')
                    time.sleep(1)
            odmrlogic.do_fit('Two Lorentzian dips')

            data.append([odmrlogic.odmr_plot_x, odmrlogic.odmr_plot_y])
            fits.append(odmrlogic.fc.current_fit_param)

            if odmrlogic.loop_breaker == True:
                breaknext = True
                break

        splitting, frequencies, means = [], [], []
        for f in fits:
            try:
                x1 = f['l0_center'].value
                x2 = f['l1_center'].value
                frequencies.append([x1, x2])
                splitting.append(np.abs(x1 - x2))
                means.append((x1 + x2) / 2.)
            except KeyError:
                frequencies.append([np.nan, np.nan])
                splitting.append(np.nan)
                means.append(np.nan)

        max_splitting_indice = np.argmax(splitting)
        best_angle = sweep_variable[max_splitting_indice]
        angle_spacing = sweep_variable[1] - sweep_variable[0]
        start = best_angle - 1.5 * angle_spacing
        stop = best_angle + 1.5 * angle_spacing

        best_angles.append(best_angle)
        allfrequencies.append(frequencies)
        allsplittings.append(splitting)
        allmeans.append(means)
        allangles.append(sweep_variable)

        if breaknext == True:
            break

    return best_angles, allangles, allfrequencies, allsplittings, allmeans


def coscan(resolution, scanner_clock_frequency, scan_size, do_refocus=True, refocus_interval=180, nap_mode=False,
           do_measurement=False, nap_height=1.e-6):
    scannerlogic.loop_breaker = False
    scannerlogic.update_do_afm(True)

    if do_measurement:
        contactmeasurement = ContactMeasurement(scanner_clock_frequency, resolution)
        noncontactmeasurement = ContactMeasurement(scanner_clock_frequency, resolution)

    contact_measurement_data = []
    non_contact_measurement_data = []

    if nap_mode:
        controltip = ControlTip()

    if do_refocus:
        confocal.refocus_clicked()
        time.sleep(0.1)
        while True:
            if not optimizerlogic.module_state() == 'locked':
                break
            else:
                time.sleep(0.1)
        time.sleep(1)

    x0, y0, z0 = confocal._scanning_logic.get_position()
    last_focus = time.time()
    coords, data = [], []
    afm_coords = []
    focus_history = []
    stepsize = scan_size / float(resolution - 1)

    # set up counter and scanner
    res1 = myafmnicard.set_up_scanner_clock(scanner_clock_frequency)
    if res1 < 0:
        print('Scanner clock error!')
    res2 = myafmnicard.set_up_scanner()
    if res2 < 0:
        print('Scanner definition error!')

    if res1 < 0 or res2 < 0:
        print("Scanner set up fail\n abort coscan")
        return

    # scan y rows

    if nap_mode:
        N_rows = int(resolution * 2)
    else:
        N_rows = int(resolution)

    for j in range(N_rows):

        if scannerlogic.loop_breaker == True:
            break

        # if doing this nap gargabe - need to raise the tip
        if nap_mode:
            if not j % 2 == 0:
                controltip.raise_tip()

        # if even line, move x in positive direction
        xi, yi, zi = confocal._scanning_logic.get_position()
        if j % 2 == 0:
            xline = np.linspace(xi, xi + scan_size, resolution)
        else:
            xline = np.linspace(xi, xi - scan_size, resolution)

        # scan the line in x
        yline = np.linspace(yi, yi, resolution)
        zline = np.linspace(zi, zi, resolution)
        line = np.vstack([xline, yline, zline])

        # need to define the afm line
        if myafmnicard._do_afm_scan:
            xline_afm = xline.copy() - x0  # only want shift from start of scan for AFM
            yline_afm = yline.copy() - y0
            zline_afm = np.zeros(xline_afm.shape)
            if nap_mode:
                if not j % 2 == 0:
                    zline_afm += 0  # height_data + nap_height

            afm_line = np.vstack([xline_afm, yline_afm, zline_afm])
            afm_coords.append(afm_line[:, ::1].T)
        else:
            afm_line = None
            afm_coords = None

        # set up pulsed measurement parameters here
        if do_measurement:
            if nap_mode:
                if j % 2 == 0:  # set up contact measurement
                    contactmeasurement.run()
                else:  # set up weak nap measurement
                    noncontactmeasurement.run()  # some sort of Ramsey measurement
            else:
                contactmeasurement.run()

        # scan the line
        linedata = myafmnicard.scan_line(line, afm_line, pixel_clock=True)

        # get pulsed measurement data
        if do_measurement:
            if nap_mode:
                if j % 2 == 0:
                    contact_measurement_data.append(contactmeasurement.stop())
                else:
                    non_contact_measurement_data.append(noncontactmeasurement.stop())
            else:
                contact_measurement_data.append(contactmeasurement.stop())
        else:
            contact_measurement_data.append(None)
            non_contact_measurement_data.append(None)

        if nap_mode:
            pass
            # height_data = linedata[:,1] # 1=sensor, 2=piezo from previous scan
        else:
            height_data = None

        data.append(linedata)
        coords.append(line.T)

        # step y - get current coords
        xc = xline[-1]
        yc = yline[-1]
        zc = zline[-1]

        if nap_mode:
            if not j % 2 == 0:
                # set Zmod offset to 0
                myafmnicard.afm_set_position(z=0)

                # put tip back down
                controltip.lower_tip()

        if nap_mode:
            if j % 2 == 0:  # at the end of an even row - stay still - odd row goes back over the top
                pass
            else:
                confocal._scanning_logic.set_position('scanner', x=xc, y=yc + stepsize, z=zc)
        else:  # without nap mode - contine as normal step y and the end of even and odd rows
            confocal._scanning_logic.set_position('scanner', x=xc, y=yc + stepsize, z=zc)

        # step y position of AFM too
        if myafmnicard._do_afm_scan:
            xc = xline_afm[-1]
            yc = yline_afm[-1]
            if nap_mode:
                if j % 2 == 0:  # at the end of an even row - stay still - odd row goes back over the top
                    pass
                else:
                    myafmnicard.afm_set_position(x=xc, y=yc + stepsize)
            else:  # without nap mode - contine as normal step y and the end of even and odd rows
                myafmnicard.afm_set_position(x=xc, y=yc + stepsize)

        if do_refocus and time.time() - last_focus > refocus_interval and not j % 2 == 0:  # for nap mode - only do refocus at end of odd line - so any shift in focus doesn't effect nap retrace

            # stop clocks and scanner to pass control to optimizer
            myafmnicard.close_scanner()
            myafmnicard.close_scanner_clock()
            time.sleep(0.5)

            # perform refocus
            confocal.refocus_clicked()
            # wait until complete
            time.sleep(0.5)
            while True:
                if not optimizerlogic.module_state() == 'locked':
                    break
                else:
                    time.sleep(0.5)
            last_focus = time.time()
            focus_history.append(confocal._scanning_logic.get_position())
            # give it time before restarting scanner and clocks
            time.sleep(1.)

            # re-set up counter and scanner
            res1 = myafmnicard.set_up_scanner_clock(scanner_clock_frequency)
            if res1 < 0:
                print('Scanner clock error!')
            res2 = myafmnicard.set_up_scanner()
            if res2 < 0:
                print('Scanner definition error!')
            if res1 < 0 or res2 < 0:
                print("Scanner set up fail\n abort coscan - refocus")
                return

    # clean up

    if len(focus_history) > 1:
        focus_history = np.array(focus_history)
        offset_from_start = focus_history[-1, :] - focus_history[0, :]
    else:
        offset_from_start = [0, 0, 0]
    myafmnicard.close_scanner()
    myafmnicard.close_scanner_clock()
    confocal._scanning_logic.set_position('scanner', x=x0 + offset_from_start[0], y=y0 + offset_from_start[1],
                                          z=z0 + offset_from_start[2])
    myafmnicard.afm_set_position(x=0, y=0, )
    # confocal.refocus_clicked()
    confocal._scanning_logic.update_do_afm(False)

    scan_results = np.array(data), np.array(coords), np.array(
        afm_coords), contact_measurement_data, non_contact_measurement_data, focus_history
    data_dict = {}
    keys = ['confocal_data', 'confocal_coords', 'afm_coords', 'contact_data', 'non_contact_data', 'focus_history']
    for i, d in enumerate(scan_results):
        data_dict.update({keys[i]: d})

    return data_dict


def save_coscan_data(data_dict, tag=None):
    path = savelogic.get_path_for_module('coscan')
    timestamp = datetime.datetime.now()
    extension = '.pkl'
    filename = '\\'.join([path, timestamp.strftime('%Y%m%d-%H%M-%S')])
    if tag is not None:
        filename += '_' + str(tag)
    filename += extension
    print(filename)
    with open(filename, 'wb') as f:
        pickle.dump(data_dict, f)


class ControlTip():
    '''
    using an Aout connected to InFastOffset on the Aslyum Crosspoint panel - a DC offset of the setpoint to a non-manageable value (eg set point of -1.0 V) will force the AFM Z PID loop to retract the tip
    seems to work OK for both AC and DC modes
    '''

    def raise_tip(self):
        aout.write_volts([-2.5])

    def lower_tip(self):
        aout.write_volts([0.0])


class Scanning_CW_ODMR():
    '''
    Performs a CW ODMR scan for a single line
    '''

    def setup(self, fstart, fstop, fstep, power, scanner_frequency=10, resolution=20):
        # set up microwaves
        smiqset = smiqmicrowave.set_sweep(fstart, fstop, fstep, power)

        # define channels on pulsestreamer
        ps_smiq_trigger_channel = 7
        ps_laser_channel = mypulser._laser_channel
        ps_sequence_channel = 0

        # define channels on timetagger
        apd_channel = 0
        tt_smiq_trigger = 5
        tt_sequence = 3
        tt_scanner_clock = 4

        # construct pulse sequences
        safety_factor = 0.99
        num_freqs = len(np.arange(smiqset[0], smiqset[1], smiqset[2]))
        num_freqs += 2  # mysterious extra two points for SMIQ
        time_per_pixel = 1. / scanner_frequency * safety_factor
        time_per_freq = time_per_pixel / num_freqs
        pulse_width = 100.e-9
        sync_channel_pulses = [(int(pulse_width * 1.e9), 1), (int((time_per_pixel - pulse_width) * 1.e9), 0)]
        smiq_trigger_pulses = [(int(time_per_freq / 2 * 1e9), 1), (int(time_per_freq / 2 * 1e9), 0)] * num_freqs
        laser_pulses = [(int(time_per_pixel * 1e9), 1)]

        # create sequencer
        self.seq = mypulser.pulse_streamer.createSequence()
        self.seq.setDigital(ps_smiq_trigger_channel, smiq_trigger_pulses)
        self.seq.setDigital(ps_laser_channel, laser_pulses)
        self.seq.setDigital(ps_sequence_channel, sync_channel_pulses)

        # set up time-tagger
        num_scanner_cycles = (resolution + 1)
        self.pixel_counter = tt.Countrate(myfastcounter._tagger, [tt_smiq_trigger, tt_scanner_clock])
        self.odmr_counter = tt.TimeDifferences(myfastcounter._tagger,
                                               apd_channel,
                                               start_channel=tt_smiq_trigger,
                                               next_channel=tt_smiq_trigger,
                                               sync_channel=tt.CHANNEL_UNUSED,
                                               binwidth=int(time_per_freq * 1e12),
                                               n_bins=1,
                                               n_histograms=num_freqs * (num_scanner_cycles * 2))

    def start(self):
        # reset and start smiq
        smiqmicrowave.reset_sweeppos()
        smiqmicrowave.sweep_on()

        # reset and start pulser and timetagger
        self.pixel_counter.clear()
        self.odmr_counter.clear()
        self.odmr_counter.start()
        self.pixel_counter.start()
        mypulser.pulse_streamer.setTrigger(
            ps.TriggerStart.HARDWARE_RISING)  # user hardware trigger (scanner clock edge)
        mypulser.pulse_streamer.stream(self.seq, n_runs=1,
                                       final=mypulser._laser_mw_on_state)  # only one frequency sweep per pixel

    def stop(self):
        self.odmr_counter.stop()
        self.pixel_counter.stop()
        smiqmicrowave.off()
        mypulser.pulse_streamer.setTrigger(ps.TriggerStart.IMMEDIATE)
        data = self.odmr_counter.getData()
        pixel_counts = self.pixel_counter.getCountsTotal()

        return data, pixel_counts


# due to the large array of possible measurements, I have put them into this class so you can change them here if you want
class ContactMeasurement():

    def __init__(self, scanner_frequency, resolution):
        self.scanning_cw_odmr = Scanning_CW_ODMR()
        self.scanner_frequency = scanner_frequency
        self.resolution = resolution

    def run(self):
        self.scanning_cw_odmr.setup(2.7e9, 3.05e9, 3.e9, -5, self.scanner_frequency,
                                    self.resolution)  # self,fstart,fstop,fstep,power,scanner_frequency=10,resolution=20
        time.sleep(0.5)
        self.scanning_cw_odmr.start()

    def stop(self):
        return self.scanning_cw_odmr.stop()


class NonContactMeasurement():

    def __init__(self, *args):
        pass

    def run(self):
        pass

    def stop(self):
        return None

