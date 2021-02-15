
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from scipy import signal
import numpy as np
import statistics
from scipy.signal import argrelextrema
from pathlib import Path
import lumicks.pylake as lk

"""define the functions of the subprocess"""


# open a folder containing raw data and lead through the analysis process
def start_subprocess(analysis_folder, timestamp, Files, input_settings, input_format, export_data, input_fitting, output_q):
    global filename_i
    global save_folder
    global start_time
    global figure1
    global Force_Distance
    global Frequency_value
    global filteredDistance_ready
    global filteredForce
    global results_F, PD_start_F, F_mm2, F_mm2_STD2_positive, F_mm2_STD2_negative
    global results_PD, PD_start_PD, PD_mm2, PD_mm2_STD2_positive, PD_mm2_STD2_negative

    save_folder = analysis_folder
    start_time = timestamp
    # iterate through the files in the selected folder
    i = 0
    total_results = []

    # proceed differently with h5 and csv files
    while i < len(Files):
        if input_format['CSV'] == 1:
            df = pd.read_csv(Files[i])
            directory_i = Path(Files[i])
            filename_i = directory_i.name[:-4]
            # access the raw data
            Force_1x = df.to_numpy()[:, 0]
            Distance_1x = df.to_numpy()[:, 1]
            # accessing the data frequency from user input
            Frequency_value = input_settings['data_frequency']
            Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)
        else:
            with h5py.File(Files[i], "r") as f:
                directory_i = Path(Files[i])
                filename_i = directory_i.name[:-3]

                # access the raw data
                if input_format['HF'] == 1:
                    Force_1x = f.get("Force HF/Force 1x")
                    Distance_1x = f.get("Distance/Piezo Distance")
                    # accessing the data frequency from the h5 file
                    Frequency_value = Force_1x.attrs['Sample rate (Hz)']
                    Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

                elif input_format['LF'] == 1:
                    load_force = f.get("Force LF/Force 1x")
                    Force_1x = load_force[:]['Value'][:]
                    load_distance = f.get("Distance/Distance 1")[:]
                    Distance_1x = load_distance['Value'][:]
                    Force_Distance = np.column_stack((Force_1x, Distance_1x))

                    # calculating the data frequency based on start- and end-time of the measurement
                    size_F_LF = len(Force_1x)
                    stop_time_F_LF = load_force.attrs['Stop time (ns)']
                    start_time_F_LF = load_force.attrs['Start time (ns)']

                    Frequency_value = size_F_LF / ((stop_time_F_LF - start_time_F_LF) / 10**9)
                    Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

        # Export down sampled and smoothened FD values
        if export_data['export_SMOOTH'] == 1:
            filename = save_folder + "/" + filename_i + "_smooth_" + start_time + ".csv"
            np.savetxt(filename, Force_Distance_um, delimiter=",")
        else:
            pass

        trim_data(Force_Distance, input_settings['F_min'])

        create_derivation(input_settings)

        """find steps based on force derivation"""
        results_F = []
        PD_start_F = []
        F_mm2 = []
        F_mm2_STD2_positive = []
        F_mm2_STD2_negative = []

        find_steps_F(
            input_settings,
            filename_i,
            Force_Distance,
            derivation_array,
            results_F,
            PD_start_F,
            F_mm2,
            F_mm2_STD2_positive,
            F_mm2_STD2_negative
            )

        """find steps based on distance derivation"""
        results_PD = []
        PD_start_PD = []
        PD_mm2 = []
        PD_mm2_STD2_positive = []
        PD_mm2_STD2_negative = []

        find_steps_PD(
            input_settings,
            filename_i,
            Force_Distance,
            derivation_array,
            results_PD,
            PD_start_PD,
            PD_mm2,
            PD_mm2_STD2_positive,
            PD_mm2_STD2_negative
            )

        results_F_list = list(results_F)

        filename_steps_F = save_folder + "/" + filename_i + "_steps_F_" + start_time + ".csv"
        if export_data['export_STEPS_F'] == 1:
            steps_results_F = pd.DataFrame(results_F_list)
            steps_results_F.to_csv(filename_steps_F, index=False, header=True)
        else:
            pass
        # except:
        #     results_F=[]
        #     err_F = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Force steps')
        #     display_output(err_F)
        #     pass

        results_PD_list = list(results_PD)

        filename_steps_PD = save_folder + "/" + filename_i + "_steps_PD_" + start_time + ".csv"
        if export_data['export_STEPS_PD'] == 1:
            steps_results_PD = pd.DataFrame(results_PD_list)
            steps_results_PD.to_csv(filename_steps_PD, index=False, header=True)

        # except:
        #     err_PD = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Distance steps')
        #     display_output(err_PD)
        #     results_PD=[]
        #     pass
        save_figure(export_data['export_PLOT'])

        try:
            common_steps = find_common_steps(results_F_list, results_PD_list)
            # print(common_steps)
            filename_common_steps = save_folder + "/" + filename_i + "_common_steps_" + start_time + ".csv"
            common_steps_results = pd.DataFrame(common_steps)
            common_steps_results.to_csv(filename_common_steps, index=False, header=True)
        except:
            err_FCS = str("Error in finding common steps" + str(filename_i) + '\n' 'There was an error in finding common steps')
            output_q.put(err_FCS)
            pass

        if export_data['export_FIT'] == 1:
            try:
                Exportname = save_folder + "/" + filename_i + "_fit_" + start_time + ".csv"
                export_fit = []
                fit = []
                start_force_ss = []
                start_distance_ss = []

                fitting_ds(input_settings, export_data, input_fitting, float(common_steps[0]['step start']), Force_Distance)
                export_fit.append(export_fit_ds)
                print(export_fit)

                if len(common_steps) > 1:
                    for n in range(0, len(common_steps) - 1):
                        fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss = fitting_ss(input_settings, export_data, input_fitting, float(common_steps[n]['step end']), float(common_steps[n + 1]['step start']), Force_Distance, 1, 1)
                        fit.append(fit_ss)
                        start_force_ss.append(f_fitting_region_ss)
                        start_distance_ss.append(d_fitting_region_ss)
                        export_fit.append(export_fit_ss)

                fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss = fitting_ss(input_settings, export_data, input_fitting, float(common_steps[len(common_steps) - 1]['step end']), max(derivation_array[:, 1]), Force_Distance, 1, 1)
                fit.append(fit_ss)
                start_force_ss.append(f_fitting_region_ss)
                start_distance_ss.append(d_fitting_region_ss)
                export_fit.append(export_fit_ss)

                df_export_fit = pd.DataFrame(export_fit)
                df_export_fit.to_csv((Exportname), index=True, header=False)
                plot_fit(fit, start_force_ss, start_distance_ss)
            except:
                print('Something went wrong with fitting')
                pass

        if i == 0:
            print('\nHard work ahead!\n')
            output_q.put('Hard work ahead!')
        elif i == int(len(Files) / 2):
            print('\nHalf way there!\n')
            output_q.put('Half way there!')
            print()
        elif i == len(Files) - 1:
            print('\nAlmost there!\n')
            output_q.put('Almost there!')
        elif i == len(Files):
            print('Analysis finished! \nProgram can be closed.')
            output_q.put('Analysis finished! \nProgram can be closed.')

        total_results.extend(results_F)
        total_results.extend(results_PD)

        i = i + 1
        print('done', i, 'from', len(Files))
        out_progress = str('Done ' + str(i) + ' from ' + str(len(Files)))
        output_q.put(out_progress)

        print(filename_i)
        output_q.put(filename_i)

    if export_data['export_TOTAL'] == 1:
        total_results_export = pd.DataFrame(total_results)
        total_results_export.to_csv((save_folder + '/total results_' + start_time + '.csv'), index=False, header=True)
    else:
        pass


def preprocess_RAW(Force, Distance, input_settings):
    global filteredDistance
    global filteredDistance_ready
    global filteredForce
    global Distance_ds
    global Force_ds

    # Downsample
    Force_ds = Force[::input_settings['downsample_value']]
    Distance_ds = Distance[::input_settings['downsample_value']]

    # Filter
    b, a = signal.butter(input_settings['filter_degree'], input_settings['filter_cut_off'])
    filteredForce = signal.filtfilt(b, a, Force_ds)
    filteredDistance = signal.filtfilt(b, a, Distance_ds)
    filteredDistance_ready = filteredDistance * 1000

    Force_Distance = np.column_stack((filteredForce, filteredDistance_ready))

    if Force_Distance[0, 1] > Force_Distance[-1, 1]:  # reverse
        Force_Distance = np.flipud(Force_Distance)

    Force_Distance_um = np.column_stack((filteredForce, filteredDistance))
    # print(Force_Distance)
    return Force_Distance, Force_Distance_um


# creates a dataset from min force threshold to max force value
def trim_data(FD_data, F_min):
    global F_trimmed
    global PD_trimmed
    global F_max
    global F_low
    global F_up

    F_trimmed = []
    PD_trimmed = []

    F_max = np.where(Force_Distance[:, 0] == max(Force_Distance[:, 0]))
    fi = F_max[0][0]

    while Force_Distance[fi, 0] < Force_Distance[fi - 10, 0]:
        fi = fi - 1

    while Force_Distance[fi, 1] < Force_Distance[fi - 10, 1]:
        fi = fi - 1
    fi0 = fi

    fi = F_max[0][0]
    # print(fi)

    while Force_Distance[fi, 0] > F_min:
        fi = fi - 1
        # print(Force_Distance[fi, 0])
    F_trimmed = FD_data[fi:fi0, 0]
    PD_trimmed = FD_data[fi:fi0, 1]
    F_low = FD_data[:fi, 0]
    F_up = FD_data[fi0:, 0]


# creates derivations for Force and Distance of the trimmed datasets
def create_derivation(input_settings):
    global derivation_array

    step1 = input_settings['downsample_value']

    d_time = 1 / Frequency_value * step1 * input_settings['step_d']

    x = input_settings['step_d']

    derivation_list = []

    while x < len(F_trimmed):
        if PD_trimmed[0] < PD_trimmed[-1]:
            F_value = (F_trimmed[x] + F_trimmed[x - input_settings['step_d']]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[x - input_settings['step_d']]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[x - input_settings['step_d']]
            delta_F = F_trimmed[x] - F_trimmed[x - input_settings['step_d']]
            F_dt = delta_F / d_time
            PD_dt = delta_PD / d_time
        else:
            F_value = (F_trimmed[x] + F_trimmed[(x - input_settings['step_d'])]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[(x - input_settings['step_d'])]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[(x - input_settings['step_d'])]
            delta_F = F_trimmed[x] - F_trimmed[(x - input_settings['step_d'])]
            F_dt = (delta_F / d_time) * (-1)
            PD_dt = (delta_PD / d_time) * (-1)

        derivation_list.append([F_value, PD_value, F_dt, PD_dt])
        x = x + input_settings['step_d']

    derivation_array = np.array(derivation_list)


# calculates the standard deviation of a dataset
def STD(input_data, column_number):
    dt_STD = statistics.pstdev(input_data[:, column_number])
    return dt_STD


# creates a median of the dataset, therefore using slices of len(window_size)
def moving_median(input_data, column_number, window_size):
    # 1st moving median
    window_number = int(len(input_data[:, column_number]) / window_size)  # this is not needed
    window_vector = range(0, (len(input_data[:, column_number]) - window_size), 1)
    mov_med = []

    for n in window_vector:
        mm = np.median(input_data[n:n + window_size, column_number])
        mov_med.append(mm)

    return mov_med


# sorting the data based on a x times STD threshold
def cut_off(input_array, column_number, mm, std, n_of_std):
    global F_dt_inside
    global PD_dt_inside
    global PD_dt_above
    global F_dt_below
    global F_dt_STD2
    global PD_dt_STD2
    global x_STD

    # sort values - inside STD region, above STD region and below STD region
    F_values_inside = []
    PD_values_inside = []
    F_dt_inside = []
    PD_dt_inside = []

    F_values_below = []
    PD_values_below = []
    F_dt_below = []
    PD_dt_below = []

    F_values_above = []
    PD_values_above = []
    F_dt_above = []
    PD_dt_above = []

    i = 0
    for n in range(0, len(input_array[:, column_number]), 1):

        if input_array[n, column_number] > mm[int(i)] + n_of_std * std:
            F_dt_above.append(input_array[n, 2])
            F_values_above.append(input_array[n, 0])
            PD_values_above.append(input_array[n, 1])
            PD_dt_above.append(input_array[n, 3])

        elif input_array[n, column_number] < mm[int(i)] - n_of_std * std:
            F_dt_below.append(input_array[n, 2])
            F_values_below.append(input_array[n, 0])
            PD_values_below.append(input_array[n, 1])
            PD_dt_below.append(input_array[n, 3])
        else:
            F_dt_inside.append(input_array[n, 2])
            F_values_inside.append(input_array[n, 0])
            PD_values_inside.append(input_array[n, 1])
            PD_dt_inside.append(input_array[n, 3])

        i = n * (len(mm) / len(input_array[:, column_number])) - 1

    Above = np.column_stack([F_values_above, PD_values_above, F_dt_above, PD_dt_above])
    Below = np.column_stack([F_values_below, PD_values_below, F_dt_below, PD_dt_below])
    Inside = np.column_stack([F_values_inside, PD_values_inside, F_dt_inside, PD_dt_inside])

    return Above, Inside, Below


# searching for minima in the force derivation to identify unfolding events
def find_steps_F(input_settings, filename_i, Force_Distance, der_arr, results_F, PD_start_F, F_mm2, F_mm2_STD2_positive, F_mm2_STD2_negative):

    STD_1 = STD(der_arr, 2)
    F_mm1 = moving_median(der_arr, 2, input_settings['window_size'])
    Above, Inside, Below = cut_off(der_arr, 2, F_mm1, STD_1, input_settings['x_STD_1'])
    std_0 = STD_1
    Inside_0 = Inside
    F_mm0 = F_mm1
    n_runs = 1

    while abs(std_0 - STD(Inside_0, 2)) / std_0 > input_settings['STD_diff']:
        # print(std_0)
        F_mm0 = moving_median(Inside_0, 2, input_settings['window_size'])
        std_0 = STD(Inside_0, 2)
        Above_0, Inside_0, Below_0 = cut_off(der_arr, 2, F_mm0, std_0, input_settings['x_STD_1'])
        n_runs = n_runs + 1
    STD_2 = std_0
    if STD_2 < 0.05:
        STD_2 = 0.05

    for i in range(len(F_mm0)):
        F_mm2.append(F_mm0[i])

    print('STD is')
    print(STD_2)
    # print(n_runs)

    Above, Inside, Below = cut_off(der_arr, 2, F_mm0, STD_2, input_settings['x_STD_2'])

    for i in range(len(F_mm2)):
        F_mm2_STD2_positive.append(F_mm2[i] + input_settings['x_STD_2'] * STD_2)
        F_mm2_STD2_negative.append(F_mm2[i] - input_settings['x_STD_2'] * STD_2)

    # find the step points
        # for those steps that cross the STD2 threshold -> find the closest 0 values prior/following to the crossing one

    # for local minima

    loc_min = argrelextrema((Below[:, 2]), np.less)
    # print('loc min is')
    # print(loc_min[0])

    n_steps = 1

    for k in loc_min[0]:
        F_dt_loc_min = Below[k, 2]
        F_index = np.where(der_arr[:, 2] == F_dt_loc_min)

        # find start and end of the step
        i_start = F_index[0][0]
        i_end = F_index[0][0]
        while der_arr[i_start, 2] < F_mm2[int(i_start * len(F_mm2) / len(der_arr[:, 2]))] and der_arr[i_start, 2] < der_arr[i_start - 1, 2] and i_start >= 1:
            i_start = i_start - 1
        if i_start == 0:
            i_start = 1
        while der_arr[i_end, 2] < F_mm2[int(i_end * len(F_mm2) / len(der_arr[:, 2]))] and der_arr[i_end, 2] < der_arr[i_end + 1, 2] and i_end<len(der_arr[:, 2])-2:
            i_end = i_end + 1
        
        PD_start_F.append(der_arr[i_start, 1])
        dict1 = {
            "filename": filename_i,
            "Derivation of": 'Force',
            'step #': n_steps,
            'F1': der_arr[i_start, 0],
            'F2': der_arr[i_end, 0],
            'Fc': (der_arr[i_start, 0] + der_arr[i_end, 0]) / 2,
            'step start': der_arr[i_start, 1],
            'step end': der_arr[i_end, 1],
            'step length': der_arr[i_end, 1] - der_arr[i_start, 1],
            'downsample rate': input_settings['downsample_value'],
            'Filter order': input_settings['filter_degree'],
            'Filter frequency': input_settings['filter_cut_off'],
            'STD1': input_settings['x_STD_1'],
            'STD2': input_settings['x_STD_2'],
            'Force min': input_settings['F_min']
            }
        # print(dict1)
        results_F.append(dict1)

        n_steps = n_steps + 1


# searching for maxima in the distance derivation to identify unfolding events
def find_steps_PD(input_settings, filename_i, Force_Distance, der_arr, results_PD, PD_start_PD, PD_mm2, PD_mm2_STD2_positive, PD_mm2_STD2_negative):
    STD_1 = STD(der_arr, 3)
    PD_mm1 = moving_median(der_arr, 3, input_settings['window_size'])

    Above, Inside, Below = cut_off(der_arr, 3, PD_mm1, STD_1, input_settings['x_STD_1'])
    std_0 = STD_1
    Inside_0 = Inside
    PD_mm0 = PD_mm1
    n_runs = 1
    while abs(std_0 - STD(Inside_0, 3)) / std_0 > input_settings['STD_diff']:
        print(std_0)
        PD_mm0 = moving_median(Inside_0, 3, input_settings['window_size'])
        std_0 = STD(Inside_0, 3)
        Above_0, Inside_0, Below_0 = cut_off(der_arr, 3, PD_mm0, std_0, input_settings['x_STD_1'])
        n_runs = n_runs + 1
    STD_2 = std_0
    if STD_2 < 0.05:
        STD_2 = 0.05

    for i in range(len(PD_mm0)):
        PD_mm2.append(PD_mm0[i])
    print('STD is')
    print(STD_2)
    # print(std_0)
    # print(n_runs)
    Above, Inside, Below = cut_off(der_arr, 3, PD_mm0, STD_2, input_settings['x_STD_2'])

    for i in range(len(PD_mm2)):
        PD_mm2_STD2_positive.append(PD_mm2[i] + input_settings['x_STD_2'] * STD_2)
        PD_mm2_STD2_negative.append(PD_mm2[i] - input_settings['x_STD_2'] * STD_2)
    # find the step points
    # for those steps that cross the 3*STD2 threshold -> find the closest 0 values prior/following to the crossing one

    # for local minima

    loc_max = argrelextrema(Above[:, 3], np.greater)
    # print(loc_max[0])

    n_steps = 1

    for k in loc_max[0]:
        PD_dt_loc_max = Above[k, 3]
        PD_index = np.where(der_arr[:, 3] == PD_dt_loc_max)

        # find start and end of the step
        i_start = PD_index[0][0]
        i_end = PD_index[0][0]

        while der_arr[i_start, 3] > PD_mm2[int(i_start * len(PD_mm2) / len(der_arr[:, 3]))] and der_arr[i_start - 1, 3] < der_arr[i_start, 3] and i_start >= 1:
            i_start = i_start - 1
        if i_start == 0:
            i_start = 1

        while der_arr[i_end, 3] > PD_mm2[int(i_end * len(PD_mm2) / len(der_arr[:, 3]))] and der_arr[i_end, 3] > der_arr[i_end + 1, 3] and i_end<len(der_arr[:, 3])-2:
            i_end = i_end + 1

        PD_start_PD.append(der_arr[i_start, 1])

        dict1 = {
            "filename": filename_i,
            "Derivation of": 'Distance',
            'step #': n_steps,
            'F1': der_arr[i_start, 0],
            'F2': der_arr[i_end, 0],
            'Fc': (der_arr[i_start, 0] + der_arr[i_end, 0]) / 2,
            'step start': der_arr[i_start, 1],
            'step end': der_arr[i_end, 1],
            'step length': der_arr[i_end, 1] - der_arr[i_start, 1],
            'downsample rate': input_settings['downsample_value'],
            'Filter order': input_settings['filter_degree'],
            'Filter frequency': input_settings['filter_cut_off'],
            'STD1': input_settings['x_STD_1'],
            'STD2': input_settings['x_STD_2'],
            'Force min': input_settings['F_min']
            }

        results_PD.append(dict1)

        n_steps = n_steps + 1


# define steps, that were found by Force- and Distance-derivation (used for fitting afterwards)
def find_common_steps(F_steps, PD_steps):
    global F_steps_dict
    global PD_steps_dict

    common_steps = []

    for n in range(0, len(F_steps)):
        F_steps_dict = F_steps[n]
        step_F_middle = (float(F_steps_dict['step start']) + float(F_steps_dict['step end'])) / 2
        for i in range(0, len(PD_steps)):
            PD_steps_dict = PD_steps[i]

            if step_F_middle > PD_steps_dict['step start'] and step_F_middle < PD_steps_dict['step end']:
                common_steps.append(PD_steps[i])

    return common_steps


def fitting_ds(input_settings, export_data, input_fitting, i_start, Force_Distance):
    global model_ds, fit_ds
    global ds_fit_dict
    global f_fitting_region_ds, d_fitting_region_ds
    global export_fit_ds

    start_step1 = np.where(derivation_array[:, 1] == i_start)
    start_step1 = start_step1[0][0]

    f_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 0]
    d_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 1]

    model_ds = lk.inverted_odijk("ds_part").subtract_independent_offset() + lk.force_offset("ds_part")

    fit_ds = lk.FdFit(model_ds)

    fit_ds.add_data("Double stranded", f_fitting_region_ds, d_fitting_region_ds)
    # Persistance length bounds
    fit_ds["ds_part/Lp"].value = input_fitting['lp_ds']
    fit_ds["ds_part/Lp"].upper_bound = 80
    fit_ds["ds_part/Lp"].lower_bound = 10

    # Force shift bounds
    fit_ds["ds_part/f_offset"].value = input_fitting['offset_f']
    fit_ds["ds_part/f_offset"].upper_bound = 3
    fit_ds["ds_part/f_offset"].lower_bound = 0
    # distance shift bounds
    fit_ds["ds_part/d_offset"].value = input_fitting['offset_d']
    fit_ds["ds_part/d_offset"].upper_bound = 300
    fit_ds["ds_part/d_offset"].lower_bound = -300
    # stiffnes
    fit_ds["ds_part/St"].value = input_fitting['ds_stiff']
    fit_ds["ds_part/St"].lower_bound = 300
    fit_ds["ds_part/St"].upper_bound = 2000
    # contour length
    Lc_initial_guess = input_fitting['lc_ds']  # nm
    Lc_range = 10
    fit_ds["ds_part/Lc"].upper_bound = Lc_initial_guess + Lc_range
    fit_ds["ds_part/Lc"].lower_bound = Lc_initial_guess - Lc_range
    fit_ds["ds_part/Lc"].value = Lc_initial_guess
    fit_ds["ds_part/Lc"].unit = 'nm'

    fit_ds.fit()
    fit_qual = fit_ds.log_likelihood()
    print(fit_ds.params)

    # export the fitting parameters
    ds_fit_dict = {
        'model': model_ds,
        'fit_model': fit_ds,
        'log_likelihood': fit_qual,
        'Lc': fit_ds["ds_part/Lc"].value,
        'Lp': fit_ds["ds_part/Lp"].value,
        'Lp_stderr': fit_ds["ds_part/Lp"].stderr,
        'St': fit_ds["ds_part/St"].value,
        'f_offset': fit_ds["ds_part/f_offset"].value,
        'd_offset': fit_ds["ds_part/d_offset"].value
        }

    export_fit_ds = pd.DataFrame.from_dict(data=ds_fit_dict, orient='index')
# figure1.savefig(plotname)
# save a figure with multiple plots to visualize the analysis


def fitting_ss(input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range):
    global model_ss
    global ss_fit_dict

    start_fitting_region = np.where(derivation_array[:, 1] == i_start)
    end_fitting_region = np.where(derivation_array[:, 1] == i_end)
    start_fitting_region = start_fitting_region[0][0]
    end_fitting_region = end_fitting_region[0][0]

    raw_f_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 0]
    raw_d_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 1]

    # downsample the data used for fitting to 200 datapoints
    if len(raw_f_fitting_region) > 100:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 100)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 100)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    model_ss = lk.odijk("DNA_2") + lk.freely_jointed_chain("RNA")

    model_ss = model_ss.invert().subtract_independent_offset() + lk.force_offset("DNA")
    fit_ss = lk.FdFit(model_ss)

    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)

    # ds part parameters
    # Persistance length bounds

    # Lp_ds_range=fit_ds["DNA/Lp"].value/10
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp'] * (1 - max_range / 100)
    # if fix==1:
    fit_ss["DNA_2/Lp"].fixed = 'True'

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset']
    fit_ss["DNA/f_offset"].fixed = 'True'

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = 'True'

    # contour length
    # Lc_ds_range=Lc_initial_guess/100 # nm
    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    # if fix==1:
    fit_ss["DNA_2/Lc"].fixed = 'True'

    # stifness

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St']
    if fix == 1:
        fit_ss["DNA_2/St"].fixed = 'True'

    # ss part parameters

    # Persistance length bounds
    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1:
        fit_ss["RNA/Lp"].fixed = 'True'

    # stiffnes
    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = 300
    fit_ss["RNA/St"].upper_bound = 1500

    # contour length
    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss'] + 200

    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    print(fit_ss.params)

    fit_qual = fit_ss.log_likelihood()


    ss_fit_dict = {
        'model': model_ss,
        'fit_model': fit_ss,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ss["DNA_2/Lc"].value,
        'Lp_ds': fit_ss["DNA_2/Lp"].value,
        'St_ds': fit_ss["DNA_2/St"].value,
        'Lc_ss': fit_ss["RNA/Lc"].value,
        'Lc_ss_stderr': fit_ss["RNA/Lc"].stderr,
        'Lp_ss': fit_ss["RNA/Lp"].value,
        'St_ss': fit_ss["RNA/St"].value,
        'f_offset': fit_ss["DNA/f_offset"].value,
        'd_offset': fit_ss["inv(DNA_2_with_RNA)/d_offset"].value
        }

    export_fit_ss = pd.DataFrame.from_dict(data=ss_fit_dict, orient='index')

    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss


def plot_fit(fit, start_force_ss, start_distance_ss):
    distance = np.arange(min(Force_Distance[:, 1]), max(Force_Distance[:, 1]) + 50, 2)
    F_ds_model = model_ds(distance, fit_ds)

    legend_elements = [
        Line2D([0], [0], color='k', lw=1, alpha=0.85),
        Line2D([0], [0], color='r', lw=1),
        Line2D([0], [0], color='gray', linestyle='dashed', lw=1)
        ]

    plt.plot(Force_Distance[:, 1], Force_Distance[:, 0], 'k', alpha=0.85)
    plt.scatter(d_fitting_region_ds, f_fitting_region_ds, color='r', s=4)
    plt.plot(distance, F_ds_model, linestyle='dashed', color='gray')
    plt.ylabel("Force [pN]")
    plt.xlabel("Distance [nm]")
    plt.legend(legend_elements, ['FD-Curve', 'Part used for fitting', 'Fitted WLC model'])

    for i in range(0, len(fit)):
        F_ss_model = model_ss(distance, fit[i])
        plt.scatter(start_distance_ss[i], start_force_ss[i], s=4)
        plt.plot(distance, F_ss_model, linestyle='dashed', color='gray')

    plotname = save_folder + "/" + filename_i + "_fit_" + start_time + ".png"

    plt.savefig(plotname, dpi=600)
    plt.clf()


def save_figure(export_PLOT):
    global Files
    global figure1

    figure1 = Figure(figsize=(10, 6), dpi=100)
    subplot1 = figure1.add_subplot(221)
    subplot2 = figure1.add_subplot(222)
    subplot3 = figure1.add_subplot(223)
    subplot4 = figure1.add_subplot(224)

    legend_elements = [
        Line2D([0], [0], color='b', lw=1),
        Line2D([0], [0], color='r', lw=1),
        Line2D([0], [0], color='k', linestyle='dashed', lw=1)
        ]

    subplot1.set_ylabel("Force (pN)")
    subplot1.scatter(Force_Distance[:, 1], Force_Distance[:, 0], marker='.', s=0.6, linewidths=None, alpha=1)

    subplot2.scatter(PD_trimmed, F_trimmed, marker='.', s=0.6, linewidths=None, alpha=1)

    for i in range(len(PD_start_F)):
        subplot2.axvline(x=PD_start_F[i], ymin=0, ymax=30, color='red', lw=0.5, alpha=0.5)
    for i in range(len(PD_start_PD)):
        subplot2.axvline(x=PD_start_PD[i], ymin=0, ymax=30, color='green', lw=0.5, alpha=0.5)

    subplot3.set_xlabel("Distance (nm)")
    subplot3.set_ylabel("delta Distance (nm/ms)")
    subplot3.plot(derivation_array[1:, 1], derivation_array[1:, 3])
    x_vector = [(max(derivation_array[1:, 1]) - min(derivation_array[1:, 1])) / len(PD_mm2) * g + min(derivation_array[1:, 1]) for g in range(0, len(PD_mm2), 1)]
    subplot3.plot(x_vector, PD_mm2)
    subplot3.fill_between(x_vector, list(PD_mm2_STD2_positive), list(PD_mm2_STD2_negative), color='black', alpha=0.30)

    subplot4.set_xlabel("Distance (nm)")
    subplot4.set_ylabel("delta Force (pN/ms)")
    subplot4.plot(derivation_array[1:, 1], derivation_array[1:, 2])
    x2_vector = [(max(derivation_array[1:, 1]) - min(derivation_array[1:, 1])) / len(F_mm2) * g + min(derivation_array[1:, 1]) for g in range(0, len(F_mm2), 1)]
    subplot4.plot(x2_vector, F_mm2)
    subplot4.fill_between(x2_vector, list(F_mm2_STD2_positive), list(F_mm2_STD2_negative), color='black', alpha=0.30)

    if export_PLOT == 1:
        plotname = save_folder + "/" + filename_i + "_plot_" + start_time + ".png"
        figure1.savefig(plotname, dpi=600)
    else:
        pass

    figure1.clf()
