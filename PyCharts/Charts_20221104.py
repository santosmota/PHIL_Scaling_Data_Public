from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

import numpy as np
import pandas as pd
import scipy.io as sio

import matplotlib.pyplot as plt
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['mathtext.fontset'] = 'cm'  # 'cm' Computer modern # 'dejavuserif', 'dejavusans'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'cmr10'  # 'https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/font_file.html

# plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
# ... https://matplotlib.org/stable/tutorials/introductory/customizing.html#matplotlibrc-sample
# DejaVu is shipped with Matplotlib and is thus guaranteed
# to be available; the other entries are
# left as examples of other possible values.
# https://stackoverflow.com/questions/58361594/matplotlib-glyph-8722-missing-from-current-font-despite-being-in-font-manager
# avoids RuntimeWarning: Glyph 8722 missing from current font. font.set_text(s, 0.0, flags=flags)
plt.rc('axes', unicode_minus=False)
#https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels

#############
# solves a warning with a previous syntax
#https://stackoverflow.com/questions/65645194/warning-set-it-to-a-single-string-instead
# plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \renewcommand{\sfdefault}{phv} \renewcommand{\rmdefault}{ptm} \renewcommand{\ttdefault}{pcr}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
from matplotlib.ticker import FormatStrFormatter
# from matplotlib.offsetbox import AnchoredText

# from: https://ctan.uib.no/macros/latex/contrib/IEEEtran/IEEEtran.cls
# % The IEEE uses Times Roman font, so we'll default to Times.
# % These three commands make up the entire times.sty package.
# \renewcommand{\sfdefault}{phv}
# \renewcommand{\rmdefault}{ptm}
# \renewcommand{\ttdefault}{pcr}


###############################################################
# According to the git below, these colors are more colorblind friendly
# https://gist.github.com/thriveth/8560036
###############################################################
color_dalt = {
    'blue':   '#377eb8',        # (55,  126, 184)
    'orange': '#ff7f00',        # (255, 127, 0),
    'green':  '#4daf4a',        # (77,  175, 74),
    'pink':   '#f781bf',        # (247, 129, 191),
    'brown':  '#a65628',        # (166, 86,  40),
    'purple': '#984ea3',        # (152, 78,  163),
    'red':    '#e41a1c',        # (228, 26,  28),
    'yellow': '#dede00',        # (222, 222, 0)
    'gray':  '#999999'         # (153, 153, 153)
}


def dtype_shape_str(x):
    """ Return string containing the dtype and shape of x."""
    return str(x.dtype) + " " + str(x.shape)


# csv_rollingsamples = 50,
# csv_00_time_shift = -0.00355,
# mat_time_shift = -9.0,
# chart_start_time_ms = 40.0,
# chart_end_time_ms = 60.0


#########################################################
# FULL TIME, NO AXIS adjustments
def plot_charts_time_fft(csv_rollingsamples=0,
                         csv_00_time_shift=-0.0036,
                         csv_01_time_shift=0.00105,
                         csv_time_window_s=0.2, # 0.5, # 0.02
                         mat_time_shift=-5.9,
                         csv_pq_mvavg_window_s=0.02,
                         rtl_time_shift=0.0,
                         figure_type='.pdf'
                         ):
    print("#####################")
    print("Function name: ", plot_charts_time_fft.__name__)

    ###############################################################
    # CASES - file names - chart limits - HARDCODED AXIS LIMITS DOWN THERE FOR POWER = 0.75 pu
    ###############################################################
    power_level = 0.8 # 0.75

    if power_level == 0.8:

        # mat_full_path = '../Simulation/version 1/Result/P_0.8_Q_0_5usstep.mat'
        mat_full_path = '../Simulation/version 1/Result/P_0.8_Q_0.375Mvar_5usstep.mat'
        # csv_full_path_00 = '../Experimental Results/20221104/C1-3KHZ-DC125.CSV'
        csv_full_path_00 = '../Experimental Results/20221104/C1.CSV '   # C1-3KHZ-DC11.CSV'
        csv_full_path_01 = '../Experimental Results/20221104/C3.CSV'
        figure_name = 'power_080' + figure_type
        figure_fft_name = 'fft_power_080' + figure_type

        p_min = 3.5
        p_max = 4.5

    else:
        print('CASE NOT VALID!')
        return

    q_min = -1.0
    q_max = 0.5

    # CREATING THE NAME OF THE FIGURE, FROM THE NAME OF THE MAT FILE
    store_in_mat_folder = False

    if store_in_mat_folder:
        index = mat_full_path.rfind('/')  # index to last character '/'
        figure_folder = mat_full_path[0:index + 1]  # selects figure folder
        figure_full_path = figure_folder + figure_name  # adds the .eps or .pdf
        figure_fft_full_path = figure_folder + figure_fft_name  # adds the .eps or .pdf
    else:
        index = csv_full_path_00.rfind('/')  # index to last character '/'
        figure_folder = csv_full_path_00[0:index + 1]  # selects figure folder
        figure_full_path = figure_folder + figure_name  # adds the .eps or .pdf
        figure_fft_full_path = figure_folder + figure_fft_name  # adds the .eps or .pdf

    ###############################################################
    # CASES - scaling between small and large converter
    ###############################################################
    hw_vac = np.array([81.0, 363.0])
    hw_iac = np.array([72.0, 72.0])

    # hw_vac = np.array([363.0, 81.0])
    # hw_iac = np.array([72.0, 72.0])

    hw_s = (3**0.5) * hw_vac * hw_iac

    sw_vac = 690.0
    sw_s = 5e6
    sw_iac = sw_s / (3**0.5 * sw_vac)

    ###############################################################
    # OSCILLOSCOPE - PROBE SCALING factors
    ###############################################################
    an_dig_raw_full_scale = np.array([2048.0, 2048.0, 2048.0, 2048.0])    # analog to digital full scale, 12 bit, +-2048
    conv_bnc_scale = np.array([1000.0, 1000.0, 1000.0, 1000.0])           # converter bnc scaling page 23 DA A scale: 1000 / 1V
    meas_full_scale = np.array([600.0, 500.0, 600.0, 500.0])              # voltage, current,app. power, app. power

    ###############################################################
    # OSCILLOSCOPE - DATA FILE - CSVs
    ###############################################################
    print('###########################')
    print('CSV ')
    print('  first file:  ', csv_full_path_00)
    # print('  second file: ', csv_full_path_01)

    # read csv files and assign column names
    csv_df_00 = pd.read_csv(csv_full_path_00, skiprows=[0], header=None)
    csv_df_00.columns = ['t', 'Va', 'Ia', 'Vb', 'Ib']
    csv_df_01 = pd.read_csv(csv_full_path_01, skiprows=[0], header=None)
    csv_df_01.columns = ['t', 'Va', 'Ia', 'Vb', 'Ib']

    # slice csv files to a time window of "csv_time_window_s"
    csv_delta_t = (csv_df_00['t'].iloc[-1] - csv_df_00['t'].iloc[0]) / csv_df_00.shape[0]
    csv_nsamples = int(csv_time_window_s / csv_delta_t) + 1

    tot = csv_df_00.shape[0]        # total number of samples
    aux = int(csv_nsamples / 2)     # half the desired number of samples
    mid = int(tot/2)                # index for the middle of the time (where zero is, +-1 or so)
    start = mid - aux               # start index
    end = mid + aux                 # end index

    csv_df_00 = csv_df_00.iloc[start:end, :]
    csv_df_01 = csv_df_01.iloc[start:end, :]

    # recalculate csv_nsamples, lazy way of solving +-1 errors in the nsamples
    csv_nsamples = csv_df_00.shape[0]

    print('Time resolution in the file:')
    print('    ', csv_full_path_00)
    print('        ', float(csv_delta_t), 's')

    print('Time window:')
    print('    start:', csv_df_00['t'].iloc[0], 's')
    print('    end:', csv_df_00['t'].iloc[-1], 's')
    print('    total:', csv_df_00['t'].iloc[-1] - csv_df_00['t'].iloc[0], 's')

    # scaling the scope to actual voltage, current
    for col, i in zip(csv_df_00.columns, range(0, 5)):
        if i > 0:
            csv_df_00[col] = csv_df_00[col] * conv_bnc_scale[i-1] * meas_full_scale[i-1] / an_dig_raw_full_scale[i-1]

    for col, i in zip(csv_df_01.columns, range(0, 5)):
        if i > 0:
            csv_df_01[col] = csv_df_01[col] * conv_bnc_scale[i-1] * meas_full_scale[i-1] / an_dig_raw_full_scale[i-1]

    # scaling the scope to the LARGE converter
    csv_df_00['Va'] = csv_df_00['Va'] * sw_vac / hw_vac[0]
    csv_df_00['Ia'] = csv_df_00['Ia'] * sw_iac / hw_iac[0]
    csv_df_00['Vb'] = csv_df_00['Vb'] * sw_vac / hw_vac[0]
    csv_df_00['Ib'] = csv_df_00['Ib'] * sw_iac / hw_iac[0]

    csv_df_01['Va'] = csv_df_01['Va'] * sw_vac / hw_vac[1]
    csv_df_01['Ia'] = csv_df_01['Ia'] * sw_iac / hw_iac[1]
    csv_df_01['Vb'] = csv_df_01['Vb'] * sw_vac / hw_vac[1]
    csv_df_01['Ib'] = csv_df_01['Ib'] * sw_iac / hw_iac[1]

    # dataframe time column
    csv_df_00['time'] = (csv_df_00['t'] + csv_00_time_shift) * 1000.0
    csv_df_01['time'] = (csv_df_01['t'] + csv_01_time_shift) * 1000.0

    csv_total_time_ms = csv_df_00['time'].iloc[-1] - csv_df_00['time'].iloc[0]
    csv_total_time_s = csv_df_00['t'].iloc[-1] - csv_df_00['t'].iloc[0]

    ###############################################################
    #Vc Ic
    ###############################################################
    csv_df_00['Vc'] = - csv_df_00['Va'] - csv_df_00['Vb']
    csv_df_00['Ic'] = - csv_df_00['Ia'] - csv_df_00['Ib']

    csv_df_00['Valpha'] = (2.0 * csv_df_00['Va'] - csv_df_00['Vb'] - csv_df_00['Vc']) / 3.0
    csv_df_00['Vbeta'] = (csv_df_00['Vb'] - csv_df_00['Vc']) / (3.0**0.5)

    csv_df_00['Ialpha'] = (2.0 * csv_df_00['Ia'] - csv_df_00['Ib'] - csv_df_00['Ic']) / 3.0
    csv_df_00['Ibeta'] = (csv_df_00['Ib'] - csv_df_00['Ic']) / (3.0 ** 0.5)

    csv_df_00['P'] = 3 * (csv_df_00['Valpha'] * csv_df_00['Ialpha'] + csv_df_00['Vbeta'] * csv_df_00['Ibeta']) / 2.0
    csv_df_00['Q'] = 3 * (csv_df_00['Vbeta'] * csv_df_00['Ialpha'] - csv_df_00['Valpha'] * csv_df_00['Ibeta']) / 2.0

    csv_df_01['Vc'] = - csv_df_01['Va'] - csv_df_01['Vb']
    csv_df_01['Ic'] = - csv_df_01['Ia'] - csv_df_01['Ib']

    csv_df_01['Valpha'] = (2.0 * csv_df_01['Va'] - csv_df_01['Vb'] - csv_df_01['Vc']) / 3.0
    csv_df_01['Vbeta'] = (csv_df_01['Vb'] - csv_df_01['Vc']) / (3.0**0.5)

    csv_df_01['Ialpha'] = (2.0 * csv_df_01['Ia'] - csv_df_01['Ib'] - csv_df_01['Ic']) / 3.0
    csv_df_01['Ibeta'] = (csv_df_01['Ib'] - csv_df_01['Ic']) / (3.0 ** 0.5)

    csv_df_01['P'] = 3 * (csv_df_01['Valpha'] * csv_df_01['Ialpha'] + csv_df_01['Vbeta'] * csv_df_01['Ibeta']) / 2.0
    csv_df_01['Q'] = 3 * (csv_df_01['Vbeta'] * csv_df_01['Ialpha'] - csv_df_01['Valpha'] * csv_df_01['Ibeta']) / 2.0

    ###############################################################
    # DATA FILE - mat
    ###############################################################

    print('###########################')
    print('MAT ')
    print('  file:  ', mat_full_path)
    mat = sio.loadmat(mat_full_path)

    ###############
    ## hardcoded stuff to understand the structure of the mat file
    # print("mat['Measurements'] is a:", dtype_shape_str(mat['Measurements']))
    # # print("mat['Measurements'][0, 0] is a:", dtype_shape_str(mat['Measurements'][0, 0]))
    # print("mat['Measurements']['time'] is a:", dtype_shape_str(mat['Measurements']['time']))
    # print("mat['Measurements']['time'][0, 0] is a:", dtype_shape_str(mat['Measurements']['time'][0, 0]))
    # print("mat['Measurements'] is a:",
    #       dtype_shape_str(mat['Measurements']))
    # print("mat['Measurements']['signals'] is a:",
    #       dtype_shape_str(mat['Measurements']['signals']))
    # print("mat['Measurements']['signals'][0, 0] is a:",
    #       dtype_shape_str(mat['Measurements']['signals'][0, 0]))
    # print("mat['Measurements']['signals'][0, 0]['values'] is a:",
    #       dtype_shape_str(mat['Measurements']['signals'][0, 0]['values']))
    # print("mat['Measurements']['signals'][0, 0]['values'][0, 0] is a:",
    #       dtype_shape_str(mat['Measurements']['signals'][0, 0]['values'][0, 0]))
    # print("mat['Measurements']['signals'][0, 0]['label'] is a:",
    #       dtype_shape_str(mat['Measurements']['signals'][0, 0]['label']))
    # print("mat['Measurements']['signals'][0, 0]['label'][0, 0] is a:",
    #       dtype_shape_str(mat['Measurements']['signals'][0, 0]['label'][0, 0]))

    mat_df = pd.DataFrame(mat['Measurements']['signals'][0, 0]['values'][0, 0],
                             columns=['Va', 'Vb', 'Vc', 'Ia', 'Ib', 'Ic', 'P', 'Q'])

    mat_df['time'] = (np.array(mat['Measurements']['time'][0, 0]) + mat_time_shift) * 1000.0

    mat_deltat = (mat_df['time'].iloc[-1] - mat_df['time'].iloc[0]) / mat_df.shape[0]

    print('Time resolutions in the file:')
    print('    ', mat_full_path)
    print('        ', float(mat_deltat)/1000.0, 's')


    ##############################################################
    # size of the figure (in inches), and linestyles - FFTs
    ##############################################################
    figsizex = 4.5 # 9
    figsizey = 3.5
    fig_fft, axes_fft = plt.subplots(2, 1, sharex=True,
                                     figsize=(figsizex, figsizey),
                                     num='Freq')

    #########################################################################
    # PLOT hardcoded colors and styles
    #########################################################################

    # mat_color = color_dalt['blue']
    csv_color_00 = color_dalt['green']
    csv_color_01 = color_dalt['red']

    mat_color = 'black'

    mat_label = 'FSC'
    csv_00_label = 'SDC Case 1'
    csv_01_label = 'SDC Case 2'

    mat_lin_style = 'solid'
    csv_00_lin_style = 'solid'
    csv_01_lin_style = 'solid'

    mat_lin_width = 0.75
    csv_lin_width = 0.5

    ##########################################################################
    # ffts
    ##########################################################################
    ###################
    # Slice mat file

    strt = mat_df.index[mat_df['time'] == 0.0].tolist()
    strt = strt[0]
    end = strt + int(csv_total_time_ms / mat_deltat)

    sig_mat_Va = mat_df.loc[strt:end, 'Va'].to_numpy()
    sig_mat_Ia = mat_df.loc[strt:end, 'Ia'].to_numpy()
    # sig_mat_P = mat_df.loc[strt:end, 'P'].to_numpy()
    # sig_mat_Q = mat_df.loc[strt:end, 'Q'].to_numpy()
    mat_slice_nsamples = sig_mat_Ia.shape[-1]

    print('Time window:')
    print('    start:', mat_df['time'].iloc[strt] / 1000.0, 's')
    print('    end:', mat_df['time'].iloc[strt + mat_slice_nsamples - 1] / 1000.0, 's')
    print('    total:', (mat_df['time'].iloc[strt + mat_slice_nsamples - 1] - mat_df['time'].iloc[strt]) / 1000.0, 's')

    ###################
    # Frequency array
    if mat_slice_nsamples % 2 == 0:  # even
        freq_mat = np.arange(start=0, stop=(mat_slice_nsamples/2 + 1), step=1) / (mat_deltat / 1000.0 * mat_slice_nsamples)
    else:  # odd
        freq_mat = np.arange(start=0, stop=((mat_slice_nsamples + 1) / 2), step=1) / (mat_deltat / 1000.0 * mat_slice_nsamples)

    if csv_nsamples % 2 == 0:  # even
        freq_csv = np.arange(start=0, stop=(csv_nsamples/2 + 1), step=1) / (csv_total_time_s)
    else:  # odd
        freq_csv = np.arange(start=0, stop=((csv_nsamples + 1) / 2), step=1) / (csv_total_time_s)

    fft_mat_Va = np.fft.rfft(sig_mat_Va) * 2 / mat_slice_nsamples
    fft_mat_Ia = np.fft.rfft(sig_mat_Ia) * 2 / mat_slice_nsamples / 1000.0
    # fft_mat_P = np.fft.rfft(sig_mat_P) * 2 / mat_slice_nsamples / 1000000.0
    # fft_mat_Q = np.fft.rfft(sig_mat_Q) * 2 / mat_slice_nsamples / 1000000.0

    fft_mat_Va[0] = fft_mat_Va[0] / 2.0
    fft_mat_Ia[0] = fft_mat_Ia[0] / 2.0
    # fft_mat_P[0] = fft_mat_P[0] / 2.0
    # fft_mat_Q[0] = fft_mat_Q[0] / 2.0

    fft_csv_00_Va = np.fft.rfft(csv_df_00['Va']) * 2 / csv_nsamples
    fft_csv_00_Ia = np.fft.rfft(csv_df_00['Ia']) * 2 / csv_nsamples / 1000.0
    # fft_csv_00_P = np.fft.rfft(csv_df_00['P']) * 2 / csv_nsamples / 1000000.0
    # fft_csv_00_Q = np.fft.rfft(csv_df_00['Q']) * 2 / csv_nsamples / 1000000.0

    fft_csv_00_Va[0] = fft_csv_00_Va[0] / 2.0
    fft_csv_00_Ia[0] = fft_csv_00_Ia[0] / 2.0
    # fft_csv_00_P[0] = fft_csv_00_P[0] / 2.0
    # fft_csv_00_Q[0] = fft_csv_00_Q[0] / 2.0

    fft_csv_01_Va = np.fft.rfft(csv_df_01['Va']) * 2 / csv_nsamples
    fft_csv_01_Ia = np.fft.rfft(csv_df_01['Ia']) * 2 / csv_nsamples / 1000.0
    # fft_csv_01_P = np.fft.rfft(csv_df_01['P']) * 2 / csv_nsamples / 1000000.0
    # fft_csv_01_Q = np.fft.rfft(csv_df_01['Q']) * 2 / csv_nsamples / 1000000.0

    fft_csv_01_Va[0] = fft_csv_01_Va[0] / 2.0
    fft_csv_01_Ia[0] = fft_csv_01_Ia[0] / 2.0
    # fft_csv_01_P[0] = fft_csv_01_P[0] / 2.0
    # fft_csv_01_Q[0] = fft_csv_01_Q[0] / 2.0

    ##################
    # DATA FRAME WITH FFTs

    fft_csv_00_df = pd.DataFrame(np.transpose(np.array([freq_csv, np.absolute(fft_csv_00_Va), np.absolute(fft_csv_00_Ia)])),
                                 columns=['freq', 'Va', 'Ia'])

    fft_csv_01_df = pd.DataFrame(np.transpose(np.array([freq_csv, np.absolute(fft_csv_01_Va), np.absolute(fft_csv_01_Ia)])),
                                 columns=['freq', 'Va', 'Ia'])

    fft_mat_df = pd.DataFrame(np.transpose(np.array([freq_mat, np.absolute(fft_mat_Va), np.absolute(fft_mat_Ia)])),
                              columns=['freq', 'Va', 'Ia'])

    ######################
    # FFT PLOTS
    # voltages
    axes_fft[0].plot(fft_csv_00_df['freq'],
                     fft_csv_00_df['Va'],
                     color=csv_color_00,
                     linewidth=csv_lin_width,
                     linestyle=csv_00_lin_style,
                     label=csv_00_label)

    axes_fft[0].plot(fft_csv_01_df['freq'],
                     fft_csv_01_df['Va'],
                     color=csv_color_01,
                     linewidth=csv_lin_width,
                     linestyle=csv_01_lin_style,
                     label=csv_01_label)

    axes_fft[0].plot(fft_mat_df['freq'],
                     fft_mat_df['Va'],
                     color=mat_color,
                     linewidth=mat_lin_width,
                     linestyle=mat_lin_style,
                     label=mat_label)

    ax_00_ins_fft = axes_fft[0].inset_axes([0.15, 0.3, 0.5, 0.6])
    x1, x2, y1, y2 = 2750.0, 3250.0, -0.10, 5.0

    ax_00_ins_fft.plot(fft_csv_00_df['freq'],
                       fft_csv_00_df['Va'],
                       color=csv_color_00,
                       linewidth=csv_lin_width,
                       linestyle=csv_00_lin_style,
                       label=csv_00_label)

    ax_00_ins_fft.plot(fft_csv_01_df['freq'],
                       fft_csv_01_df['Va'],
                       color=csv_color_01,
                       linewidth=csv_lin_width,
                       linestyle=csv_01_lin_style,
                       label=csv_01_label)

    ax_00_ins_fft.plot(fft_mat_df['freq'],
                       fft_mat_df['Va'],
                       color=mat_color,
                       linewidth=mat_lin_width,
                       linestyle=mat_lin_style,
                       label=mat_label)

    # ax_10_ins_fft.set_xscale('symlog',
    #                          linthresh=100,
    #                          subs=np.arange(2, 10),
    #                          linscale=0.35)

    ax_00_ins_fft.set_xlim(x1, x2)
    ax_00_ins_fft.set_ylim(y1, y2)
    # ax_10_ins_fft.set_xticklabels([])

    axes_fft[0].indicate_inset_zoom(ax_00_ins_fft, edgecolor="black")

    ######################
    # FFT
    # currents
    axes_fft[1].plot(fft_csv_00_df['freq'],
                     fft_csv_00_df['Ia'],
                     color=csv_color_00,
                     linewidth=csv_lin_width,
                     linestyle=csv_00_lin_style,
                     label=csv_00_label)

    axes_fft[1].plot(fft_csv_01_df['freq'],
                     fft_csv_01_df['Ia'],
                     color=csv_color_01,
                     linewidth=csv_lin_width,
                     linestyle=csv_01_lin_style,
                     label=csv_01_label)

    axes_fft[1].plot(fft_mat_df['freq'],
                     fft_mat_df['Ia'],
                     color=mat_color,
                     linewidth=mat_lin_width,
                     linestyle=mat_lin_style,
                     label=mat_label)

    ax_10_ins_fft = axes_fft[1].inset_axes([0.2, 0.3, 0.75, 0.6])
    x1, x2, y1, y2 = 2750.0, 3250, -0.025, 0.3

    ax_10_ins_fft.plot(fft_csv_00_df['freq'],
                       fft_csv_00_df['Ia'],
                       color=csv_color_00,
                       linewidth=csv_lin_width,
                       linestyle=csv_00_lin_style,
                       label=csv_00_label)

    ax_10_ins_fft.plot(fft_csv_01_df['freq'],
                       fft_csv_01_df['Ia'],
                       color=csv_color_01,
                       linewidth=csv_lin_width,
                       linestyle=csv_01_lin_style,
                       label=csv_01_label)

    ax_10_ins_fft.plot(fft_mat_df['freq'],
                       fft_mat_df['Ia'],
                       color=mat_color,
                       linewidth=mat_lin_width,
                       linestyle=mat_lin_style,
                       label=mat_label)

    # ax_10_ins_fft.set_xscale('symlog',
    #                          linthresh=100,
    #                          subs=np.arange(2, 10),
    #                          linscale=0.35)

    ax_10_ins_fft.set_xlim(x1, x2)
    ax_10_ins_fft.set_ylim(y1, y2)
    # ax_10_ins_fft.set_xticklabels([])

    # ax_10_ins_fft.legend(loc='best', frameon=False, prop={'size': 9})

    axes_fft[1].indicate_inset_zoom(ax_10_ins_fft, edgecolor="black")

    ##########################################################################
    # Grouping the interharmonics
    ##########################################################################
    # ffts_grouped_df = pd.DataFrame( {'h': range(1, 101)})



    # for i in range()

    ###############################################################
    # size of the figure (in inches), and linestyles - TIME DOMAIN
    ###############################################################
    figsizex = 9
    figsizey = 3.5
    fig, axes = plt.subplots(2, 2, sharex=True,
                             figsize=(figsizex, figsizey),
                             num='Charts')

    #######################
    # moving average
    if csv_pq_mvavg_window_s > 0.0:
        n_mvag = int(csv_pq_mvavg_window_s / csv_delta_t)

        csv_df_00['P'] = csv_df_00['P'].rolling(n_mvag).mean()
        csv_df_00['Q'] = csv_df_00['Q'].rolling(n_mvag).mean()

        csv_df_01['P'] = csv_df_01['P'].rolling(n_mvag).mean()
        csv_df_01['Q'] = csv_df_01['Q'].rolling(n_mvag).mean()

        n_mvag = int(csv_pq_mvavg_window_s *1000.0 / mat_deltat)

        mat_df['P'] = mat_df['P'].rolling(n_mvag).mean()

        mat_df['Q'] = mat_df['Q'].rolling(n_mvag).mean()


    #######################
    # Voltages
    # smoothing oscilloscope data
    if csv_rollingsamples > 1:
        csv_df_00['Va'] = csv_df_00['Va'].rolling(csv_rollingsamples).mean()
        csv_df_00['Ia'] = csv_df_00['Ia'].rolling(csv_rollingsamples).mean()
        # csv_df_00['P'] = csv_df_00['P'].rolling(csv_rollingsamples).mean()
        # csv_df_00['Q'] = csv_df_00['Q'].rolling(csv_rollingsamples).mean()

        csv_df_01['Va'] = csv_df_01['Va'].rolling(csv_rollingsamples).mean()
        csv_df_01['Ia'] = csv_df_01['Ia'].rolling(csv_rollingsamples).mean()
        # csv_df_01['P'] = csv_df_01['P'].rolling(csv_rollingsamples).mean()
        # csv_df_01['Q'] = csv_df_01['Q'].rolling(csv_rollingsamples).mean()
        #

    #######################
    # Voltages
    axes[0, 0].plot(csv_df_00['time'],
                    csv_df_00['Va'],
                    color=csv_color_00,
                    linewidth=csv_lin_width,
                    linestyle=csv_00_lin_style,
                    label=csv_00_label)

    axes[0, 0].plot(csv_df_01['time'],
                    csv_df_01['Va'],
                    color=csv_color_01,
                    linewidth=csv_lin_width,
                    linestyle=csv_01_lin_style,
                    label=csv_01_label)

    axes[0, 0].plot(mat_df['time'],
                    mat_df['Va'],
                    color=mat_color,
                    linewidth=mat_lin_width,
                    linestyle=mat_lin_style,
                    label=mat_label)

    ax_00_ins = axes[0, 0].inset_axes([0.19, 0.075, 0.28, 0.5])
    x1, x2, y1, y2 = 52.0, 54.0, -500.0, -100.0

    ax_00_ins.plot(csv_df_00['time'],
                   csv_df_00['Va'],
                   color=csv_color_00,
                   linewidth=csv_lin_width,
                   linestyle=csv_00_lin_style,
                   label=csv_00_label)

    ax_00_ins.plot(csv_df_01['time'],
                   csv_df_01['Va'],
                   color=csv_color_01,
                   linewidth=csv_lin_width,
                   linestyle=csv_01_lin_style,
                   label=csv_01_label)

    ax_00_ins.plot(mat_df['time'],
                   mat_df['Va'],
                   color=mat_color,
                   linewidth=mat_lin_width,
                   linestyle=mat_lin_style,
                   label=mat_label)

    ax_00_ins.set_xlim(x1, x2)
    ax_00_ins.set_ylim(y1, y2)
    ax_00_ins.set_xticklabels([])
    # ax_00_ins.set_yticklabels([])
    axes[0, 0].indicate_inset_zoom(ax_00_ins, edgecolor="black")

    #######################
    # Currents
    axes[1, 0].plot(csv_df_00['time'],
                    csv_df_00['Ia']/1e3,
                    color=csv_color_00,
                    linewidth=csv_lin_width,
                    linestyle=csv_00_lin_style,
                    label=csv_00_label)

    axes[1, 0].plot(csv_df_01['time'],
                    csv_df_01['Ia'] / 1e3,
                    color=csv_color_01,
                    linewidth=csv_lin_width,
                    linestyle=csv_01_lin_style,
                    label=csv_01_label)

    axes[1, 0].plot(mat_df['time'],
                    mat_df['Ia']/1e3,
                    color=mat_color,
                    linewidth=mat_lin_width,
                    linestyle=mat_lin_style,
                    label=mat_label)

    ax_10_ins = axes[1, 0].inset_axes([0.19, 0.05, 0.28, 0.5])
    x1, x2, y1, y2 = 50.0, 52.0, -2.5, 2.5

    ax_10_ins.plot(csv_df_00['time'],
                   csv_df_00['Ia']/1e3,
                   color=csv_color_00,
                   linewidth=csv_lin_width,
                   linestyle=csv_00_lin_style,
                   label=csv_00_label)

    ax_10_ins.plot(csv_df_01['time'],
                   csv_df_01['Ia'] / 1e3,
                   color=csv_color_01,
                   linewidth=csv_lin_width,
                   linestyle=csv_01_lin_style,
                   label=csv_01_label)

    ax_10_ins.plot(mat_df['time'],
                   mat_df['Ia']/1e3,
                   color=mat_color,
                   linewidth=mat_lin_width,
                   linestyle=mat_lin_style,
                   label=mat_label)

    ax_10_ins.set_xlim(x1, x2)
    ax_10_ins.set_ylim(y1, y2)
    ax_10_ins.set_xticklabels([])
    # ax_10_ins.set_yticklabels([])
    axes[1, 0].indicate_inset_zoom(ax_10_ins, edgecolor="black")

    #######################
    # Active power
    axes[0, 1].plot(csv_df_00['time'],
                    csv_df_00['P']/1e6,
                    color=csv_color_00,
                    linewidth=csv_lin_width,
                    linestyle=csv_00_lin_style,
                    label=csv_00_label)

    axes[0, 1].plot(csv_df_01['time'],
                    csv_df_01['P'] / 1e6,
                    color=csv_color_01,
                    linewidth=csv_lin_width,
                    linestyle=csv_01_lin_style,
                    label=csv_01_label)

    axes[0, 1].plot(mat_df['time'],
                    mat_df['P']/1e6,
                    color=mat_color,
                    linewidth=mat_lin_width,
                    linestyle=mat_lin_style,
                    label=mat_label)
    #
    # #######################
    # # Reactive power, not done yet
    axes[1, 1].plot(csv_df_00['time'],
                    csv_df_00['Q']/1e6,
                    color=csv_color_00,
                    linewidth=csv_lin_width,
                    linestyle=csv_00_lin_style,
                    label='SDC Case 1 Experimental')  # csv_00_label)

    axes[1, 1].plot(csv_df_01['time'],
                    csv_df_01['Q'] / 1e6,
                    color=csv_color_01,
                    linewidth=csv_lin_width,
                    linestyle=csv_01_lin_style,
                    label='SDC Case 2 Experimental')  # csv_01_label)

    axes[1, 1].plot(mat_df['time'],
                    mat_df['Q']/1e6,
                    color=mat_color,
                    linewidth=mat_lin_width,
                    linestyle=mat_lin_style,
                    label='FSC Simulation')  # mat_label)


    ##########################################################################
    # axis limits
    ##########################################################################
    axes[1, 0].set_xticks(np.arange(40, 80, 2))
    axes[1, 0].set_xlim([40, 60])
    # axes[1, 0].set_xlim([0, 60])

    axes[0, 0].set_yticks(np.arange(-600.0, 800.0, 200.0))
    axes[0, 0].set_ylim([-650, 650])
    # axes[0][0].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    axes[1, 0].set_yticks(np.arange(-6.0, 8.0, 2.0))
    axes[1, 0].set_ylim([-6.0, 6.0])

    axes[0, 1].set_yticks(np.arange(3.5, 4.5, 0.1))
    axes[0, 1].set_ylim([3.8, 4.2])
    # axes[0, 1].set_ylim([2, 6])

    axes[1, 1].set_yticks(np.arange(-1.0, 1.0, 0.25))
    axes[1, 1].set_ylim([-0.75, 0.25])

    # axes_fft[1, 0].set_xticks(np.arange(40, 80, 2))
    # axes_fft[1, 0].set_xlim([0, 8000])

    axes_fft[1].set_xscale('symlog',
                              linthresh=100,
                              subs=np.arange(2, 10),
                              linscale=0.35)
    axes_fft[1].set_xlim([0, 10000])
    axes_fft[1].set_xticks([0, 50, 100, 1000, 10000])
    axes_fft[1].set_xticklabels([0, 50, 100, 1000, 10000])

    # axes_fft[1, 0].set_yscale('log')
    # axes_fft[1, 0].set_ylim([1, 5000])

    ##########################################################################
    # axis names
    ##########################################################################
    axes[1][0].set_xlabel(r'Time (ms)')
    axes[1][1].set_xlabel(r'Time (ms)')

    axes[0][0].set_ylabel(r'Voltage (V)')
    axes[1][0].set_ylabel(r'Current (kA)')
    axes[0][1].set_ylabel(r'Active Power (MW)')
    axes[1][1].set_ylabel(r'Reactive (Mvar)')

    axes_fft[1].set_xlabel(r'Time (ms)')
    # axes_fft[1][1].set_xlabel(r'Time (ms)')

    axes_fft[00].set_ylabel(r'Voltage (V)')
    axes_fft[1].set_ylabel(r'Current (kA)')
    # axes_fft[0][1].set_ylabel(r'Active Power (MW)')
    # axes_fft[1][1].set_ylabel(r'Reactive (Mvar)')

    ##########################################################################
    # chart identification - legend - abcdefghi
    ##########################################################################
    # https://matplotlib.org/stable/gallery/color/named_colors.html
    # colors lightgray gray aliceblue whitesmoke
    corlegenda = 'whitesmoke'
    #
    axes[0, 0].annotate(r'a', xy=(0.05, 0.85), xycoords='axes fraction',
                        bbox=dict(boxstyle='circle', fc=corlegenda))

    axes[0, 1].annotate(r'b', xy=(0.05, 0.85), xycoords='axes fraction',
                        bbox=dict(boxstyle='circle', fc=corlegenda))

    axes[1, 0].annotate(r'c', xy=(0.05, 0.85), xycoords='axes fraction',
                        bbox=dict(boxstyle='circle', fc=corlegenda))

    axes[1, 1].annotate(r'd', xy=(0.05, 0.85), xycoords='axes fraction',
                        bbox=dict(boxstyle='circle', fc=corlegenda))

    axes_fft[0].annotate(r'a', xy=(0.025, 0.85), xycoords='axes fraction',
                         bbox=dict(boxstyle='circle', fc=corlegenda))

    axes_fft[1].annotate(r'b', xy=(0.025, 0.85), xycoords='axes fraction',
                         bbox=dict(boxstyle='circle', fc=corlegenda))

    ##########################################################################
    # axis legends
    ##########################################################################
    # axes[0, 0].legend(loc='best', frameon=False, prop={'size': 9})
    # # axes[1, 0].legend(loc='best', frameon=False, prop={'size': 9})
    # # axes[0, 1].legend(loc='best', frameon=False, prop={'size': 9})
    axes[1, 1].legend(loc='best', frameon=False, prop={'size': 9})
    #
    axes_fft[0].legend(loc='best', frameon=False, prop={'size': 9})

    ##########################################################################
    # align, tighten, shown and save
    ##########################################################################
    fig.align_ylabels(axes[:, :])
    fig.tight_layout()
    # fig.show()

    # fig_fft.align_ylabels(axes_fft[:, :])
    fig_fft.align_ylabels(axes_fft[:])
    fig_fft.tight_layout()

    if figure_type == '.pdf':
        fig.savefig(figure_full_path, format="pdf", bbox_inches="tight")
        fig_fft.savefig(figure_fft_full_path, format="pdf", bbox_inches="tight")
    elif figure_type == '.eps':
        fig.savefig(figure_full_path, format='eps')
        fig_fft.savefig(figure_fft_full_path, format='eps')

    if figure_type == '.pdf' or figure_type == '.eps':
        # fft_csv_00_df.to_csv(figure_folder + 'fft_csv_00_20ms.csv')
        # fft_csv_01_df.to_csv(figure_folder + 'fft_csv_01_20ms.csv')
        # fft_mat_df.to_csv(figure_folder + 'fft_mat_20ms.csv')
        fft_csv_00_df.to_csv(figure_folder + 'fft_csv_00.csv')
        fft_csv_01_df.to_csv(figure_folder + 'fft_csv_01.csv')
        fft_mat_df.to_csv(figure_folder + 'fft_mat.csv')

    plt.show()



def main():
    print("#####################")
    print("Function name: ", main.__name__)

    figure_type = '.eps'
    # figure_type = '.pdf'
    # figure_type = 'none'


    plot_charts_time_fft(figure_type=figure_type)


if __name__ == '__main__':
    main()
