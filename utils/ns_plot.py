from utils.tools import save_load, chkdir
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats
import os.path

def moments_calc(n, var_ls):
    '''
    Calculate the 0th raw, 1st raw, 2nd central, and 3rd standardized moment of a given distribution.
    '''
    n = np.array(n)
    m0 = np.sum(n) # 0th unitless raw moment: integration
    m1 = np.sum(var_ls)/m0
    m2_c = stats.moment(var_ls, moment=2) # 2nd central moment: variance
    std = m2_c**0.5 # standard deviation
    m3_c = stats.moment(var_ls, moment=3) # 3rd central moment
    m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
    return n, m1, std, m3_s


def SL_ns(calc_func, data, varname, vtime=False):
    '''
    Load the calculated result of a measurement.
    If not exists, calculate using ns_func, and then save.
    :calc_func: a calc_func that calculate the desired measurement.
    :data: in which the descriptions of nanostars (trajectory) are stored.
    :varname: codename of the measurement.
    '''
    arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
    traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
    savepath = f'data/composed_traj/{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    plotpath = f'results/{arms}arms{conf_suffix}/{varname}/{varname}_hist-{label}.png'
    if vtime == True:
        plotpath = f'results/{arms}arms{conf_suffix}/{varname}/{varname}_vtime-{label}.png'
    if varname == 'sys' or varname == 'ns':
        var_path = f'{savepath}.{varname}'
    else:
        var_path = f'{savepath}.{varname}tp'
    var_ls_results = save_load(var_path, None)
    if var_ls_results == False:
        var_ls_results = calc_func(top_path, traj_path, arms, dims_ls, f'{savepath}.ns', f'{savepath}.sys')
        var_ls_results = save_load(var_path, var_ls_results)
    return var_ls_results, label, plotpath


def ns_plot(data_process_func, results, plot_confs, data, varname):
    '''
    Plot the histogram (distribution) of a single trajectory.
    Plot: frequency vs value.
    :data_process_func: a function that process the obtained measurements into desired float-valued numbers.
    :results: a dictionary/nested list that contains the obtained measurements
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :varname: codename of the measurement.
    '''
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results

    var_ls = data_process_func(var_ls_results, data)
    if min(var_ls) < x_lim[0] or max(var_ls) > x_lim[1]:
        print(f'ns_plot out of range: {min(var_ls)} < {x_lim[0]} : {min(var_ls) < x_lim[0]} or {max(var_ls)} > {x_lim} : {max(var_ls) > x_lim[1]}')
        # assert 0 == 1
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

    # moments
    n, m1, std, m3_s = moments_calc(n, var_ls)

    ax = plt.axes()
    n_sum = np.sum(n)
    plt.errorbar(bin_centers,np.divide(np.array(n),n_sum),yerr = np.divide(np.sqrt(np.array(n)),n_sum),marker = '.',drawstyle = 'steps-mid',c='#4994FF')
    plt.title(label, fontsize=18)
    # customization
    txt = (
        rf'$N={np.sum(n):.2f}$'
        '\n'
        rf'$\mu={m1:.2f}$'
        '\n'
        rf'$\sigma={std:.2f}$'
        '\n'
        rf'$\gamma={m3_s:.2f}$'
    )
    text_x, text_y = text_pos
    plt.text(text_x, text_y, txt, c='#4994FF', fontsize=18) # 28%
    plt.xlabel(x_var, fontsize=18)
    plt.ylabel('Frequencies', fontsize=18)
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    chkdir(os.path.dirname(plotpath))
    plt.savefig(plotpath)
     # plt.show()
    plt.close()
    plt.clf()
    return m1,std,m3_s

def ns_time_plot(data_process_func, results, plot_confs, data, varname):
    '''
    Plot the value (patch angle vtime) vs. time plot of a single trajectory.
    Plot: value vs time.
    :data_process_func: a function that process the obtained measurements into desired float-valued numbers.
    :results: a dictionary/nested list that contains the obtained measurements
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :varname: codename of the measurement.
    # Arm color schema: Arm_0 ~ #4994FF (blue), Arm_1 ~ #FCC555 (yellow), Arm_2 ~ #7AA77A, Arm_3 ~ #E55050
    # Patch_color_list order: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    Strand-sharing: Arm_0 ~ strands: 1 & 4 ; Arm_1 ~ strands: 4 & 3 ; Arm_2 ~ strands: 3 & 2 ; Arm_3 ~ strands: 2 & 1
    '''
    time_window_width = 25
    assert time_window_width % 2 == 1 # must be odd
    patch_color_list = ['#A3ADAA','#629EBD','#9772A8','#BBB668','#F18B53','#B07C65']
    compare_list = [((0,1),(2,3)),((0,3),(1,2)),((0,2),(1,3)),((0,1),(0,3)),((0,1),(0,2)),((0,3),(0,2))] # former 3 expected to be similar; latter 3 are 3c2=3
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    # only for patch angle vtime.
    assert varname in ['pa','pj']
    var_dic = data_process_func(var_ls_results, data, vtime = True)
    row,col = (len((var_dic.items())),4)
    fig = plt.figure(figsize=(3*3*col+1,3.5*row))
    gs = fig.add_gridspec(row, col, hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
    axs = gs.subplots(sharex='col')
    for i, ((ia1,ia2), ang_vtime_dic) in enumerate(var_dic.items()):
        # i: loop of arm patches.
        ang_ls = [ang_vtime_dic[t][0] for t in ang_vtime_dic if type(t) == int]
        time_ls = [t for t in ang_vtime_dic if type(t) == int]
        time_idx = time_ls[time_window_width//2:-time_window_width//2+1]
        ang_runningavg = [np.sum(ang_ls[j:j+time_window_width-1])/time_window_width for j in range(len(time_idx))]
        ang_runningstd = [np.std(ang_ls[j:j+time_window_width-1]) for j in range(len(time_idx))]
        axs[i,0].plot(time_idx,ang_runningavg, label = f'Arm{ia1}~Arm{ia2}') # ,c=patch_color_list[i]
        axs[i,0].plot([min(time_idx),max(time_idx)],[90,90],c='#000000')
        axs[i,0].plot([min(time_idx),max(time_idx)],[180,180],c='#000000')
        axs[i,0].set_xlabel('Time') #, fontsize=18
        axs[i,0].set_ylabel(f'{x_var} (Avg, wid:{time_window_width})') #, fontsize=18
        axs[i,0].legend(loc = 'lower right')
        axs[i,1].plot(time_idx,ang_runningstd, label = f'Arm{ia1}~Arm{ia2}') # ,c=patch_color_list[i]
        axs[i,1].set_xlabel('Time') #, fontsize=18
        axs[i,1].set_ylabel(f'{x_var} (STD, wid:{time_window_width})') #, fontsize=18
        axs[i,1].legend(loc = 'lower right')
    for i, ((ia1,ia2), ang_vtime_dic) in enumerate(var_dic.items()):
        ang_ls = [ang_vtime_dic[t][0] for t in ang_vtime_dic if type(t) == int]
        time_ls = [t for t in ang_vtime_dic if type(t) == int]
        axs[i,2].scatter(time_ls,ang_ls, label = f'Arm{ia1}~Arm{ia2}') # ,c=patch_color_list[i]
        axs[i,2].plot([min(time_ls),max(time_ls)],[90,90],c='#000000')
        axs[i,2].plot([min(time_ls),max(time_ls)],[180,180],c='#000000')
        axs[i,2].set_xlabel('Time') #, fontsize=18
        axs[i,2].set_ylabel(x_var) #, fontsize=18
        axs[i,2].legend(loc = 'lower right')
    for i, (idx1, idx2) in enumerate(compare_list):
        ang_vtime1 = var_dic[idx1]
        ang_vtime2 = var_dic[idx2]
        time_ls = [t for t in ang_vtime1 if type(t) == int]
        ang_diff_ls = [ang_vtime1[t][0] - ang_vtime2[t][0] for t in time_ls]
        axs[i,3].plot(time_ls,ang_diff_ls, label = f'Plot{idx1} - Plot{idx2}')
        axs[i,3].set_xlabel('Time') #, fontsize=18
        axs[i,3].set_ylabel(f'{x_var}, Diff') #, fontsize=18
        axs[i,3].legend(loc = 'lower right')
    # plt.xlim(x_lim)
    # plt.ylim(y_lim)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    # plotpath = plotpath.split('.')
    # plotpath[-2] = plotpath[-2]+'_vtime'
    # plotpath = '.'.join(plotpath)
    chkdir(os.path.dirname(plotpath))
    plt.savefig(plotpath,dpi=400)
     # plt.show()
    plt.close()
    return