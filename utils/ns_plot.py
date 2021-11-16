from utils.tools import save_load, chkdir
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats
import os.path

def moments_calc(n, var_ls):
    n = np.array(n)
    m0 = np.sum(n) # 0th unitless raw moment: integration
    m1 = np.sum(var_ls)/m0
    m2_c = stats.moment(var_ls, moment=2) # 2nd central moment: variance
    std = m2_c**0.5 # standard deviation
    m3_c = stats.moment(var_ls, moment=3) # 3rd central moment
    m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
    return n, m1, std, m3_s


def SL_ns(calc_func, data, varname):
    arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
    traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
    savepath = f'data/composed_traj/{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    plotpath = f'results/{arms}arms{conf_suffix}/{varname}/{varname}_hist-{label}.png'
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
    plt.clf()
    return m1,std,m3_s