from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_kj_plot import ns_kj_plot
import os.path


def special_tasks(axs, data, task_list):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    # 1. a horizontal dashed line at 5/32 for all mean
    for i in range(len(arm_num_list)):
        axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
    return axs


def summ_plot_kj(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix=''):
    '''
    Summary plot (design fixed) of k2 of the junction.
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'kj'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (0, 0.3)
    ylim_std = (0, 0.3)
    ylim_skw = (0.4, 2.0)
    y_var = 'k2 of Junction'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix
    # load data
    summary_dic, savepath = SL(ns_kj_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, special_tasks, sp_suffix=sp_suffix)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True
