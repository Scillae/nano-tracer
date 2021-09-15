from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_k2_plot import ns_k2_plot
import os.path


def special_tasks(axs, data, task_list):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    # 1. a horizontal dashed line at 5/32 for all mean
    for i in range(len(arm_num_list)):
        axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
    return axs


def summ_plot_k2(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'k2'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (0, 0.3)
    ylim_std = (0, 0.3)
    ylim_skw = (0.4, 2.0)
    y_var = 'k2 values'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list
    # load data
    summary_dic, savepath = SL(ns_k2_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, special_tasks)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True
