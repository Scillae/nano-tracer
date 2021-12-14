from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_js_plot import ns_js_plot
import os.path


# def special_tasks(axs, data, task_list):
#     conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
#     # 1. a horizontal dashed line at 5/32 for all mean
#     for i in range(len(arm_num_list)):
#         axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
#     return axs


def summ_plot_js(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    '''
    Summary plot (design fixed) of junction shift: distance between CoM of the whole nanostar and the CoM of the junction.
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'js'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (0, 2)
    ylim_std = (0, 2)
    ylim_skw = (0.4, 2.0)
    y_var = 'Junction Shift (nm)'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list
    # load data
    summary_dic, savepath = SL(ns_js_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True