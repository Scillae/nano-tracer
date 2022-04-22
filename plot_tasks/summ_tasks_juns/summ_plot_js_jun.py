from utils.tools import chkdir
from utils.summ_plot import summ_plot_jun, SL_jun
from plot_tasks.ns_plots.ns_js_plot import ns_js_plot
import os.path


# def special_tasks(axs, data, task_list):
#     jun_list, dims_ls, temp_list, arm_num_list, sp_suffix = data
#     # 1. a horizontal dashed line at 5/32 for all mean
#     for i in range(len(arm_num_list)):
#         axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
#     return axs

def summ_plot_js_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix=''):
    '''
    Summary plot (#unpaired at junction varied) of junction shift: distance between CoM of the whole nanostar and the CoM of the junction.
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    assert len(color_list) == len(marker_list) == len(temp_list)
    varname = 'js'
    #### plot confs ####
    xlim = (-1, 11)
    ylim_avg = (0, 1.5)
    ylim_std = (0, 1.5)
    ylim_skw = (0.4, 2.0)
    y_var = 'Junction Shift (nm)'
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    #### conf ends ####
    # packing
    data = (jun_list, dims_ls, temp_list, arm_num_list, sp_suffix)
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    # load
    jun_summ_dic, savepath = SL_jun(ns_js_plot, data, conc_list, varname)
    #### Plot Summaries ####
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # data: conc ~ temperature ~ (jun~summ)
    for conc in conc_list:
        plt = summ_plot_jun(jun_summ_dic, plot_confs, data, conc, task_list, color_list, marker_list)
        chkdir(os.path.dirname(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png'))
        plt.savefig(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png',dpi=500)
        plt.clf()
    return True