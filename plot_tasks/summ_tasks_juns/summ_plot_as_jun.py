from utils.tools import chkdir
from utils.summ_plot import summ_plot_jun, SL_jun
from plot_tasks.ns_plots.ns_as_plot import ns_as_plot
import os.path


def summ_plot_as_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    assert len(color_list) == len(marker_list) == len(temp_list)
    varname = 'as'
    #### plot confs ####
    xlim = (-1, 11)
    ylim_avg = (0., 0.5)
    ylim_std = (0, 0.3)
    ylim_skw = (1, 2.5)
    y_var = 'Arm Stiffnesses'
    #### conf ends ####
    # packing
    data = (jun_list, dims_ls, temp_list, arm_num_list)
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    # load
    jun_summ_dic, savepath = SL_jun(ns_as_plot, data, conc_list, varname)
    #### Plot Summaries ####
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # data: conc ~ temperature ~ (jun~summ)
    for conc in conc_list:
        plt = summ_plot_jun(jun_summ_dic, plot_confs, data, conc, task_list, color_list, marker_list)
        chkdir(os.path.dirname(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png'))
        plt.savefig(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png',dpi=500)
        plt.clf()
    return True
