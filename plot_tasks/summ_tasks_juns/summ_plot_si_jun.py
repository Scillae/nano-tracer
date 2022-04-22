from utils.tools import chkdir
from utils.summ_plot import summ_plot_jun, SL_jun
from plot_tasks.ns_plots.ns_si_plot import ns_si_plot
import os.path


def summ_plot_pa_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix=''):
    '''
    Summary plot (#unpaired at junction varied) of ...
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    # assert len(color_list) == len(marker_list) == len(temp_list)
    varname = 'si'
    #### plot confs ####
    xlim = (-1, 11)
    ylim_avg = (60, 150)
    ylim_std = (0, 90)
    ylim_skw = (-1.3, 1.3)
    y_var = 'Unavailable'
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    #### conf ends ####
    # packing
    data = (jun_list, dims_ls, temp_list, arm_num_list, sp_suffix)
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    # load
    jun_summ_dic, savepath = SL_jun(ns_si_plot, data, conc_list, varname)
    #### Plot Summaries ####
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # data: conc ~ temperature ~ (jun~summ)
    for conc in conc_list:
        plt = summ_plot_jun(jun_summ_dic, plot_confs, data, conc, task_list, color_list, marker_list)
        chkdir(os.path.dirname(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png'))
        plt.savefig(f'{savepath}-{varname}-jun{jun_list}-{conc}M.png',dpi=500)
        plt.clf()
    return True