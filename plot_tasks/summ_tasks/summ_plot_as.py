from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_as_plot import ns_as_plot
import os.path


def summ_plot_as(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix='', flag_suffix=''):
    '''
    Summary plot (design fixed) of arm stiffness: breaking the arms into 2 (customizable) vectors, and summarizing their angles (normalized).
    AS == 0: straight. AS == 1: 180°
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'as'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (0,0.4)
    ylim_std = (0,0.4)
    ylim_skw = (-1,2.5)
    y_var = 'Arm Stiffnesses'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix, flag_suffix
    # load data
    summary_dic, savepath = SL(ns_as_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, sp_suffix=sp_suffix, flag_suffix=flag_suffix)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True
