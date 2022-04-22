from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_k2_plot import ns_k2_plot
import os.path


def special_tasks(axs, data, task_list):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix = data
    # 1. a horizontal dashed line at 5/32 for all mean
    for i in range(len(arm_num_list)):
        axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
    return axs


def summ_plot_k2(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix=''):
    '''
    Summary plot (design fixed) of k2 ==1-27*λ1*λ2*λ3/(λ1+λ2+λ3)^3, where λi is the principal value of the tensor of inertia along an axis.
    k2 ---> 1: anisotropic; ---> 0: isotropic.
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    # assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'k2'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (0, 0.3)
    ylim_std = (0, 0.3)
    ylim_skw = (0.4, 2.0)
    y_var = 'k2 values'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix
    # load data
    summary_dic, savepath = SL(ns_k2_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, special_tasks, sp_suffix=sp_suffix)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True
