from utils.tools import chkdir
from utils.summ_plot import summ_plot, SL
from plot_tasks.ns_plots.ns_pa_plot import ns_pa_plot
import os.path

def special_tasks(axs, data, task_list):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix = data
    # special tasks
    # 1. a horizontal dashed line at 109.5° for the 4 arm
    axs[0,1].plot((15,55), (109.5,109.5),c='#1AA555',ls=':')
    # 2. one at 90° for the 6 arm patch angle plot. 
    axs[0,3].plot((15,55), (90,90),c='#1AA555',ls=':')
    # 3. one at 120° for the 3 arm patch angle plot
    axs[0,0].plot((15,55), (120,120),c='#1AA555',ls=':')
    # 4. a solid horizontal black line at skew=0° in the bottom row of plots.
    for i in range(len(arm_num_list)):
        axs[2,i].plot((15,55), (0,0),c='#000000')
    # special tasks ends
    axs[0,len(arm_num_list)-1].legend()
    return axs


def summ_plot_pa(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list, sp_suffix=''):
    '''
    Summary plot (design fixed) of patch angle: angle of two arms that share a common strand.
    Set varname and plot confs.
    Define special tasks to customize the plot.
    '''
    # assert len(conc_list) == len(color_list) == len(marker_list)
    varname = 'pa'
    #### plot confs ####    
    xlim = (15,55)
    ylim_avg = (60, 150)
    ylim_std = (0, 90)
    ylim_skw = (-1.3, 1.3)
    y_var = rf'Patch Angles ($^\circ$)'
    #### conf ends ####
    plot_confs = (xlim, ylim_avg, ylim_std, ylim_skw, y_var)
    data = conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix
    # load data
    summary_dic, savepath = SL(ns_pa_plot, data, varname)
    # plot
    plt = summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, special_tasks, sp_suffix=sp_suffix)
    chkdir(os.path.dirname(f'{savepath}-{varname}.png'))
    plt.savefig(f'{savepath}-{varname}.png',dpi=500)
    plt.clf()
    return True
