from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from plot_tasks.ns_plots import ns_as_plot as ns_as_plot
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os.path


def summ_plot_as(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    assert len(conc_list) == len(color_list) == len(marker_list)
    #### plot confs ####
    xlim=(15,55)
    ylim_avg=(0.1,0.4)
    ylim_std=(0.1,0.4)
    ylim_skw=(1,2.5)
    #### conf ends ####
    #### SL ####
    savepath = f'summary/{arm_num_list}Arms{conf_suffix}/{temp_list}C-{conc_list}M'
    su_path = f'{savepath}-as.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
    su_dic_results = save_load(su_path, None)
    if su_dic_results == False:
        summary_dic = {} # {(keys):(mn)}
        for arm_num in arm_num_list:
            for conc in conc_list:
                for temp in temp_list:
                    m1, std, m3_s = ns_as_plot(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls=dims_ls)
                    summary_dic[(arm_num,conc,temp)] = (m1,std,m3_s)
        su_dic_results = save_load(su_path, summary_dic)      
    summary_dic = su_dic_results
    #### SL ends ####
    fig = plt.figure(figsize=(3*len(arm_num_list),3*len(task_list)))
    gs = fig.add_gridspec(len(task_list), len(arm_num_list), hspace=0.3, wspace=0.1)
    axs = gs.subplots(sharey='row')
    for i,task in enumerate(task_list):
        for j,arm_num in enumerate(arm_num_list):
            # separate plots
            for conc, color, marker in zip(conc_list, color_list, marker_list):
                # separate series
                m_ls = []
                for temp in temp_list:
                    m_ls.append(summary_dic[(arm_num,conc,temp)][i]) # m
                # draw
                axs[i,j].scatter(temp_list, m_ls,c=color,marker=marker, label=f'{conc} M NaCl')
                if i == 0:
                    axs[i,j].set_title(rf'{arm_num} arms') # , fontsize=16
                axs[i,j].tick_params(bottom=True, top=True, left=True, right=True)
                axs[i,j].tick_params(axis="x", direction="in")
                axs[i,j].tick_params(axis="y", direction="in")
                axs[i,j].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                if j == 0: # first subplot of a row
                    # set limits of axes.
                    axs[i, j].set_xlim(xlim)
                    if i == 0:
                        axs[i, j].set_ylim(ylim_avg)
                    elif i == 1:
                        axs[i, j].set_ylim(ylim_std)
                    elif i == 2:
                        axs[i, j].set_ylim(ylim_skw)
                    axs[i, j].set_xlabel(r'Temperature ($^\circ$C)')
                    axs[i, j].set_ylabel(rf'{task} of Arm Stiffnesses')
                else: # sync x,y axes across a row
                    axs[i, j].sharex(axs[i, 0])
                    axs[i, j].sharey(axs[i, 0])
    axs[0,len(arm_num_list)-1].legend()
    # special tasks
    # special tasks ends
    chkdir(os.path.dirname(f'{savepath}-as.png'))
    plt.savefig(f'{savepath}-as.png',dpi=500)
    plt.clf()
    return True
