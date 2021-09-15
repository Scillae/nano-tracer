from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from plot_tasks.ns_plots import ns_pa_plot as ns_pa_plot
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os.path


def summ_plot_pa(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    assert len(conc_list) == len(color_list) == len(marker_list)
    #### SL ####
    savepath = f'summary/{arm_num_list}Arms{conf_suffix}/{temp_list}C-{conc_list}M'
    su_path = f'{savepath}-pa.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
    su_dic_results = save_load(su_path, None)
    if su_dic_results == False:
        summary_dic = {} # {'summ':{(keys):(mn)}}
        for arm_num in arm_num_list:
            for conc in conc_list:
                for temp in temp_list:
                    m1, std, m = ns_pa_plot(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls= dims_ls)
                    summary_dic[(arm_num,conc,temp)] = (m1,std,m)
        su_dic_results = save_load(su_path, summary_dic)      
    summary_dic = su_dic_results
    #### SL ends ####
    #### Plot Summaries ####
    fig = plt.figure(figsize=(3*len(arm_num_list),3*len(task_list)))
    gs = fig.add_gridspec(len(task_list), len(arm_num_list), hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
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
                axs[i,j].scatter(temp_list, m_ls,c=color,marker=marker, label=f'{conc} M NaCl') # s=128
                if i == 0:
                    axs[i,j].set_title(rf'{arm_num} arms') # , fontsize=16
                axs[i,j].tick_params(bottom=True, top=True, left=True, right=True)
                axs[i,j].tick_params(axis="x", direction="in")
                axs[i,j].tick_params(axis="y", direction="in")
                axs[i,j].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                if j == 0: # first subplot of a row
                    # set limits of axes.
                    axs[i, j].set_xlim((15,55))
                    if i == 0:
                        axs[i, j].set_ylim((60, 130))
                    elif i == 1:
                        axs[i, j].set_ylim((0, 70))
                    elif i == 2:
                        axs[i, j].set_ylim((-1.3, 1.3))
                    axs[i, j].set_xlabel(r'Temperature ($^\circ$C)')
                    axs[i, j].set_ylabel(rf'{task} of Patch Angles ($^\circ$)')
                else: # sync x,y axes across a row
                    axs[i, j].sharex(axs[i, 0])
                    axs[i, j].sharey(axs[i, 0])
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
    chkdir(os.path.dirname(f'{savepath}-pa.png'))
    plt.savefig(f'{savepath}-pa.png',dpi=500)
    # plt.show()
    plt.clf()
    return True
