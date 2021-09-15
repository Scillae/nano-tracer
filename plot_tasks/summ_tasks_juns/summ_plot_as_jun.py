from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from plot_tasks.ns_plots import ns_as_plot as ns_as_plot
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os.path


def summ_plot_as_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    assert len(color_list) == len(marker_list) == len(temp_list)
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # assume saved, read corr. dics
    jun_summ_dic = {} # {jun:{(keys):(mn)}}
    for jun in jun_list:
        if jun == 2:
            conf_suffix = ''
        else:
            conf_suffix = f'-jun_{jun}'
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
        jun_summ_dic[jun] = su_dic_results
    #### Plot Summaries ####
    # reorganize: conc ~ temperature ~ (jun~summ)
    for conc in conc_list:
        # a figure
        fig = plt.figure(figsize=(3*len(arm_num_list),3*len(task_list)))
        gs = fig.add_gridspec(len(task_list), len(arm_num_list), hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
        axs = gs.subplots(sharey='row')
        for i,task in enumerate(task_list):
            # a row
            for j,arm_num in enumerate(arm_num_list):
                # a cols
                for temp,color,marker in zip(temp_list, color_list, marker_list):
                    # a series
                    m_ls = []
                    for jun in jun_list:
                        # a point
                        m_ls.append(jun_summ_dic[jun][(arm_num,conc,temp)][i])
                    # draw
                    axs[i,j].scatter(jun_list, m_ls,c=color,marker=marker, label=rf'{temp} ($^\circ$C)') # s=128
                    if i == 0:
                        axs[i,j].set_title(rf'{arm_num} arms') # , fontsize=16
                    axs[i,j].tick_params(bottom=True, top=True, left=True, right=True)
                    axs[i,j].tick_params(axis="x", direction="in")
                    axs[i,j].tick_params(axis="y", direction="in")
                    axs[i,j].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                    axs[i,j].xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
                    if j == 0: # first subplot of a row
                        # set limits of axes.
                        axs[i, j].set_xlim((-1,11))
                        if i == 0:
                            axs[i, j].set_ylim((0, 0.3))
                        elif i == 1:
                            axs[i, j].set_ylim((0, 0.3))
                        elif i == 2:
                            axs[i, j].set_ylim((0.4, 2.0))
                        axs[i, j].set_xlabel(r'Number of junction bases')
                        axs[i, j].set_ylabel(rf'{task} of Arm Stiffnesses')
                    else: # sync x,y axes across a row
                        axs[i, j].sharex(axs[i, 0])
                        axs[i, j].sharey(axs[i, 0])
        axs[0,len(arm_num_list)-1].legend()
        plt.suptitle(f'{conc}M group')
        # special tasks
        # 1. a horizontal dashed line at 5/32 for all mean
        for i in range(len(arm_num_list)):
            axs[0,i].plot((-1,11), (5/32,5/32),c='#1AA555',ls=':')
        # special tasks ends        
        chkdir(os.path.dirname(f'{savepath}-as-jun{jun_list}-{conc}M.png'))
        plt.savefig(f'{savepath}-as-jun{jun_list}-{conc}M.png',dpi=500)
        plt.clf()
    return True
