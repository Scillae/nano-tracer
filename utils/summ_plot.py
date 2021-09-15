from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from plot_tasks.ns_plots import ns_as_plot as ns_as_plot
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os.path

def summ_plot(summary_dic, plot_confs, conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list):
    xlim, ylim_avg, ylim_std, ylim_skw = plot_confs # unpack configurations
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
    return plt