from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from plot_tasks.ns_plots import ns_as_plot as ns_as_plot
from utils.tools import save_load, chkdir
from utils.summ_plot import summ_plot
from collections import OrderedDict
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
    params = (xlim, ylim_avg, ylim_std, ylim_skw)
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
    plt = summ_plot(summary_dic, params, conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # special tasks
    # special tasks ends
    chkdir(os.path.dirname(f'{savepath}-as.png'))
    plt.savefig(f'{savepath}-as.png',dpi=500)
    plt.clf()
    return True
