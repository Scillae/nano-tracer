from utils.tools import save_load
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collections import OrderedDict
from joblib import Parallel, delayed

def parallel_ns_func_decorator(ns_func,arm_num,conc,temp,data):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix, flag_suffix = data
    cond = (arm_num,conc,temp)
    m1, std, m3_s = ns_func(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls=dims_ls, sp_suffix=sp_suffix, flag_suffix=flag_suffix)
    return cond, (m1, std, m3_s)

def SL(ns_func, data, varname):
    '''
    Load the result of summary.
    If not exists, calculate using ns_func, and then save.
    :ns_func: a ns_func that calculate and summarize the distribution of the desired measurement.
    :data: in which the descriptions of nanostars (trajectory) are stored.
    :varname: codename of the measurement.
    '''
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix, flag_suffix = data
    savepath = f'summary/{arm_num_list}Arms{conf_suffix}{flag_suffix}{sp_suffix}/{temp_list}C-{conc_list}M'
    su_path = f'{savepath}-{varname}.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
    su_dic_results = save_load(su_path, None)
    if su_dic_results == False:
        summary_dic = OrderedDict() # {(keys):(mn)}
        # for arm_num in arm_num_list:
        #     for conc in conc_list:
        #         for temp in temp_list:
        #             par_cond_list.append((arm,temp,conc))
        #             m1, std, m3_s = ns_func(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls=dims_ls, sp_suffix=sp_suffix, flag_suffix=flag_suffix)
        #             summary_dic[(arm_num,conc,temp)] = (m1,std,m3_s)
        par_cond_list = [(arm_num,conc,temp) for arm_num in arm_num_list for conc in conc_list for temp in temp_list]
        r_ls = Parallel(n_jobs=12)(delayed(parallel_ns_func_decorator)(ns_func, arm, conc, temp, data) for arm, conc, temp in par_cond_list)
        for cond, dic in r_ls:
            summary_dic[cond] = dic
        su_dic_results = save_load(su_path, summary_dic)      
    return su_dic_results, savepath

def SL_jun(ns_func, data, conc_list, varname):
    '''
    Load the result of summary. Design of nanostars varied (#unpaired junction bases).
    If not exists, calculate using ns_func, and then save.
    :ns_func: a calc_func that calculate the desired measurement.
    :conc_list: desired concentrations. Each conc one plot, so provided separately.
    :data: in which the descriptions of nanostars (trajectory) are stored.
    :varname: codename of the measurement.
    '''
    jun_list, dims_ls, temp_list, arm_num_list, sp_suffix, flag_suffix = data
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # assume saved, read corr. dics
    jun_summ_dic = OrderedDict() # {jun:{(keys):(mn)}}
    assert len(jun_list) > 0
    for jun in jun_list:
        if jun == 2:
            conf_suffix = ''
        else:
            conf_suffix = f'-jun_{jun}'
            dims_ls[1] = jun
        savepath = f'summary/{arm_num_list}Arms{conf_suffix}{flag_suffix}{sp_suffix}/{temp_list}C-{conc_list}M'
        su_path = f'{savepath}-{varname}.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
        su_dic_results = save_load(su_path, None)
        if su_dic_results == False:
            summary_dic = OrderedDict() # {(keys):(mn)}
            for arm_num in arm_num_list:
                for conc in conc_list:
                    for temp in temp_list:
                        m1, std, m3_s = ns_func(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, flag_suffix=flag_suffix, dims_ls=dims_ls)
                        summary_dic[(arm_num,conc,temp)] = (m1,std,m3_s)
            su_dic_results = save_load(su_path, summary_dic)      
        jun_summ_dic[jun] = su_dic_results
    savepath = f'summary/{arm_num_list}Arms{jun_list}/{temp_list}C-{conc_list}M'
    return jun_summ_dic, savepath

def summ_plot(summary_dic, plot_confs, data, task_list, color_list, marker_list, special_tasks=None, sp_suffix='', flag_suffix=''):
    '''
    Plot the summary of a measurement.
    Plot: len(task_list) * len(arm_num_list), var vs temp, series: conc
    :summary_dic: summary of measurement.
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :task_list: tasks to be plotted.
    ...
    :special_tasks: a function that modifies the plots.
    '''
    xlim, ylim_avg, ylim_std, ylim_skw, y_var = plot_confs # unpack configurations
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, sp_suffix, flag_suffix = data
    # plot
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
                    axs[i, j].set_ylabel(rf'{task} of {y_var}')
                else: # sync x,y axes across a row
                    axs[i, j].sharex(axs[i, 0])
                    axs[i, j].sharey(axs[i, 0])
    axs[0,len(arm_num_list)-1].legend()
    if special_tasks is not None:
        axs = special_tasks(axs, data, task_list)
    return plt

def summ_plot_jun(jun_summ_dic, plot_confs, data, conc, task_list, color_list, marker_list, special_tasks=None, sp_suffix=''):
    '''
    Plot the summary of a measurement. Design of nanostars varied (#unpaired junction bases).
    Plot: len(task_list) * len(arm_num_list), var vs #unpaired bases, series: temp
    :summary_dic: summary of measurement.
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :task_list: tasks to be plotted.
    :conc: concentration. Note that different conc ~ different summary.
    ...
    :special_tasks: a function that modifies the plots.
    '''
    xlim, ylim_avg, ylim_std, ylim_skw, y_var = plot_confs # unpack configurations
    jun_list, dims_ls, temp_list, arm_num_list, sp_suffix, flag_suffix = data
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
                axs[i,j].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=2)) # TODO
                if j == 0: # first subplot of a row
                    # set limits of axes.
                    axs[i, j].set_xlim(xlim)
                    if i == 0:
                        axs[i, j].set_ylim(ylim_avg)
                    elif i == 1:
                        axs[i, j].set_ylim(ylim_std)
                    elif i == 2:
                        axs[i, j].set_ylim(ylim_skw)
                    axs[i, j].set_xlabel(r'Number of junction bases')
                    axs[i, j].set_ylabel(rf'{task} of {y_var}')
                else: # sync x,y axes across a row
                    axs[i, j].sharex(axs[i, 0])
                    axs[i, j].sharey(axs[i, 0])
    axs[0,len(arm_num_list)-1].legend()
    plt.suptitle(f'{conc}M group')
    if special_tasks is not None:
        axs = special_tasks(axs, data, task_list)
    return plt