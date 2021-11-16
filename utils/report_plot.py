from utils.tools import save_load
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from scipy import optimize

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

def data_process_func_js(js_ls_res, data):
    # arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    js_ls = [i[1]*10*0.8518 for i in js_ls_res]
    return js_ls

def SL(ns_func, data, varname):
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    savepath = f'summary/{arm_num_list}Arms{conf_suffix}/{temp_list}C-{conc_list}M'
    su_path = f'{savepath}-{varname}.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
    su_dic_results = save_load(su_path, None)
    if su_dic_results == False:
        summary_dic = {} # {(keys):(mn)}
        for arm_num in arm_num_list:
            for conc in conc_list:
                for temp in temp_list:
                    m1, std, m3_s = ns_func(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls=dims_ls)
                    summary_dic[(arm_num,conc,temp)] = (m1,std,m3_s)
        su_dic_results = save_load(su_path, summary_dic)      
    return su_dic_results, savepath

def SL_jun(ns_func, data, conc_list, varname):
    jun_list, dims_ls, temp_list, arm_num_list = data
    # plot: conc ~ {x: jun_nums, y: summaries, series: temperature}
    # assume saved, read corr. dics
    jun_summ_dic = {} # {jun:{(keys):(mn)}}
    assert len(jun_list) > 0
    for jun in jun_list:
        if jun == 2:
            conf_suffix = ''
        else:
            conf_suffix = f'-jun_{jun}'
            dims_ls[1] = jun
        savepath = f'summary/{arm_num_list}Arms{conf_suffix}/{temp_list}C-{conc_list}M'
        su_path = f'{savepath}-{varname}.sudic' # (p_angs_dic, ns_tm, ns_last, sys)
        su_dic_results = save_load(su_path, None)
        if su_dic_results == False:
            summary_dic = {} # {(keys):(mn)}
            for arm_num in arm_num_list:
                for conc in conc_list:
                    for temp in temp_list:
                        m1, std, m3_s = ns_func(single=True, arms=arm_num, temp=temp, conc=conc, conf_suffix=conf_suffix, dims_ls=dims_ls)
                        summary_dic[(arm_num,conc,temp)] = (m1,std,m3_s)
            su_dic_results = save_load(su_path, summary_dic)      
        jun_summ_dic[jun] = su_dic_results
    savepath = f'summary/{arm_num_list}Arms{jun_list}/{temp_list}C-{conc_list}M'
    return jun_summ_dic, savepath

def pa_3d_report_plot(data):
    z_arm_id = 0
    x_lim = (-10,10)
    y_lim = (-10,10)
    z_lim = (-10,10)
    from utils.ns_plot import SL_ns
    from calc_tasks.k2_calc import k2_calc
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    color_list = ['#4994FF','#E55050','#FCC555','#7AA77A']
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d', xlim=x_lim, ylim=y_lim, zlim=z_lim)
    for j, conc in enumerate(conc_list):
        for i, temp in enumerate(temp_list):
            k2_arr_ls = []
            for arm_num in arm_num_list:
                # savepath = f'data/composed_traj/{arm_num}arms{conf_suffix}/{loose_lbl}/{label}'
                data = (arm_num, temp, conc, '', conf_suffix, dims_ls)
                ns_tm, label, plotpath = SL_ns(None, data, 'ns')
                vec_dic = {} # {t:[(ia,(x,y,z),(shareia,shareia2))]}
                for t_stamp, ns in ns_tm.time_capsule.items():
                    # print(t_stamp)
                    vec_dic[t_stamp] = nanostar_vectorize(ns,z_arm_id)
                # prepare for plotting
                arms_seq_dic = {} # {ia:{seq_i}}
                seq_x = []
                seq_y = []
                seq_z = []
                for t_stamp, vec_arms_ls in vec_dic.items():
                    for (ia,(x,y,z),(shareia,shareia2), is_sharing_sel) in vec_arms_ls:
                        if ia not in arms_seq_dic.keys():
                            arms_seq_dic[ia] = {}
                            arms_seq_dic[ia]['x'] = []
                            arms_seq_dic[ia]['y'] = []
                            arms_seq_dic[ia]['z'] = []
                        arms_seq_dic[ia]['x'].append(x)
                        arms_seq_dic[ia]['y'].append(y)
                        arms_seq_dic[ia]['z'].append(z)
                        arms_seq_dic[ia]['is_sharing_sel'] = is_sharing_sel
                    #     ax.scatter(x,y,z)
                    #     print(f"Arm:{ia}, coord:({x:.5f},{y:.5f},{z:.5f})")
                    # plt.show()
                # additional data processing goes here
                for ia, seq_dic in arms_seq_dic.items():
                    if ia == z_arm_id:
                        ax.scatter(seq_dic['x'],seq_dic['y'],seq_dic['z'],c=color_list[0], s=5)
                    else:
                        if seq_dic['is_sharing_sel'] == 1:
                            ax.scatter(seq_dic['x'],seq_dic['y'],seq_dic['z'],c=color_list[1], s=5)
                            assert 1==1
                        elif seq_dic['is_sharing_sel'] == 2:
                            ax.scatter(seq_dic['x'],seq_dic['y'],seq_dic['z'],c=color_list[2], s=5)
                            assert 1==1
                        elif seq_dic['is_sharing_sel'] == False:
                            ax.scatter(seq_dic['x'],seq_dic['y'],seq_dic['z'],c=color_list[3], s=5)
                            assert 1==1
                        else:
                            assert 0==1
    plt.show()
    return plt

def junction_CoM(ns):
    # computing CoM of Junction
    CoM_pos = np.zeros(3)
    base_cnt = 0
    center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
    center_ls.extend([base for arm in ns.arms.values() for base in list(arm.base_pairs.values())[-1]])
    for base in center_ls:
        CoM_pos = np.add(CoM_pos, np.array(base.position))
        base_cnt += 1
    CoM_pos = np.divide(CoM_pos, base_cnt)
    return CoM_pos

def coord_rotate(ns, z_arm_id):
    z_pair = [base for base in list(ns.arms[z_arm_id].base_pairs.values())[-1]] # start
    z_end_pair = [base for base in list(ns.arms[z_arm_id].base_pairs.values())[0]] # end
    # y = np.array(z_pair[0].position) - np.array(z_pair[1].position)
    # y_unit = y/np.linalg.norm(y)
    # y normal 2
    z_arm = ns.arms[z_arm_id]
    s_arm = False
    for arm in ns.arms.values(): # this iteration always goes in the same order.
        if arm.strand_id_0 == z_arm.strand_id_0 or arm.strand_id_0 == z_arm.strand_id_1 or arm.strand_id_1 == z_arm.strand_id_0 or arm.strand_id_1 == z_arm.strand_id_1:
            s_arm = arm
    s_pair = [base for base in list(s_arm.base_pairs.values())[-1]]
    y_unit = np.cross(np.array(z_pair[0].normal), np.array(s_pair[0].normal))
    z = (np.array(z_end_pair[1].position) + np.array(z_end_pair[0].position))/2 - (np.array(z_pair[1].position) + np.array(z_pair[0].position))/2
    z_unit = z/np.linalg.norm(z)
    x_unit = np.cross(y_unit,z_unit)
    rot_mat = np.vstack((x_unit,y_unit,z_unit)).T
    rot_mat = np.linalg.inv(rot_mat)
    return rot_mat

def nanostar_vectorize(ns, sel_arm_id):
    vec_ns_ls = []
    jun_CoM = junction_CoM(ns)
    rot_mat = coord_rotate(ns, sel_arm_id)
    for ia, arm in ns.arms.items():
        end_base_pair = [base for base in list(arm.base_pairs.values())[0]] # -1 is center
        start_base_pair = [base for base in list(arm.base_pairs.values())[-1]] # -1 is center
        # pair_pos = tuple(sum(x)/2 for x in zip(end_base_pair[0].position, end_base_pair[1].position))
        end_pos = (np.array(end_base_pair[1].position) + np.array(end_base_pair[0].position))/2
        start_pos = (np.array(start_base_pair[1].position) + np.array(start_base_pair[0].position))/2
        vec = end_pos - start_pos
        vec = np.matmul(rot_mat,vec)
        strands_using = (arm.strand_id_0, arm.strand_id_1)
        is_sharing_sel = ns.arms[sel_arm_id].strand_id_0 in strands_using or ns.arms[sel_arm_id].strand_id_1 in strands_using
        if is_sharing_sel:
            is_sharing_sel = 1 if ns.arms[sel_arm_id].strand_id_0 in strands_using else 2
        # print(end_base_pair[0].base_id, end_base_pair[1].base_id)
        vec_ns_ls.append((ia, tuple(vec), strands_using, is_sharing_sel))
    return vec_ns_ls

def k2_report_plot(summary_dic, plot_confs, data, color_list, marker_list, special_tasks=None):
    from utils.ns_plot import SL_ns
    from calc_tasks.k2_calc import k2_calc
    # mean only
    msize_list = [9,6]
    ebsize_list = [8,5]
    xlim, ylim_avg, ylim_std, ylim_skw, y_var = plot_confs # unpack configurations
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    # plot
    fig = plt.figure(figsize=(3*len(conc_list), 3*1+0.5)) # figsize=(3*len(conc_list), 3*1)
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Arm Number', labelpad= 0)
    plt.ylabel(r'$k_2$', rotation=0, fontsize=15, usetex = True, labelpad= 10)
    gs = fig.add_gridspec(1, len(conc_list), hspace=0, wspace=0) # 
    axs = gs.subplots(sharey='row')
    for j, conc in enumerate(conc_list):
        for i, temp in enumerate(temp_list):
            k2_arr_ls = []
            for arm_num in arm_num_list:
                # m_ls.append(summary_dic[(arm_num,conc,temp)][0]) # m
                # m_std.append(summary_dic[(arm_num,conc,temp)][1])   
                data = (arm_num, temp, conc, '', conf_suffix, dims_ls)
                k2_ls, label, plotpath = SL_ns(k2_calc, data, 'k2')
                k2_arr = np.array([tu[1] for tu in k2_ls])
                k2_arr_ls.append(k2_arr)
            # axs[j].violinplot(k2_arr_ls, positions=[a-0.3/2+0.3*i for a in arm_num_list], split=True)

            # boxplot
            bp = axs[j].boxplot(k2_arr_ls, positions=[a-0.3/2+0.3*i for a in arm_num_list], widths = 0.2)
            for whisker in bp['whiskers']:
                whisker.set(color =color_list[i], linewidth = 1)        
            for cap in bp['caps']:
                cap.set(color =color_list[i], xdata=cap.get_xdata() + (-0.1,+0.1), linewidth = 1)            
            for box in bp['boxes']:
                box.set(color =color_list[i])
            for median in bp['medians']:
                median.set(color =color_list[i])
            for flier in bp['fliers']:
                flier.set(markersize=3, markeredgecolor='None')
                flier.set_markerfacecolor(color_list[i])
                flier.set_marker('.')

            # # errorbar
            # axs[j].errorbar(arm_num_list, m_ls, yerr=m_std, c=color_list[i],markersize=msize_list[i], marker=marker_list[i], label=f'{temp} ℃', ls='none', capsize=ebsize_list[i])
            axs[j].tick_params(axis="y", left=True, right=True, direction='in',which='both')
        # print(f'Plot[{j}], Arm_num_ls:{arm_num_list}, Mean_ls:{m_ls}, STD_ls:{m_std}. T={temp}C,[NaCl]={conc}M')
        axs[j].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=0.05))
        axs[j].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.1))
        axs[j].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))
        axs[j].set_xticklabels((0,3,4,5,6))
        axs[j].tick_params(axis='x', length=0)
        axs[j].set_xlim(xlim)
        axs[j].set_ylim(ylim_avg)
        axs[j].set_title(f'{conc} M [NaCl]')
        for arm_num in arm_num_list[:-1]:
            axs[j].axvline(arm_num+0.5,c='black', linewidth=0.5, linestyle='--')
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=color_list[i], lw=2) for i in range(len(temp_list))]
    legend = axs[0].legend(custom_lines, [f'{temp}℃' for temp in temp_list], loc='upper left')
    legend.get_frame().set_alpha(None)
    legend.get_frame().set_facecolor((1, 1, 1, 1))
    legend.get_frame().set_linewidth(0.0)
    return plt


def js_report_plot():
    color_list = ['#4994FF','#E55050','#FCC555','#7AA77A']
    #### plot confs ####
    x_var = r'Junction Shift $\textup{~\AA}$'
    x_lim = (0,5.5)
    y_lim = (0,0.10)
    bin_num = 50
    text_pos = (0.7*(x_lim[1]-x_lim[0])+x_lim[0], (0.215/0.3)*(y_lim[1]-y_lim[0]))
    from utils.ns_plot import SL_ns
    arm_nums = [6,5,4,3]
    conc = 0.5
    temp = 20
    from calc_tasks.jun_shift_calc import jun_shift_calc
    for fig_idx, arms in enumerate(arm_nums):
        data = (arms, temp, conc, '', '', [20,2,7])

        result_val, _, _ = SL_ns(jun_shift_calc,data,'js')
    
        var_ls = data_process_func_js(result_val, data)
        print(max(var_ls))
        if min(var_ls) < x_lim[0] or max(var_ls) > x_lim[1]:
            print(f'ns_plot out of range: {min(var_ls)} < {x_lim[0]} : {min(var_ls) < x_lim[0]} or {max(var_ls)} > {x_lim} : {max(var_ls) > x_lim[1]}')
            # assert 0 == 1
        n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

        # moments
        from utils.ns_plot import moments_calc
        n, m1, std, m3_s = moments_calc(n, var_ls)

        # ax = plt.axes()
        n_sum = np.sum(n)
        cut_off = 0.015
        if True:
            from copy import deepcopy
            var_trun = deepcopy(var_ls)
            var_tmp = []
            n_max = max(n)
            n_t = []
            for v,bin_lb,bin_ub in zip(n,bin_edges[:-1],bin_edges[1:]):
                if v < cut_off*1000: # n_max*0.1
                    var_tmp = [var for var in var_trun if (var>bin_ub or var<bin_lb) ]
                    var_trun = deepcopy(var_tmp)
                else:
                    n_t.append(v)
            n_t = np.array(n_t)
            n_t, m1_t, std_t, m3_s_t = moments_calc(n_t, var_trun)
            from scipy.stats import norm
            fit_x = np.linspace(x_lim[0], x_lim[1], len(var_trun))
            fit_y = norm.pdf(fit_x,m1_t,std_t)            
        else:
            from scipy.stats import norm
            fit_x = np.linspace(x_lim[0], x_lim[1], len(var_ls))
            fit_y = norm.pdf(fit_x,m1,std)
        '''
        x = np.linspace(x_lim[0], x_lim[1], len(var_ls))
        popt, _ = optimize.curve_fit(gaussian, x, var_ls,maxfev=5000)cao
        y = gaussian(x, *popt)
        '''
        plt.plot(fit_x, fit_y*(x_lim[1]-x_lim[0])/bin_num, c=color_list[fig_idx]) # , label=rf'$\mu$={m1:.2f}, $\sigma$={std:.2f}'
        print(np.sum(np.divide(np.array(n),n_sum)))
        print(arms, m1, std)
        plt.errorbar(bin_centers,np.divide(np.array(n),n_sum),yerr = np.divide(np.sqrt(np.array(n)),n_sum),marker = '.',drawstyle = 'steps-mid',c=color_list[fig_idx],label=f'{arms} arms', ls='none')
    # plt.title(f'{temp}C,{conc}M', fontsize=18)
    plt.plot([x_lim[0],x_lim[1]], [cut_off,cut_off], c='#000000', ls='--') # : {0.008}
    plt.legend()
    # # customization
    # txt = (
    #     rf'$\mu={m1:.2f}$'
    #     '\n'
    #     rf'$\sigma={std:.2f}$'
    #     '\n'
    #     rf'$\gamma={m3_s:.2f}$'
    # )
    # text_x, text_y = text_pos
    # plt.text(text_x, text_y, txt, c='#4994FF', fontsize=18) # 28%
    plt.xlabel(x_var, fontsize=18, usetex=True)
    plt.ylabel('Frequency', fontsize=18,labelpad=5)
    plt.xlim(x_lim) # x_lim
    plt.ylim(y_lim)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    from utils.tools import chkdir
    from os.path import dirname
    plotpath = 'report/js.png'
    chkdir(dirname(plotpath))
    plt.savefig(plotpath,dpi=800)
    plt.show()
    plt.clf()
    return True


def k2_report_plot_old(summary_dic, plot_confs, data, color_list, marker_list, special_tasks=None):
    # mean only
    msize_list = [9,6]
    ebsize_list = [8,5]
    xlim, ylim_avg, ylim_std, ylim_skw, y_var = plot_confs # unpack configurations
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    # plot
    fig = plt.figure(figsize=(3*len(conc_list), 3*1+0.5)) # figsize=(3*len(conc_list), 3*1)
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Arm Number')
    plt.ylabel(r'$k_2$', rotation=0, fontsize=15, usetex = True)
    gs = fig.add_gridspec(1, len(conc_list), hspace=0, wspace=0) # 
    axs = gs.subplots(sharey='row')
    task = 'Mean'
    for j, conc in enumerate(conc_list):
        for i, temp in enumerate(temp_list):
            m_ls = []
            m_std = []
            for arm_num in arm_num_list:
                m_ls.append(summary_dic[(arm_num,conc,temp)][0]) # m
                m_std.append(summary_dic[(arm_num,conc,temp)][1])   
            axs[j].errorbar(arm_num_list, m_ls, yerr=m_std, c=color_list[i],markersize=msize_list[i], marker=marker_list[i], label=f'{temp} ℃', ls='none', capsize=ebsize_list[i])
            axs[j].tick_params(axis="y", left=True, right=True, direction='in',which='both')
        # print(f'Plot[{j}], Arm_num_ls:{arm_num_list}, Mean_ls:{m_ls}, STD_ls:{m_std}. T={temp}C,[NaCl]={conc}M')
        axs[j].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=0.05))
        axs[j].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.1))
        axs[j].set_xlim(xlim)
        axs[j].set_ylim(ylim_avg)
        axs[j].set_title(f'{conc} M [NaCl]')
    handles, labels = axs[0].get_legend_handles_labels()
    handles = [h[0] for h in handles]
    axs[0].legend(reversed(handles), reversed(labels), loc='upper left')
    # axs[0].set_ylabel(r'$k_2$', rotation=0, fontsize=15,fontname='serif',usetex = True,position=(100,0))
    # for i, conc in enumerate(conc_list):
    #     axs[i].xaxis.set_visible(True)
    #     axs[i].set_xlabel('Arm Number',fontname='serif')
    '''
    for j,arm_num in enumerate(arm_num_list):
        # separate plots
        for conc, color, marker in zip(conc_list, color_list, marker_list):
            # separate series
            m_ls = []
            m_std = []
            for temp in temp_list:
                m_ls.append(summary_dic[(arm_num,conc,temp)][0]) # m
                m_std.append(summary_dic[(arm_num,conc,temp)][1])
            # draw
            # axs[j].scatter(temp_list, m_ls,c=color,marker=marker, label=f'{conc} M NaCl')
            axs[j].errorbar(temp_list, m_ls, yerr=m_std, c=color, marker=marker, label=f'{conc} M NaCl')
            axs[j].set_title(rf'{arm_num} arms') # , fontsize=16
            axs[j].tick_params(bottom=True, top=True, left=True, right=True)
            axs[j].tick_params(axis="x", direction="in")
            axs[j].tick_params(axis="y", direction="in")
            axs[j].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            if j == 0: # first subplot of a row
                # set limits of axes.
                axs[j].set_xlim(xlim)
                axs[j].set_ylim(ylim_avg)
                axs[j].set_xlabel(r'Temperature ($^\circ$C)')
                axs[j].set_ylabel(rf'{task} of {y_var}')
            else: # sync x,y axes across a row
                axs[j].sharex(axs[0])
                axs[j].sharey(axs[0])
    axs[len(arm_num_list)-1].legend()
    if special_tasks is not None:
        axs = special_tasks(axs, data, None)
    '''
    return plt


def k2_report_plot_seaborn(summary_dic, plot_confs, data, color_list, marker_list, special_tasks=None):
    from utils.ns_plot import SL_ns
    from calc_tasks.k2_calc import k2_calc
    import seaborn as sns
    import pandas as pd
    from copy import deepcopy
    # mean only
    msize_list = [9,6]
    ebsize_list = [8,5]
    xlim, ylim_avg, ylim_std, ylim_skw, y_var = plot_confs # unpack configurations
    conf_suffix, dims_ls, conc_list, temp_list, arm_num_list = data
    # plot
    fig = plt.figure(figsize=(3*len(conc_list), 3*1+0.5)) # figsize=(3*len(conc_list), 3*1)
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Arm Number', labelpad= 0)
    plt.ylabel(r'$k_2$', rotation=0, fontsize=15, usetex = True, labelpad= 10)
    gs = fig.add_gridspec(1, len(conc_list), hspace=0, wspace=0) # 
    axs = gs.subplots(sharey='row')
    for j, conc in enumerate(conc_list):
        for i, temp in enumerate(temp_list):
            df_ls = []
            for arm_num in arm_num_list:
                df = pd.DataFrame()
                data = (arm_num, temp, conc, '', conf_suffix, dims_ls)
                k2_ls, label, plotpath = SL_ns(k2_calc, data, 'k2')
                k2_arr = np.array([tu[1] for tu in k2_ls])
                df[arm_num] = k2_arr
                df['Temp'] = temp
                df_ls.append(deepcopy(df))
            df = pd.concat([df for df in df_ls])
            sns.violinplot(x='Scenario', y='LMP', hue='Temp', split=True, data=df, ax = axs[j])
            # # boxplot
            # bp = axs[j].boxplot(k2_arr_ls, positions=[a-0.3/2+0.3*i for a in arm_num_list], widths = 0.2)
            # for whisker in bp['whiskers']:
            #     whisker.set(color =color_list[i], linewidth = 1)        
            # for cap in bp['caps']:
            #     cap.set(color =color_list[i], xdata=cap.get_xdata() + (-0.1,+0.1), linewidth = 1)            
            # for box in bp['boxes']:
            #     box.set(color =color_list[i])
            # for median in bp['medians']:
            #     median.set(color =color_list[i])
            # for flier in bp['fliers']:
            #     flier.set(markersize=3, markeredgecolor='None')
            #     flier.set_markerfacecolor(color_list[i])
            #     flier.set_marker('.')

            # # errorbar
            # axs[j].errorbar(arm_num_list, m_ls, yerr=m_std, c=color_list[i],markersize=msize_list[i], marker=marker_list[i], label=f'{temp} ℃', ls='none', capsize=ebsize_list[i])
            axs[j].tick_params(axis="y", left=True, right=True, direction='in',which='both')
        # print(f'Plot[{j}], Arm_num_ls:{arm_num_list}, Mean_ls:{m_ls}, STD_ls:{m_std}. T={temp}C,[NaCl]={conc}M')
        axs[j].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=0.05))
        axs[j].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.1))
        axs[j].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))
        axs[j].set_xticklabels((0,3,4,5,6))
        axs[j].tick_params(axis='x', length=0)
        axs[j].set_xlim(xlim)
        axs[j].set_ylim(ylim_avg)
        axs[j].set_title(f'{conc} M [NaCl]')
        for arm_num in arm_num_list[:-1]:
            axs[j].axvline(arm_num+0.5,c='black', linewidth=0.5, linestyle='--')
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=color_list[i], lw=2) for i in range(len(temp_list))]
    legend = axs[0].legend(custom_lines, [f'{temp}℃' for temp in temp_list], loc='upper left')
    legend.get_frame().set_alpha(None)
    legend.get_frame().set_facecolor((1, 1, 1, 1))
    legend.get_frame().set_linewidth(0.0)
    return plt
