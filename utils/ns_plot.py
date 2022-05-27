from utils.tools import save_load, chkdir
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats
import os.path
from copy import deepcopy
import pickle

def moments_calc(n, var_ls):
    '''
    Calculate the 0th raw, 1st raw, 2nd central, and 3rd standardized moment of a given distribution.
    '''
    n = np.array(n)
    m0 = np.sum(n) # 0th unitless raw moment: integration
    m1 = np.sum(var_ls)/m0
    m2_c = stats.moment(var_ls, moment=2) # 2nd central moment: variance
    std = m2_c**0.5 # standard deviation
    m3_c = stats.moment(var_ls, moment=3) # 3rd central moment
    m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
    return n, m1, std, m3_s


def SL_ns(calc_func, data, varname, vtime=False):
    '''
    Load the calculated result of a measurement.
    If not exists, calculate using ns_func, and then save.
    :calc_func: a calc_func that calculate the desired measurement.
    :data: in which the descriptions of nanostars (trajectory) are stored.
    :varname: codename of the measurement.
    '''
    arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls = data
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{flag_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}{flag_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
    traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}{flag_suffix}/{loose_lbl}/trajectory.dat'
    savepath = f'data/composed_traj/{arms}arms{conf_suffix}{flag_suffix}/{loose_lbl}/{label}'
    plotpath = f'results/{arms}arms{conf_suffix}{flag_suffix}/{varname}/{varname}_hist-{label}.png'
    if vtime == True:
        plotpath = f'results/{arms}arms{conf_suffix}{flag_suffix}/{varname}/{varname}_vtime-{label}.png'
    if varname == 'sys' or varname == 'ns':
        var_path = f'{savepath}.{varname}'
    else:
        var_path = f'{savepath}.{varname}tp'
    var_ls_results = save_load(var_path, None)
    if var_ls_results == False:
        var_ls_results = calc_func(top_path, traj_path, arms, dims_ls, f'{savepath}.ns', f'{savepath}.sys')
        var_ls_results = save_load(var_path, var_ls_results)
    return var_ls_results, label, plotpath


def ns_plot(data_process_func, results, plot_confs, data, varname):
    '''
    Plot the histogram (distribution) of a single trajectory.
    Plot: frequency vs value.
    :data_process_func: a function that process the obtained measurements into desired float-valued numbers.
    :results: a dictionary/nested list that contains the obtained measurements
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :varname: codename of the measurement.
    '''
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results

    var_ls = data_process_func(var_ls_results, data)
    if min(var_ls) < x_lim[0] or max(var_ls) > x_lim[1]:
        print(f'ns_plot out of range: {min(var_ls)} < {x_lim[0]} : {min(var_ls) < x_lim[0]} or {max(var_ls)} > {x_lim} : {max(var_ls) > x_lim[1]}')
        # assert 0 == 1
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

    # moments
    n, m1, std, m3_s = moments_calc(n, var_ls)

    ax = plt.axes()
    n_sum = np.sum(n)
    plt.errorbar(bin_centers,np.divide(np.array(n),n_sum),yerr = np.divide(np.sqrt(np.array(n)),n_sum),marker = '.',drawstyle = 'steps-mid',c='#4994FF')
    plt.title(label, fontsize=18)
    # customization
    txt = (
        rf'$N={np.sum(n):.2f}$'
        '\n'
        rf'$\mu={m1:.2f}$'
        '\n'
        rf'$\sigma={std:.2f}$'
        '\n'
        rf'$\gamma={m3_s:.2f}$'
    )
    text_x, text_y = text_pos
    plt.text(text_x, text_y, txt, c='#4994FF', fontsize=18) # 28%
    plt.xlabel(x_var, fontsize=18)
    plt.ylabel('Frequencies', fontsize=18)
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    chkdir(os.path.dirname(plotpath))
    plt.savefig(plotpath)
     # plt.show()
    plt.close()
    plt.clf()
    return m1,std,m3_s


def ns_time_js_plot(data_process_func, results, plot_confs, data, varname):
    time_window_width = 11 # 15
    assert time_window_width % 2 == 1 # must be odd
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    var_ls = data_process_func(var_ls_results, data, vtime = True)
    row,col = (1,3)
    fig = plt.figure(figsize=(3*3*col+1,3.5*row))
    gs = fig.add_gridspec(row, col, hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
    axs = gs.subplots(sharex='col')
    # loop over row should be here
    # loop col
    v_ls = [var_tp[1] for var_tp in var_ls]
    time_ls = [var_tp[0] for var_tp in var_ls]
    time_idx = time_ls[time_window_width//2:-time_window_width//2+1]
    ang_runningavg = [np.sum(v_ls[j:j+time_window_width-1])/time_window_width for j in range(len(time_idx))]
    ang_runningstd = [np.std(v_ls[j:j+time_window_width-1]) for j in range(len(time_idx))]
    c = '#4994FF'
    axs[0].scatter(time_idx,ang_runningavg, label = f'running_avg', s=4, c=c) # ,c=patch_color_list[i]
    axs[0].set_xlabel('Time') #, fontsize=18
    axs[0].set_ylabel(f'{x_var} (Avg, wid:{time_window_width})') #, fontsize=18
    l = axs[0].legend(loc = 'lower right')
    axs[1].scatter(time_idx,ang_runningstd, label = f'running_rmsd', s=4, c=c) # ,c=patch_color_list[i]
    axs[1].set_xlabel('Time') #, fontsize=18
    axs[1].set_ylabel(f'{x_var} (STD, wid:{time_window_width})') #, fontsize=18
    l = axs[1].legend(loc = 'lower right')
    axs[2].scatter(time_ls,v_ls, label = f'raw value', s=4) # ,c=patch_color_list[i]
    axs[2].set_xlabel('Time') #, fontsize=18
    axs[2].set_ylabel(x_var) #, fontsize=18
    axs[2].legend(loc = 'lower right')
    chkdir(os.path.dirname(plotpath))
    plt.savefig(plotpath,dpi=400)
    # plt.show()
    plt.close()    
    return True

def ns_time_k2_plot(data_process_func, results, plot_confs, data, varname):
    assert varname in ['k2']
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    var_dic = data_process_func(var_ls_results, data, vtime = True)
    # 0: d1<l2, d2<l2; 1: d1>l2, d2<l2; 2: d1<l2, d2>l2; 3: d1>l2, d2>l2
    shape_state_ls = [1 if (var_dic['deltas'][i][0]>var_dic['lambdas'][i][1]*0.5) and not (var_dic['deltas'][i][1]>var_dic['lambdas'][i][1]*0.5) else 2 if (var_dic['deltas'][i][1]>var_dic['lambdas'][i][1]*0.5) and not (var_dic['deltas'][i][0]>var_dic['lambdas'][i][1]*0.5) else 3 if (var_dic['deltas'][i][0]>var_dic['lambdas'][i][1]*0.5) and (var_dic['deltas'][i][1]>var_dic['lambdas'][i][1]*0.5) else 0 for i in range(len(var_dic['deltas']))]
    lbds0 = [var_dic['lambdas'][i][0] for i in range(len(var_dic['deltas']))]
    avg0 = sum(lbds0)/len(lbds0)
    lbds1 = [var_dic['lambdas'][i][1] for i in range(len(var_dic['deltas']))]
    avg1 = sum(lbds1)/len(lbds1)
    lbds2 = [var_dic['lambdas'][i][2] for i in range(len(var_dic['deltas']))]
    avg2 = sum(lbds2)/len(lbds2)
    std2 = np.std(np.array(lbds2) / avg0)
    std1 = np.std(np.array(lbds1) / avg0)
    std0 = np.std(np.array(lbds0) / avg0)    
    avg2 /= avg0
    avg1 /= avg0
    avg0t = avg0
    avg0 /= avg0
    # with open('tmp\\tmp.txt','a') as f:
    #     f.write(f'{label}. avg_lbd0:{avg0:.1f} ± {std0:.1f}, avg_lbd1:{avg1:.1f} ± {std1:.1f}, avg_lbd2:{avg2:.1f} ± {std2:.1f} , divider lbd0: {avg0t:.4f} \n')
    return True

def ns_heatmap_k2_plot(data_process_func, results, plot_confs, data, varname):
    assert varname in ['k2']
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    time_window_width, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd, ns_struc = get_params(int(label[0]))    
    var_dic = data_process_func(var_ls_results, data, vtime = True)
    # read stacking
    path = plotpath.split('/')
    p = path[3].split('_')
    p[0] = path[2] = 'pj'
    path[3] = '_'.join(p)
    pj_path = '/'.join(path)
    with open(os.path.splitext(pj_path)[0]+'.stack','rb') as f:
        nonstacking_vtime_dic, stacking_vtime_dic = pickle.load(f)
    # determine if the NS is stacking? NO. Plot everything.
    # load data
    is_high_prop_stacking = True
    lbd_axs_ls = var_dic['axes']
    l2_arr = np.array([conf_lbd[1][0]/conf_lbd[0][0] for conf_lbd in lbd_axs_ls])
    l3_arr = np.array([conf_lbd[2][0]/conf_lbd[1][0] for conf_lbd in lbd_axs_ls])
    draw_k2_heatmap(is_high_prop_stacking, True, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath)
    # draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=0)
    # draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=1)
    # draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=2)
    # draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=0, is_reverse_select_number_stack=True)
    # determine if the NS is stacking
    # time_idx = stacking_vtime_dic[(0,1)]['t']
    # stacking_state_arr = np.zeros(len(time_idx))
    # for idx in ns_struc['linked_PA']:
    #     stacking_state_arr += np.array(stacking_vtime_dic[idx]['bool'], dtype=int)
    # num_stacking = len(stacking_state_arr) - sum([1 if b == 0 else 0 for b in stacking_state_arr])
    # num_double_stacking = sum([1 if b == 2 else 0 for b in stacking_state_arr])
    # is_high_prop_stacking = True if num_stacking/len(stacking_state_arr) > 0.122 else False
    # # load data
    # lbd_axs_ls = var_dic['axes']
    # lbd1_ls = [conf_lbd[0][0] for conf_lbd in lbd_axs_ls]
    # l2_ls = [conf_lbd[1][0]/conf_lbd[0][0] for conf_lbd in lbd_axs_ls]
    # l3_ls = [conf_lbd[2][0]/conf_lbd[1][0] for conf_lbd in lbd_axs_ls]
    # # for conf_lbd in lbd_axs_ls:
    # #     (lbd1,axs1),(lbd2,axs2),(lbd3,axs3) = conf_lbd
    # #     lbd1_ls.append(lbd1)
    # #     l2_ls.append(lbd2/lbd1)
    # #     l3_ls.append(lbd3/lbd2)
    # l2_arr = np.array(l2_ls)
    # l3_arr = np.array(l3_ls)
    # if is_high_prop_stacking:
    #     draw_k2_heatmap(is_high_prop_stacking, True, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath)
    #     draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=0)
    #     draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=1)
    #     draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=2)
    #     draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=0, is_reverse_select_number_stack=True)
    # else:
    #     draw_k2_heatmap(is_high_prop_stacking, True, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath)
    #     draw_k2_heatmap(is_high_prop_stacking, False, l2_arr, l3_arr, stacking_vtime_dic, ns_struc, plotpath, number_stack=0)
    return True

def draw_k2_heatmap(is_high_prop_stacking, is_draw_all, l2_arr, l3_arr, stack_dict, ns_struc, plotpath, number_stack = 2, is_reverse_select_number_stack = False):
    # title_high_prop = 'Stacking' if is_high_prop_stacking else 'No-stacking'
    title_high_prop = 'FWHM'
    if is_draw_all:
        is_containing_stack = np.ones(len(stack_dict[(0,1)]['bool']),dtype=bool) # draw everything
    else:
        is_containing_stack = np.zeros(len(stack_dict[(0,1)]['bool']),dtype=int)
        for arm_idx in ns_struc['PA']: # not 'linked_PA': no longer assume linked
            is_containing_stack += np.array(stack_dict[arm_idx]['bool'],dtype=int)
        is_containing_stack = (is_containing_stack==number_stack)
    if len(is_containing_stack) != len(l2_arr): # if stack_dict comes from pj_protocol, its length is [5:-5] instead of [:]
        # l2_arr = l2_arr[5:-5]
        # l3_arr = l3_arr[5:-5]
        is_containing_stack = np.append(is_containing_stack,[is_containing_stack[-1],is_containing_stack[-1],is_containing_stack[-1],is_containing_stack[-1],is_containing_stack[-1]])
        is_containing_stack = np.insert(is_containing_stack,0,[is_containing_stack[0],is_containing_stack[0],is_containing_stack[0],is_containing_stack[0],is_containing_stack[0]])
    if not is_reverse_select_number_stack:
        l2_arr = l2_arr[is_containing_stack]
        l3_arr = l3_arr[is_containing_stack]
        title_not_reverse = ''
    else:
        l2_arr = l2_arr[~is_containing_stack]
        l3_arr = l3_arr[~is_containing_stack]        
        title_not_reverse = 'NOT'
    if sum(is_containing_stack) < 10:
        return False
    import scipy.stats as stat
    xmin = 1
    xmax = 3
    ymin = 1
    ymax = 2
    X, Y = np.mgrid[xmin:xmax:3000j, ymin:ymax:3000j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([l2_arr, l3_arr])
    try:
        kernel = stat.gaussian_kde(values,bw_method=0.2)
    except:
        if is_draw_all:
            savefig_draw_num_stackings = 'All'
        else:
            savefig_draw_num_stackings = number_stack
        f=open(os.path.splitext(plotpath)[0]+'-Heatmap'+f'-{title_high_prop}'+f'-{title_not_reverse}{savefig_draw_num_stackings}'+'.txt','w')
        f.write("Woops! Singular matrix! Maybe too few points.")
        f.close()
        return False
    Z = np.reshape(kernel(positions).T, X.shape)
    Z = np.rot90(Z)    
    Z_68 = (generate_Z_layer(Z,0.68,is_FWHM=True).astype(bool)) # Z_FWHM
    Z_rings = Z_68
    Z_draw = Z*((~Z_rings).astype(int))    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.imshow(Z_draw, cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
    # ax.scatter(l2_arr, l3_arr, s=0.4,c='#FFFFFF')
    ax.set_xlim([1, 3])
    ax.set_ylim([1, 2])
    if is_draw_all:
        ax.set_title(f'{title_high_prop} Confs, all confs, {title_not_reverse} {is_containing_stack.sum()} pts')
        savefig_draw_num_stackings = 'All'
    else:
        ax.set_title(f'{title_high_prop}, Only confs w/ {title_not_reverse} {number_stack} stackings, {is_containing_stack.sum()} pts')
        savefig_draw_num_stackings = number_stack
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.tick_params(bottom=True,top=True,left=True,right=True,direction='in')
    ax.tick_params(bottom=True,top=True,left=True,right=True,direction='in',which='minor')
    plt.savefig(os.path.splitext(plotpath)[0]+'-Heatmap'+f'-{title_high_prop}'+f'-{title_not_reverse}{savefig_draw_num_stackings}'+'.png',dpi=400)    
    return True

def generate_Z_layer(Z,prob,delta=0.05,is_FWHM=False):
    if is_FWHM:
        Z_threshold = (np.max(Z))/2 # delta at center now
    else:
        prob = 1-prob
        Z_sum = np.sum(Z)
        Z_thresholds_arr = np.linspace(np.min(Z),np.max(Z),1000)
        Z_threshold = 0
        for Z_th in Z_thresholds_arr:
            if np.sum(Z*(Z>Z_th)) / Z_sum <= prob: # integral of values above Z_th drops to 1-prob
                Z_threshold = Z_th
                break
    # Produce Z
    Z_sub = Z*(Z>Z_threshold-delta/2)
    Z_sub_ring = Z_sub*(~(Z_sub>Z_threshold+delta/2))
    return Z_sub_ring

def ns_time_si_plot(data_process_func, results, plot_confs, data, varname):
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    # obtain stacking_vtime
    stacking_vtime_dic, nonstacking_vtime_dic = data_process_func(var_ls_results, data, vtime=True)
    path = plotpath.split('/')
    if not os.path.exists('/'.join(path[:-1])):
        os.makedirs('/'.join(path[:-1]))
    p = path[3].split('_')
    p[0] = path[2] = 'pj'
    path[3] = '_'.join(p)
    pj_path = '/'.join(path)
    with open(os.path.splitext(pj_path)[0]+'.stack','wb') as f:
        print('si write to: ' + os.path.splitext(pj_path)[0]+'.stack')
        pickle.dump((nonstacking_vtime_dic,stacking_vtime_dic), f) # legacy: pj 'nonstacking' --> pa 'stacking'; si corrected this issue
    return True

def ns_time_pa_plot(data_process_func, results, plot_confs, data, varname):
    '''
    Plot the value (patch angle vtime) vs. time plot of a single trajectory.
    Plot: value vs time.
    :data_process_func: a function that process the obtained measurements into desired float-valued numbers.
    :results: a dictionary/nested list that contains the obtained measurements
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    :varname: codename of the measurement.
    # Arm color schema: Arm_0 ~ #4994FF (blue), Arm_1 ~ #FCC555 (yellow), Arm_2 ~ #7AA77A, Arm_3 ~ #E55050
    # Patch_color_list order: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    Strand-sharing: Arm_0 ~ strands: 1 & 4 ; Arm_1 ~ strands: 4 & 3 ; Arm_2 ~ strands: 3 & 2 ; Arm_3 ~ strands: 2 & 1
    '''
    # compare_list = [((0,1),(2,3)),((0,3),(1,2)),((0,2),(1,3)),((0,1),(0,3)),((0,1),(0,2)),((0,3),(0,2))] # former 3 expected to be similar; latter 3 are 3c2=3
    varname, x_var, x_lim, y_lim, text_pos, bin_num = plot_confs
    var_ls_results, label, plotpath = results
    time_window_width, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd, ns_struc = get_params(int(label[0]))
    stacking_vtime_dic = {} # only for PJ!
    nonstacking_vtime_dic = {} # only for PJ!    
    # only for patch angle vtime.
    assert varname in ['pa','pj']
    var_dic = data_process_func(var_ls_results, data, vtime = True)
    row,col = (len((var_dic.items())),3)
    fig = plt.figure(figsize=(3*3*col+1,3.5*row))
    gs = fig.add_gridspec(row, col, hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
    axs = gs.subplots(sharex='col')
    # reading is_stacking from pj folder. overriding the crit settings.
    if varname == 'pa':
        path = plotpath.split('/')
        p = path[3].split('_')
        p[0] = path[2] = 'pj'
        path[3] = '_'.join(p)
        pj_path = '/'.join(path)
        with open(os.path.splitext(pj_path)[0]+'.stack','rb') as f:
            nonstacking_vtime_dic, stacking_vtime_dic = pickle.load(f) # legacy: pj 'nonstacking' --> pa 'stacking'
    elif varname == 'pj':
        with open(os.path.splitext(plotpath)[0]+'.stack','rb') as f:
            nonstacking_vtime_dic, stacking_vtime_dic = pickle.load(f) # legacy: pj 'nonstacking' --> pa 'stacking'
    # reading ends
    for i, ((ia1,ia2), ang_vtime_dic) in enumerate(var_dic.items()):
        # i: loop of arm patches.
        if (ia1,ia2) in ns_struc['linked_PA']:
            legend_color='#F7AED0'# #E55050
        else:
            legend_color='#4994FF'
        ang_ls = [ang_vtime_dic[t][0] for t in ang_vtime_dic if type(t) == int]
        time_ls = [t for t in ang_vtime_dic if type(t) == int]
        time_idx = time_ls[time_window_width//2:-time_window_width//2+1]
        ang_runningavg = [np.sum(ang_ls[j:j+time_window_width])/time_window_width for j in range(len(time_idx))]
        ang_runningstd = [np.std(ang_ls[j:j+time_window_width]) for j in range(len(time_idx))]
        if varname == 'pj':
            # avg_running_is_stacking, avg_running_is_nonstacking = identifying_stacking(ang_runningavg, ang_runningstd, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd)
            # stacking_vtime_dic[(ia1,ia2)] = {'bool':avg_running_is_stacking,'t':time_idx,'val':ang_runningavg}
            # nonstacking_vtime_dic[(ia1,ia2)] = {'bool':avg_running_is_nonstacking,'t':time_idx,'val':ang_runningavg}
            avg_running_is_stacking = stacking_vtime_dic[(ia1,ia2)]['bool']
            avg_running_is_nonstacking = nonstacking_vtime_dic[(ia1,ia2)]['bool']
            if time_idx != stacking_vtime_dic[(ia1,ia2)]['t']:
                # coming from 'si'
                avg_running_is_stacking = stacking_vtime_dic[(ia1,ia2)]['bool'][5:-5]
                avg_running_is_nonstacking = nonstacking_vtime_dic[(ia1,ia2)]['bool'][5:-5]                                  
            stacking_vtime_dic[(ia1,ia2)]['val'] = ang_runningavg
            nonstacking_vtime_dic[(ia1,ia2)]['val'] = ang_runningavg
            stacking_vtime_dic[(ia1,ia2)]['raw'] = ang_ls[time_window_width//2:-time_window_width//2+1]
            nonstacking_vtime_dic[(ia1,ia2)]['raw'] = ang_ls[time_window_width//2:-time_window_width//2+1]
        elif varname == 'pa':
            avg_running_is_stacking = stacking_vtime_dic[(ia1,ia2)]['bool']
            avg_running_is_nonstacking = nonstacking_vtime_dic[(ia1,ia2)]['bool']
            if time_idx != stacking_vtime_dic[(ia1,ia2)]['t']:
                # coming from 'si'
                avg_running_is_stacking = stacking_vtime_dic[(ia1,ia2)]['bool'][5:-5]
                avg_running_is_nonstacking = nonstacking_vtime_dic[(ia1,ia2)]['bool'][5:-5]                                  
            stacking_vtime_dic[(ia1,ia2)]['val'] = ang_runningavg
            nonstacking_vtime_dic[(ia1,ia2)]['val'] = ang_runningavg
            stacking_vtime_dic[(ia1,ia2)]['raw'] = ang_ls[time_window_width//2:-time_window_width//2+1]
            nonstacking_vtime_dic[(ia1,ia2)]['raw'] = ang_ls[time_window_width//2:-time_window_width//2+1]
        # transiting_vtime_dic[(ia1,ia2)] = list(filter(None,[time_idx[i] if not(a[i] or b[i]) else None for i in range(len(c))]))
        # plot
        c = ['#00FF00' if is_stacking == True else '#000000' for is_stacking in avg_running_is_stacking] # stacking: dirty
        # c = ['#FF0000' if avg_running_is_stacking[i] == True else c[i] for i in range(len(c))] # nonstacking
        axs[i,0].scatter(time_idx,ang_runningavg, label = f'Arm{ia1}~Arm{ia2}', s=4, c=c) # ,c=patch_color_list[i]
        # axs[i,0].plot([min(time_idx),min(time_idx)+nonstacking_min_length*100000],[90,90],c='#FF0000', linewidth=1.5)
        # axs[i,0].plot([min(time_idx),max(time_idx)],[90,90],c='#000000', linewidth=0.5)
        # axs[i,0].plot([min(time_idx),max(time_idx)],[nonstacking_crit_ang,nonstacking_crit_ang],c='#FF0000', linewidth=0.5)
        # axs[i,0].set_xlabel('Time') #, fontsize=18
        axs[i,0].set_ylim((0,180))
        axs[i,0].set_ylabel(f'{x_var} (Avg, wid:{time_window_width})') #, fontsize=18
        l = axs[i,0].legend(loc = 'lower right')
        for text in l.get_texts():
            text.set_color(legend_color)
        axs[i,1].scatter(time_idx,ang_runningstd, label = f'Arm{ia1}~Arm{ia2}', s=4, c=c) # ,c=patch_color_list[i]
        axs[i,1].set_xlabel('Time') #, fontsize=18
        axs[i,1].plot([min(time_idx),max(time_idx)],[nonstacking_crit_rmsd,nonstacking_crit_rmsd],c='#FF0000', linewidth=0.5)
        axs[i,1].set_ylabel(f'{x_var} (STD, wid:{time_window_width})') #, fontsize=18
        l = axs[i,1].legend(loc = 'lower right')
        for text in l.get_texts():
            text.set_color(legend_color)
    for i, ((ia1,ia2), ang_vtime_dic) in enumerate(var_dic.items()):
        ang_ls = [ang_vtime_dic[t][0] for t in ang_vtime_dic if type(t) == int]
        time_ls = [t for t in ang_vtime_dic if type(t) == int]
        axs[i,2].scatter(time_ls,ang_ls, label = f'Arm{ia1}~Arm{ia2}', s=4, c= ['#00FF00' if is_stacking == True else '#000000' for is_stacking in stacking_vtime_dic[(ia1,ia2)]['bool']]) # ,c=patch_color_list[i]
        # axs[i,2].plot([min(time_ls),max(time_ls)],[90,90],c='#000000')
        # axs[i,2].plot([min(time_ls),max(time_ls)],[180,180],c='#000000')
        # axs[i,2].plot([min(time_ls),max(time_ls)],[130,130],c='#000000')
        # axs[i,2].set_xlabel('Time') #, fontsize=18
        axs[i,2].set_ylabel(x_var) #, fontsize=18
        axs[i,2].set_ylim((0,180))
        l = axs[i,2].legend(loc = 'lower right')
        for text in l.get_texts():
            text.set_color(legend_color)
    # for i, (idx1, idx2) in enumerate(compare_list):
    #     ang_vtime1 = var_dic[idx1]
    #     ang_vtime2 = var_dic[idx2]
    #     time_ls = [t for t in ang_vtime1 if type(t) == int]
    #     ang_diff_ls = [ang_vtime1[t][0] - ang_vtime2[t][0] for t in time_ls]
    #     axs[i,3].plot(time_ls,ang_diff_ls, label = f'Plot{idx1} - Plot{idx2}')
    #     axs[i,3].set_xlabel('Time') #, fontsize=18
    #     axs[i,3].set_ylabel(f'{x_var}, Diff') #, fontsize=18
    #     axs[i,3].legend(loc = 'lower right')
    # plt.xlim(x_lim)
    # plt.ylim(y_lim)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    # plotpath = plotpath.split('.')
    # plotpath[-2] = plotpath[-2]+'_vtime'
    # plotpath = '.'.join(plotpath)
    total_count = counting_stacking(stacking_vtime_dic, 100, True, False, ns_struc) # pj: nonstacking--->stacking
    stacking_0_count = counting_stacking(stacking_vtime_dic, 0, False, False, ns_struc)
    stacking_not0_count = counting_stacking(stacking_vtime_dic, 0, False, True, ns_struc)
    stacking_1_count = counting_stacking(stacking_vtime_dic, 1, False, False, ns_struc)
    stacking_2_count = counting_stacking(stacking_vtime_dic, 2, False, False, ns_struc)
    stacking_3_count = counting_stacking(stacking_vtime_dic, 3, False, False, ns_struc)
    plt.suptitle(f'Total:{total_count}, 0_stk:{stacking_0_count}, Not0_stk:{stacking_not0_count}, 1_stk: {stacking_1_count}, 2_stk: {stacking_2_count}, 3_stk: {stacking_3_count}')
    chkdir(os.path.dirname(plotpath))
    plt.savefig(plotpath,dpi=400)
     # plt.show()
    plt.close()
    # with open(os.path.splitext(plotpath)[0]+'.stack','wb') as f:
    #     # if varname == 'pj': # definition of stacking in pj is the reverse.varname
    #     #     pickle.dump((nonstacking_vtime_dic,stacking_vtime_dic), f)
    #     # else:
    #     #     pass
    #     pickle.dump((stacking_vtime_dic,nonstacking_vtime_dic), f)
    return True

def identifying_stacking(ang_runningavg, ang_runningstd, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd):
    avg_running_is_stacking = [all(ang > stacking_crit_ang for ang in ang_runningavg[i:i+stacking_min_length]) and all(rmsd < stacking_crit_rmsd for rmsd in ang_runningstd[i:i+stacking_min_length]) for i in range(len(ang_runningavg)-stacking_min_length+1)] # 4arm only!
    avg_running_is_stacking2 = [any(cond==True for cond in avg_running_is_stacking[i:i+stacking_min_length]) for i in range(len(avg_running_is_stacking)-stacking_min_length+1)]
    avg_running_is_stacking2[0:0] = [avg_running_is_stacking[i] if True not in avg_running_is_stacking[0:i] else True for i in range(stacking_min_length-1)]
    avg_running_is_stacking2.extend([any(cond==True for cond in avg_running_is_stacking[-stacking_min_length+i+1:]) for i in range(stacking_min_length-1)])
    # non-stacking
    avg_running_is_nonstacking = [all(ang < nonstacking_crit_ang for ang in ang_runningavg[i:i+nonstacking_min_length]) and all(rmsd < nonstacking_crit_rmsd for rmsd in ang_runningstd[i:i+stacking_min_length]) for i in range(len(ang_runningavg)-nonstacking_min_length+1)] # 4arm only!
    avg_running_is_nonstacking2 = [any(cond==True for cond in avg_running_is_nonstacking[i:i+nonstacking_min_length]) for i in range(len(avg_running_is_nonstacking)-nonstacking_min_length+1)]
    avg_running_is_nonstacking2[0:0] = [avg_running_is_nonstacking[i] if True not in avg_running_is_nonstacking[0:i] else True for i in range(nonstacking_min_length-1)]
    avg_running_is_nonstacking2.extend([any(cond==True for cond in avg_running_is_nonstacking[-nonstacking_min_length+i+1:]) for i in range(nonstacking_min_length-1)])
    return avg_running_is_stacking2, avg_running_is_nonstacking2


def counting_stacking(stack_dict, number_stack, is_count_all, is_reverse_select_number_stack, ns_struc):
    if is_count_all:
        # is_containing_stack = np.array(stack_dict[(0,1)]['bool']) + np.array(stack_dict[(1,2)]['bool']) + np.array(stack_dict[(2,3)]['bool']) + np.array(stack_dict[(0,3)]['bool']) # selecting w/ stacking
        is_containing_stack_arr = np.ones(len(stack_dict[(0,1)]['bool']),dtype=bool) # draw everything
    else:
        is_containing_stack_arr = np.zeros(len(stack_dict[(0,1)]['bool']),dtype=int)
        for iaidx in ns_struc['PA']: # not 'linked_PA': no longer assume linked
            is_containing_stack_arr += np.array(stack_dict[iaidx]['bool'],dtype=int)
        is_containing_stack_arr = is_containing_stack_arr == number_stack
    if not is_reverse_select_number_stack:
        return np.sum(is_containing_stack_arr)
    else:
        return np.sum(~is_containing_stack_arr)     


    stacking_result = {}
    # stacking in one series
    for iaidx in ns_struc['linked_PA']:
        s = stacking_vtime_dic[iaidx]['bool']
        stacking_result[iaidx] = [sum(s),len(s),sum(s)/len(s)]
    # simutaneous stacking
    if ns_struc['#arm'] == 4:
        for iaidx1, iaidx2 in ns_struc['pairing_linked']:
            s1 = stacking_vtime_dic[iaidx1]['bool']
            s2 = stacking_vtime_dic[iaidx2]['bool']
            compare_result = [True if s1[i]==True and s2[i]==True else False for i in range(len(s1))]
            stacking_result[(iaidx1,iaidx2)] = [sum(compare_result),len(compare_result),sum(compare_result)/len(compare_result)]
    return stacking_result

def get_params(arm_num):
    # time_window_width = 11 # 15
    # stacking_min_length = 60 # 4arm only!
    # stacking_crit_ang = 120 # 4arm only!
    # stacking_vtime_dic = {} # only for PA!
    # nonstacking_min_length = 60 # 4arm only!
    # nonstacking_crit_ang = 120 # 4arm only!
    # nonstacking_vtime_dic = {} # only for PA! 
    # ns_struc = ns_struc_3arm if int(label[0]) == 3 else ns_struc_4arm if int(label[0]) == 4 else ns_struc_5arm if int(label[0]) == 5 else ns_struc_6arm if int(label[0]) == 6 else False   
    # PJ use only!
    # if arm_num == 3:
    #     time_window_width = 7 # 15
    #     stacking_min_length = 10 # 4arm only!
    #     stacking_crit_ang = 120 # 4arm only!
    #     stacking_crit_rmsd = 9 # 4arm only!
    #     nonstacking_min_length = 10 # 4arm only!
    #     nonstacking_crit_ang = 105 # 4arm only!
    #     nonstacking_crit_rmsd = 12 # 4arm only!
    #     ns_struc = {'#arm':3, 'linked_PA': [(0,1),(0,2),(1,2)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    # elif arm_num == 4:
    #     time_window_width = 11 # 15
    #     stacking_min_length = 1 # 4arm only!
    #     stacking_crit_ang = 110 # 4arm only!
    #     stacking_crit_rmsd = 13 # 4arm only!
    #     nonstacking_min_length = 1 # 4arm only!
    #     nonstacking_crit_ang = 65 # 4arm only!
    #     nonstacking_crit_rmsd = 10 # 4arm only!
    #     ns_struc = {'#arm':4, 'linked_PA': [(0,1),(0,3),(1,2),(2,3)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    # elif arm_num == 5:
    #     time_window_width = 11 # 15
    #     stacking_min_length = 15 # 4arm only!
    #     stacking_crit_ang = 145 # 4arm only! # should be 130?
    #     stacking_crit_rmsd = 7 # 4arm only!
    #     nonstacking_min_length = 15 # 4arm only!
    #     nonstacking_crit_ang = 65 # 4arm only!
    #     nonstacking_crit_rmsd = 10 # 4arm only!
    #     ns_struc = {'#arm':5, 'linked_PA': [(0,1),(0,4),(1,2),(2,3),(3,4)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    # elif arm_num == 6:
    #     time_window_width = 11 # 15
    #     stacking_min_length = 15 # 4arm only!
    #     stacking_crit_ang = 110 # 4arm only!
    #     stacking_crit_rmsd = 7 # 4arm only!
    #     nonstacking_min_length = 15 # 4arm only!
    #     nonstacking_crit_ang = 50 # 4arm only!
    #     nonstacking_crit_rmsd = 10 # 4arm only!
    #     ns_struc = {'#arm':6, 'linked_PA': [(0,1),(0,5),(1,2),(2,3),(3,4),(4,5)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    # else:
    #     assert 0==1
    if arm_num == 3:
        time_window_width = 11 # 15
        stacking_min_length = 10 # 
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 105 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'#arm':3, 'linked_PA': [(0,1),(0,2),(1,2)], 'PA' : [(0,1),(0,2),(1,2)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 4:
        time_window_width = 11 # 15
        stacking_min_length = 10 # 
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 65 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'#arm':4, 'linked_PA': [(0,1),(0,3),(1,2),(2,3)], 'PA': [(0,1),(0,3),(1,2),(2,3),(0,2),(1,3)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 5:
        time_window_width = 11 # 15
        stacking_min_length = 10 # useful when CoM is at the center of PJ.
        stacking_crit_ang = 120 #  # should be 130?
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 65 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'#arm':5, 'linked_PA': [(0,1),(0,4),(1,2),(2,3),(3,4)], 'PA' : [(0,1),(0,4),(1,2),(2,3),(3,4),(0,2),(1,3),(2,4),(0,3),(1,4)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 6:
        time_window_width = 11 # 15
        stacking_min_length = 10 # useful when CoM is at the center of PJ.
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 50 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'#arm':6, 'linked_PA': [(0,1),(0,5),(1,2),(2,3),(3,4),(4,5)], 'PA' : [(0,1),(0,5),(1,2),(2,3),(3,4),(4,5),(0,2),(1,3),(2,4),(3,5),(0,3),(1,4),(2,5),(0,4),(1,5)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    else:
        assert 0==1    
    assert time_window_width % 2 == 1 # must be odd
    return time_window_width, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd, ns_struc