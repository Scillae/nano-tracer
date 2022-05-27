from calc_tasks import stacking_local_identify_calc
from utils.ns_plot import SL_ns, ns_plot, ns_time_si_plot
from utils.tools import dims_adjust


def data_process_func(stacking_res, data, vtime = False): # should be trimmed.
    import numpy as np
    if vtime:
        # get stack_dic and nonstack_dic
        # stacking_res: OrderedDict{t_stamp: is_stacking, stack_ls}
        # stack_dic: 
        arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls = data
        # init
        stacking_vtime_dic = {}
        nonstacking_vtime_dic = {}
        for ia1,ia2 in create_ia_idx_ls(arms):
            stacking_vtime_dic[(ia1,ia2)] = {}
            nonstacking_vtime_dic[(ia1,ia2)] = {}
            stacking_vtime_dic[(ia1,ia2)]['bool'] = np.zeros(len(stacking_res.keys()),dtype=bool)
            nonstacking_vtime_dic[(ia1,ia2)]['bool'] = np.zeros(len(stacking_res.keys()),dtype=bool)
            stacking_vtime_dic[(ia1,ia2)]['t'] = np.zeros(len(stacking_res.keys()),dtype=int)
            nonstacking_vtime_dic[(ia1,ia2)]['t'] = np.zeros(len(stacking_res.keys()),dtype=int)
            stacking_vtime_dic[(ia1,ia2)]['adj_bps'] = np.zeros(len(stacking_res.keys()),dtype=bool)
        # fill data in
        for i, (t, (is_stacking, stack_ls)) in enumerate(stacking_res.items()):
            
            if is_stacking == False:
                continue
            for adj_bp, adj_bp2, arm_id, arm2_id in stack_ls:
                stacking_vtime_dic[(arm_id,arm2_id)]['bool'][i] = True
                stacking_vtime_dic[(arm_id,arm2_id)]['adj_bps'][i] = adj_bp, adj_bp2 # deprecated
        # smoothing
        window_hw = 5 # window_hw *2 +1 == window_width
        for ia1,ia2 in create_ia_idx_ls(arms):
            bool_ls = stacking_vtime_dic[(ia1,ia2)]['bool']
            smooth_ls = [True if sum(bool_ls[i:i+window_hw*2 +1]) > window_hw+1 else False for i in range(len(bool_ls)-window_hw*2)] # threshold: 50%
            r_ls = [True if sum(bool_ls[0:i+window_hw +1]) > window_hw or bool_ls[i] else False for i in range(window_hw)] # still 50% of full window width?
            r_ls.extend(smooth_ls)
            r_ls.extend([True if sum(bool_ls[i-window_hw:]) > window_hw or bool_ls[i] else False for i in range(-window_hw,0)])
            assert len(r_ls) == len(bool_ls)
            stacking_vtime_dic[(ia1,ia2)]['bool'] = r_ls
            stacking_vtime_dic[(ia1,ia2)]['t'] = list(stacking_res) # filling in time indices as well
        return stacking_vtime_dic, nonstacking_vtime_dic
    return True

def create_ia_idx_ls(arm_num):
    return [(ia1,ia2) for ia1 in range(arm_num) for ia2 in range(ia1+1,arm_num)]

def ns_si_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', flag_suffix='', dims_ls= [20,2,7]):
    varname = 'si'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls)
    results = SL_ns(stacking_local_identify_calc, data, varname)
    results_vtime = SL_ns(stacking_local_identify_calc, data, varname, vtime=True)
    # return (0,0,0)
    #### plot confs ####
    x_var = rf'Patch Angles ($^\circ$)'
    x_lim = (0,360)
    y_lim = (0,0.17)
    bin_num = 50
    text_pos = (10, 0.12)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    ns_time_si_plot(data_process_func, results_vtime, plot_confs, data, varname)
    # ns_plot does save figures.
    return (0,0,0)
    return summary