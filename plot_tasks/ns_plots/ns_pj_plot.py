from calc_tasks.patch_jun_calc import patch_jun_calc
from utils.ns_plot import SL_ns, ns_plot, ns_time_pa_plot
from utils.tools import dims_adjust


def data_process_func(p_ang_ls_res, data, vtime = False): # should be trimmed.
    from collections import OrderedDict
    import numpy as np
    p_angs_dic, arms_idx = p_ang_ls_res
    # tracing of specific p-angles setup
    angle_dic = OrderedDict() # {(ia1, ia2):{t_stamp:angle}} 
    for idx_1, ia1 in enumerate(arms_idx):
        for idx_2 in range(idx_1+1, len(arms_idx)):
            ia2 = arms_idx[idx_2]
            angle_dic[(ia1, ia2)] = OrderedDict()
    # hist of p-angle
    drop_t_ls = []
    for t_stamp, ang_ls in p_angs_dic.items():
        ## sanity chk: because of old pairing method. TODO: look into the dropped frames
        s = [1 for i in ang_ls if i[1] == True]
        if len(s) != len(arms_idx): # n arms should yield n angles.
            print(t_stamp)
            drop_t_ls.append(t_stamp)
            continue
        ## chk ends
        for ang, is_sharing, ia_tp in ang_ls:
            # {(ia1, ia2):{t_stamp:angle, 'is_sharing':is_sharing}}
            angle_dic[ia_tp][t_stamp] = [ang] # tracing of specific p-angles
            if 'is_sharing' not in angle_dic[ia_tp].keys():
                angle_dic[ia_tp]['is_sharing'] = is_sharing
    print(f'Total time steps dropped: {len(drop_t_ls)}')
    t_ls = [t for t in angle_dic[(0,1)] if type(t) == int]
    ang_ls = [t_dic[t][0] for t_dic in angle_dic.values() if t_dic['is_sharing']==True for t in t_ls] # pooled for is_sharing
    if vtime == False:
        return ang_ls
    else:
        return angle_dic


def ns_pj_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    varname = 'pj'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(patch_jun_calc, data, varname)
    results_vtime = SL_ns(patch_jun_calc, data, varname, vtime=True)
    #### plot confs ####
    x_var = rf'Patch Angles of Junction ($^\circ$)'
    x_lim = (0,360)
    y_lim = (0,0.17)
    bin_num = 50
    text_pos = (10, 0.12)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    ns_time_pa_plot(data_process_func, results_vtime, plot_confs, data, varname)
    # ns_plot does save figures.
    return summary
        # moments ends
        # # TODO: plot single_batch angle here:
        # '''
        # for ia_tp in angle_dic.keys():
        #     t_ls = list(angle_dic[ia_tp].keys())
        #     a_ls = list(angle_dic[ia_tp].values())
        #     print(f'Angle:{ia_tp}, time: {t_ls} vs angle: {a_ls}')
        # '''
