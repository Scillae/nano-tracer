from calc_tasks import patch_angle_calc
from utils.ns_plot import SL_ns, ns_plot


def data_process_func(p_ang_ls_res, data): # should be trimmed.
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
    angle_ls = [] # historical approach.
    t_ls = []        
    drop_t_ls = []
    for t_stamp, ang_ls in p_angs_dic.items():
        ## sanity chk: because of old pairing method.
        s = [1 for i in ang_ls if i[1] == True]
        if len(s) != len(arms_idx): # n arms should yield n angles.
            print(t_stamp)
            drop_t_ls.append(t_stamp)
            continue
        ## chk ends
        for ang, is_sharing, ia_tp in ang_ls:
            if is_sharing == True:  # patch angle: sharing a strand
                angle_ls.append(ang)
                t_ls.append(t_stamp)
                # {(ia1, ia2):{t_stamp:angle}}
                angle_dic[ia_tp][t_stamp] = ang # tracing of specific p-angles
    print(f'Total time steps dropped: {len(drop_t_ls)}')
    ang_arr = np.array(angle_ls)
    ang_deg = np.degrees(ang_arr)
    return ang_deg
        # # TODO: plot single_batch angle:
        # '''
        # for ia_tp in angle_dic.keys():
        #     t_ls = list(angle_dic[ia_tp].keys())
        #     a_ls = list(angle_dic[ia_tp].values())
        #     print(f'Angle:{ia_tp}, time: {t_ls} vs angle: {a_ls}')
        # '''

        # plt.xticks(np.linspace(0,180,num=10))
        # plt.yticks(np.linspace(0,0.16,num=9))


def ns_pa_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
   
    varname = 'pa'
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(patch_angle_calc, data, varname)
    #### plot confs ####
    x_var = rf'Patch Angles ($^\circ$)'
    x_lim = (0,180)
    y_lim = (0,0.17)
    bin_num = 50
    text_pos = (10, 0.12)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
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
