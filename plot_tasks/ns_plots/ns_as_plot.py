from calc_tasks import arm_stiffness_calc
from utils.ns_plot import SL_ns, ns_plot


def data_process_func(as_results, data):
    as_dic, arms_idx = as_results
    from collections import OrderedDict
    stf_dic = OrderedDict() # {arm_id:{t:armstf}}
    for ia in arms_idx:
        stf_dic[ia] = OrderedDict()
    for t_stamp, as_ls in as_dic.items():
        for armstf, ia in as_ls:
            stf_dic[ia][t_stamp] = armstf
    stf_ls = [stf for s_a_dic in list(stf_dic.values()) for stf in list(s_a_dic.values())] # s_a_dic: {t:armstf} of a single arm; stf_ls:[stf]
    return stf_ls



def ns_as_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):

    varname = 'as'
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(arm_stiffness_calc, data, varname) # results: var_ls_results, label, plotpath
    #### plot confs ####
    x_var = 'Arm Stiffness'
    x_lim = (0,1)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7, 0.215)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    # ns_plot does save figures.
    return summary