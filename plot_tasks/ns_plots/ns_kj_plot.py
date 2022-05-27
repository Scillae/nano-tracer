from calc_tasks import k2_jun_calc
from utils.ns_plot import SL_ns, ns_plot, ns_time_js_plot
from utils.tools import dims_adjust


def data_process_func(k2_ls_res, data, vtime = False):
    # arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls = data
    k2_ls = [i[1] for i in k2_ls_res]
    if vtime == False:
        return k2_ls
    else:
        return k2_ls_res


def ns_kj_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', flag_suffix='', dims_ls= [20,2,7]):
    varname = 'kj'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls)
    results = SL_ns(k2_jun_calc, data, varname) # results: var_ls_results, label, plotpath
    results_vtime = SL_ns(k2_jun_calc, data, varname, vtime=True)
    #### plot confs ####
    x_var = 'k2 of Junction'
    x_lim = (0,1)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7, 0.215)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    ns_time_js_plot(data_process_func, results_vtime, plot_confs, data, varname)
    # ns_plot does save figures.
    return summary