from calc_tasks import k2_calc
from utils.ns_plot import SL_ns, ns_plot, ns_time_k2_plot
from utils.tools import dims_adjust


def data_process_func(k2_ls_res, data, vtime=False):
    # arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    if vtime:
        t_ls = [i[0] for i in k2_ls_res]
        k2_ls = [i[1] for i in k2_ls_res]
        lbd_ls = [sorted(i[2]) for i in k2_ls_res]
        del_ls = [(i[1]-i[0],i[2]-i[1]) for i in lbd_ls]
        axs_ls = [sorted(zip(i[2],i[3])) for i in k2_ls_res] # [((lbds),(axes))] instead of [(axes)]
        return {'t':t_ls, 'k2':k2_ls, 'lambdas':lbd_ls, 'deltas':del_ls, 'axes':axs_ls}
    else:
        k2_ls = [i[1] for i in k2_ls_res]
        return k2_ls


def ns_k2_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    varname = 'k2'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(k2_calc, data, varname) # results: var_ls_results, label, plotpath
    results_vtime = SL_ns(k2_calc, data, varname, vtime=True)
    #### plot confs ####
    x_var = 'k2 value'
    x_lim = (0,1)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7, 0.215)
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    ns_time_k2_plot(data_process_func, results_vtime, plot_confs, data, varname)
    # ns_plot does save figures.
    return summary