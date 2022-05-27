from calc_tasks import RoG_jun_calc
from utils.ns_plot import SL_ns, ns_plot, ns_time_js_plot
from utils.tools import dims_adjust


def data_process_func(RoG_ls_res, data, vtime = False):
    # arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls = data
    RoG_ls = [i[1]*0.8518 for i in RoG_ls_res]
    if vtime == False:
        return RoG_ls
    else:
        return RoG_ls_res


def ns_rj_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', flag_suffix='', dims_ls= [20,2,7]):
    varname = 'rj'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls)
    results = SL_ns(RoG_jun_calc, data, varname) # results: var_ls_results, label, plotpath
    results_vtime = SL_ns(RoG_jun_calc, data, varname, vtime=True)
    #### plot confs ####
    x_var = 'R of Gyr of the Junction (nm)'
    x_lim = (0,4.5)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7*(x_lim[1]-x_lim[0])+x_lim[0], (0.215/0.3)*(y_lim[1]-y_lim[0]))
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    ns_time_js_plot(data_process_func, results_vtime, plot_confs, data, varname)
    # ns_plot does save figures.
    return summary