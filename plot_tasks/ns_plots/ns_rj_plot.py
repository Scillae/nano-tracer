from calc_tasks import RoG_jun_calc
from utils.ns_plot import SL_ns, ns_plot
from utils.tools import dims_adjust


def data_process_func(RoG_ls_res, data):
    # arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    RoG_ls = [i[1] for i in RoG_ls_res]
    return RoG_ls


def ns_rj_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    varname = 'rj'
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix)
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(RoG_jun_calc, data, varname) # results: var_ls_results, label, plotpath
    #### plot confs ####
    x_var = 'RoG of the Junction'
    x_lim = (3,12)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7*(x_lim[1]-x_lim[0])+x_lim[0], (0.215/0.3)*(y_lim[1]-y_lim[0]))
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    # ns_plot does save figures.
    return summary