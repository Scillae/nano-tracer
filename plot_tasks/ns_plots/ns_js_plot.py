from calc_tasks import jun_shift_calc
from utils.ns_plot import SL_ns, ns_plot


def data_process_func(js_ls_res, data):
    # arms, temp, conc, sp_suffix, conf_suffix, dims_ls = data
    js_ls = [i[1]*0.8518 for i in js_ls_res]
    return js_ls


def ns_js_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    varname = 'js'
    dims_ls[1] = int(conf_suffix.split('_')[1])
    data = (arms, temp, conc, sp_suffix, conf_suffix, dims_ls)
    results = SL_ns(jun_shift_calc, data, varname) # results: var_ls_results, label, plotpath
    #### plot confs ####
    x_var = 'Junction Shift (nm)'
    x_lim = (0,2)
    y_lim = (0,0.3)
    bin_num = 50
    text_pos = (0.7*(x_lim[1]-x_lim[0])+x_lim[0], (0.215/0.3)*(y_lim[1]-y_lim[0]))
    #### conf ends ####
    plot_confs = varname, x_var, x_lim, y_lim, text_pos, bin_num
    summary = ns_plot(data_process_func, results, plot_confs, data, varname) # summary: m1,std,m3_s
    # ns_plot does save figures.
    return summary