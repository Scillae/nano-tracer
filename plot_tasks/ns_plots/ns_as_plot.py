from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def ns_as_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    # path & str define
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    if single == True:
        top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
        traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
    else:
        print('TODO! RETURNING...')
        return
    # path & str end
    # saveload
    savepath = f'data/composed_traj/{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    as_path = f'{savepath}.astp' # (arm_stf_vtime_dic, ns_tm, ns_last, sys), TODO: reduce saved data
    as_results = save_load(as_path, None)
    if as_results == False:
        as_results = arm_stiffness_calc(top_path, traj_path, f'{savepath}.ns', f'{savepath}.sys', arm_num=arms) # ns: last_conf
        as_results = save_load(as_path, as_results)      
    as_dic, arms_idx = as_results
    # saveload end

    # data preparing
    stf_dic = OrderedDict() # {arm_id:{t:armstf}}
    for ia in arms_idx:
        stf_dic[ia] = OrderedDict()
    for t_stamp, as_ls in as_dic.items():
        for armstf, ia in as_ls:
            stf_dic[ia][t_stamp] = armstf
    stf_ls = [stf for s_a_dic in list(stf_dic.values()) for stf in list(s_a_dic.values())] # s_a_dic: {t:armstf} of a single arm; stf_ls:[stf]
    # t_ls = list(as_dic.keys()) # tracking one arm
    # data preparing ends
    # pre-plot
    n,bin_edges = np.histogram(stf_ls,bins=50, range = (0,1))
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
    # moments
    from scipy import stats
    n = np.array(n)
    m0 = np.sum(n) # 0th unitless raw moment: integration
    m1 = np.sum(stf_ls)/m0
    m2_c = stats.moment(stf_ls, moment=2) # 2nd central moment: variance
    std = m2_c**0.5 # standard deviation
    m3_c = stats.moment(stf_ls, moment=3) # 3rd central moment
    m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
    # moments ends
    # pre_plot ends
    # plot
    ax = plt.axes()
    n_sum = np.sum(n)
    plt.errorbar(bin_centers,np.divide(np.array(n),n_sum),yerr = np.divide(np.sqrt(np.array(n)),n_sum),marker = '.',drawstyle = 'steps-mid',c='#4994FF')
    plt.title(label, fontsize=18)
    # customization
    txt = (
        rf'$\mu={m1:.2f}$'
        '\n'
        rf'$\sigma={std:.2f}$'
        '\n'
        rf'$\gamma={m3_s:.2f}$'
    )
    plt.text(0.7, 0.215, txt, c='#4994FF', fontsize=18) # 28%
    plt.xlabel('Stiffness', fontsize=18)
    plt.ylabel('Frequencies', fontsize=18)
    plt.xlim((0,1))
    plt.ylim((0,0.3))
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_position((0.12501,0.12501,0.8078,0.8078))
    # customization ends
    chkdir(f'results/{arms}arms{conf_suffix}/')
    plt.savefig(f'results/{arms}arms{conf_suffix}/as_hist-{label}.png')
    # plt.show()
    plt.clf()
    # plot ends
    return m1,std,m3_s
