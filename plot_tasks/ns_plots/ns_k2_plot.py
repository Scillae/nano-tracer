from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def ns_k2_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    if single == True:
        top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
        traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
    else:
        if conf_suffix == '-x20':
            top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
            traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
        else:
            assert 0 == 1
    savepath = f'data/composed_traj/{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    
    if single == True:
        # #### plot ####
        k2_path = f'{savepath}.k2tp' # (k2_ls)
        k2_ls_results = save_load(k2_path, None)
        if k2_ls_results == False:
            k2_ls_results = k2_calc(top_path, traj_path, f'{savepath}.sys')
            k2_ls_results = save_load(k2_path, k2_ls_results)
        k2_tp_ls = k2_ls_results

        #strip off time
        k2_ls = [i[1] for i in k2_tp_ls]

        n,bin_edges = np.histogram(k2_ls,bins = 50, range = (0,1))
        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

        # moments
        from scipy import stats
        n = np.array(n)
        m0 = np.sum(n) # 0th unitless raw moment: integration
        m1 = np.sum(k2_ls)/m0
        m2_c = stats.moment(k2_ls, moment=2) # 2nd central moment: variance
        std = m2_c**0.5 # standard deviation
        m3_c = stats.moment(k2_ls, moment=3) # 3rd central moment
        m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
        # moments ends
        # plt.figure(figsize=(3,3))
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
        plt.xlabel('k2 value', fontsize=18)
        plt.ylabel('Frequencies', fontsize=18)
        plt.xlim((0,1))
        plt.ylim((0,0.3))
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # ax.set_xticklabels(fontsize=14) # [int(i) for i in np.linspace(0,1,num=10)],
        # ax.set_yticklabels(fontsize=14) # np.linspace(0,0.16,num=9),
        ax.set_position((0.12501,0.12501,0.8078,0.8078))
        # customization ends
        chkdir(f'results/{arms}arms{conf_suffix}/')
        plt.savefig(f'results/{arms}arms{conf_suffix}/k2_hist-{label}.png')
        # plt.show()
        plt.clf()
        return m1,std,m3_s
    assert 0 == 1