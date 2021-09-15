from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star, arm_stiffness_calc
from utils.tools import save_load, chkdir
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def ns_pa_plot(single=True, arms=4, temp=30, conc=0.5, sp_suffix='', conf_suffix='', dims_ls= [20,2,7]):
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
        pa_path = f'{savepath}.patp' # (p_angs_dic, ns_tm, ns_last, sys)
        p_angs_results = save_load(pa_path, None)
        if p_angs_results == False:
            p_angs_results = patch_angle_calc(top_path, traj_path, f'{savepath}.ns', f'{savepath}.sys', arm_num=arms)
            p_angs_results = save_load(pa_path, p_angs_results)      
        p_angs_dic, arms_idx = p_angs_results
        # #### plot ####
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
        n,bin_edges = np.histogram(ang_deg,bins=36, range = (0,180))
        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
        # moments
        from scipy import stats
        n = np.array(n)
        m0 = np.sum(n) # 0th unitless raw moment: integration
        m1 = np.sum(ang_deg)/m0
        m2_c = stats.moment(ang_deg, moment=2) # 2nd central moment: variance
        std = m2_c**0.5 # standard deviation
        m3_c = stats.moment(ang_deg, moment=3) # 3rd central moment
        m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
        # moments ends
        # # TODO: plot single_batch angle here:
        # '''
        # for ia_tp in angle_dic.keys():
        #     t_ls = list(angle_dic[ia_tp].keys())
        #     a_ls = list(angle_dic[ia_tp].values())
        #     print(f'Angle:{ia_tp}, time: {t_ls} vs angle: {a_ls}')
        # '''
        n_sum = np.sum(n)
        ax = plt.axes()
        plt.errorbar(bin_centers,np.divide(np.array(n),n_sum),yerr = np.divide(np.sqrt(np.array(n)),n_sum),marker = '.',drawstyle = 'steps-mid')
        plt.title(label, fontsize=18)
        txt = (
            rf'$\mu={m1:.2f}$'
            '\n'
            rf'$\sigma={std:.2f}$'
            '\n'
            rf'$\gamma={m3_s:.2f}$'
        )
        plt.text(10, 0.12, txt, c='b', fontsize=18)
        plt.xlabel('Patch angle (degrees)', fontsize=18)
        plt.ylabel('Frequencies', fontsize=18)
        plt.xlim((0,180))
        plt.ylim((0,0.17))
        plt.xticks(np.linspace(0,180,num=10))
        plt.yticks(np.linspace(0,0.16,num=9))
        x_lb_tmp = [int(i) for i in np.linspace(0,180,num=10)] #temporary
        ax.set_xticklabels(x_lb_tmp,fontsize=14)
        ax.set_yticklabels(np.linspace(0,0.16,num=9),fontsize=14)
        ax.set_position((0.12501,0.12501,0.8078,0.8078))
        chkdir(f'results/{arms}arms{conf_suffix}/')
        plt.savefig(f'results/{arms}arms{conf_suffix}/pa_hist-{label}.png')
        plt.clf()
        return m1, std, m3_s
    assert 0 == 1
