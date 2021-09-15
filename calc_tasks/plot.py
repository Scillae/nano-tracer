from readers import Reader, NanoConstructor
from calc_tasks import patch_angle_calc, k2_calc, x20_star
from collections import OrderedDict
import pickle
import os.path
import numpy as np

# 08/21-temporary
def save_load(p, obj):
    print(f'TM path: {p}')
    if p is None:
        print('TM skipping!')
        return obj
    if (obj is not None) and (not os.path.isfile(p)):
        if os.path.isdir(os.path.split(p)[0]) == False:
            os.makedirs(os.path.split(p)[0])
        print('TM saving!')
        pickle.dump(obj, open(p,"wb"))
        r_obj = obj
    elif (obj is None) and os.path.isfile(p):
        print('TM loading!')
        r_obj = pickle.load(open(p, "rb"))
    elif (obj is not None) and os.path.isfile(p):
        print('TM updating savepoint!')
        pickle.dump(obj, open(p,"wb"))
        r_obj = obj
    else:
        print('TM save_load both empty')
        r_obj = False
    return r_obj
# temporary ends.


# Noella: Accepted_sci
# call patch_angle_calc in main.py
if __name__ == '__main__':
    # switches here
    single = True
    arms = 4
    temp = 30
    conc = 0.5
    sp_suffix = '' # '-redo'
    conf_suffix = '' # '-x20'
    # switches above

    # # switches here
    # single = False
    # arms = 4
    # temp = 30
    # conc = 0.1
    # sp_suffix = '' # '-redo'
    # conf_suffix = '-x20'
    # # switches above

    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'

    if single == True:
        top_path = f'../../ox-sync/simul-inputs-{arms}arms/trial_{arms}arm_rod.top'
        traj_path = f'../../ox-sync/simul-inputs-{arms}arms/{loose_lbl}/trajectory.dat'
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
            p_angs_results = patch_angle_calc(top_path, traj_path, f'{savepath}.ns', f'{savepath}.sys', arm_num=arms) # ns: last_conf
            p_angs_results = save_load(pa_path, p_angs_results)      
        p_angs_dic, ns_tm, ns_last, sys = p_angs_results
        #### plot import ####
        import matplotlib.pyplot as plt
        import numpy as np
        # #### plot ####
        # tracing of specific p-angles setup
        arms_idx = list(list(ns_tm.time_capsule.values())[0].arms.keys())
        angle_dic = OrderedDict() # {(ia1, ia2):{t_stamp:angle}} arms_idx = list(list(ns_tm.time_capsule.values())[0].arms.keys())
        for idx_1, ia1 in enumerate(arms_idx):
            for idx_2 in range(idx_1+1, len(arms_idx)):
                ia2 = arms_idx[idx_2]
                angle_dic[(ia1, ia2)] = OrderedDict()
        # hist of p-angle
        angle_ls = [] # historical approach.
        t_ls = []        
        for t_stamp, ang_ls in p_angs_dic.items():
            ## sanity chk: because of old pairing method.
            s = [1 for i in ang_ls if i[1] == True]
            if len(s) != len(ns_last.arms): # n arms should yield n angles.
                print(t_stamp)
                continue
            ## chk ends
            for ang, is_sharing, ia_tp in ang_ls:
                if is_sharing == True:  # patch angle: sharing a strand
                    angle_ls.append(ang)
                    t_ls.append(t_stamp)
                    # {(ia1, ia2):{t_stamp:angle}}
                    angle_dic[ia_tp][t_stamp] = ang # tracing of specific p-angles
        ang_arr = np.array(angle_ls)
        ang_deg = np.degrees(ang_arr)
        # n, bins, patches = plt.hist(ang_deg,bins = 180, range = (0,180))
        # n, bins, patches = plt.hist(ang_deg,bins = 36, range = (0,180))
        n,bin_edges = np.histogram(ang_deg,bins=36, range = (0,180))
        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
        # moments
        from scipy import stats
        n = np.array(n)
        # m0 = np.sum(n*(np.array(bins[0:len(bins)-1])+0.5)) # 0th raw moment: integration
        # m0 = np.sum(n*(np.array(bin_edges[0:len(bin_edges)-1])+0.5)) # 0th raw moment: integration
        m0 = np.sum(n) # 0th unitless raw moment: integration
        # m1 = m0/np.sum(n) # 1st raw moment: mean
        m1 = np.sum(ang_deg)/m0
        # n = n/m0 # normalization, integration m0 == 1 now
        # m2_c = stats.moment(n, moment=2) # 2nd central moment: variance
        m2_c = stats.moment(ang_deg, moment=2) # 2nd central moment: variance
        std = m2_c**0.5 # standard deviation
        # m3_c = stats.moment(n, moment=3) # 3rd central moment
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
        # plt.title(f'Histogram of patch angles, {label}')
        plt.title(label, fontsize=18)
        txt = '\n'.join((
            r'$\mu=%d$' % (m1, ),
            r'$\sigma=%d$' % (std, ),
            r'$\gamma=%.1f$' % (m3_s, ))) # was  rf'N={m0}, $\mu$={m1}, $\sigma$={std}, $\gamma$={m3_s}'
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
        plt.savefig(f'results/{arms}arms{conf_suffix}/pa_hist-{label}.png')
        plt.show()
        # # line of p-angle
        # plt.scatter(t_ls,ang_deg)
        # plt.title(f'Patch angles vs time, {label}')
        # plt.savefig(f'results/{arms}arms{conf_suffix}/pa_v_t-{label}.png')
        # plt.show()
        # print(len(ang_arr))
        # #### plot end ####
        # k2_ls = k2_calc(top_path, traj_path, sys)f'{savepath}.sys'
        # k2_path = f'{savepath}.k2tp' # (k2_ls)
        # k2_ls_results = save_load(k2_path, None)
        # if k2_ls_results == False:
        #     k2_ls_results = k2_calc(top_path, traj_path, f'{savepath}.sys')
        #     k2_ls_results = save_load(k2_path, k2_ls_results)
        # k2_ls = k2_ls_results
        # #### plot ####
        # plt.clf()
        # # line of k2        
        # t = [k2_tp[0] for k2_tp in k2_ls]
        # k2 = [k2_tp[1] for k2_tp in k2_ls]
        # # plt.plot(t,k2)
        # # plt.title(f'k2 vs time, {label}')
        # #plt.savefig(f'results/{arms}arms{conf_suffix}/k2_v_t-{label}.png')
        # # plt.show()
        # # hist of k2
        # # n, bins, patches = plt.hist(k2,bins = 100, range = (0,1))
        # n, bins, patches = plt.hist(k2,bins = 50, range = (0,1))
        # # # moments
        # # from scipy import stats
        # # n = np.array(n)
        # # m0 = np.sum(n*(np.array(bins[0:len(bins)-1])+0.5)) # 0th raw moment: integration
        # # m1 = m0/np.sum(n) # 1st raw moment: mean
        # # n = n/m0 # normalization, integration m0 == 1 now
        # # m2_c = stats.moment(n, moment=2) # 2nd central moment: variance
        # # std = m2_c**0.5 # standard deviation
        # # m3_c = stats.moment(n, moment=3) # 3rd central moment
        # # m3_s = m3_c / (std**3) # 3rd standardized moment: skewness
        # # # moments ends
        # plt.title(f'Histogram of k2, {label}')
        # plt.savefig(f'results/{arms}arms{conf_suffix}/k2_hist-{label}.png')
        # plt.show()
        ##### plot end ####
    else:
        p_angs_dic, nm_tm, sys = x20_star(top_path, traj_path, f'{savepath}.ns', f'{savepath}.sys', arm_n=arms)
        #### plot import ####
        import matplotlib.pyplot as plt
        import numpy as np
        #### plot ####       
        # hist of p-angle
        angle_ls = []
        t_ls = []
        for t_stamp, ang_dic in p_angs_dic.items():
            for nsid, ang_ls in ang_dic.items():
                for ang, is_sharing in ang_ls:
                    if is_sharing == True: # patch angle: sharing a strand
                        angle_ls.append(ang)
                        t_ls.append(t_stamp)
        ang_deg = np.degrees(np.array(angle_ls))
        plt.hist(ang_deg,bins = 180, range = (0,180))
        plt.title(f'Histogram of patch angles, {label}')
        plt.savefig(f'results/{arms}arms{conf_suffix}/pa_hist-{label}.png')
        plt.show()
        # line of p-angle
        plt.scatter(t_ls,ang_deg)
        plt.title(f'Patch angles vs time, {label}')
        #plt.savefig(f'results/{arms}arms{conf_suffix}/pa_v_t-{label}.png')
        plt.show()
        print(len(angle_ls))
        #### plot end ####