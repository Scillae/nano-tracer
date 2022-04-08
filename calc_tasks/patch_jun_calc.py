from readers import Reader, NanoConstructor
import numpy as np
from collections import OrderedDict

def obtain_cos(v1,v2):
    return np.degrees(np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))) # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs

def patch_jun_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_dic = reader.read_data(p=sys_input)
    else:
        strands_dic = reader.read_data(obj=sys_input)
    # savepoint loading: nano-stars
    nc = NanoConstructor(strands_dic, dims_ls, arm_num)
    if type(ns_input) is str or ns_input == None: 
        ns_tm = nc.construct(p=ns_input)
    else:
        ns_tm = nc.construct(obj=ns_input)
    # finish savepoint loading
    
    # TODO: change the definition of CENTER to: CoM of all junction bases.
    p_angs_vtime_dic = OrderedDict() #{t_stamp: angle_results_ls}
    for t_stamp, ns in ns_tm.time_capsule.items():
        center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
        center_ls.extend([base for arm in ns.arms.values() for base in arm.base_pairs[1]]) # 
        junction_pos = np.average(np.array([base.position for base in center_ls]), 0) # assuming mass of nucs are the same.
        arms_idx = list(ns.arms.keys())
        is_sharing_strand = True
        angle_results_ls = []   # [(ang_rad, is_sharing)]
        for idx_1, ia1 in enumerate(arms_idx):
            arm1 = ns.arms[ia1]
            for idx_2 in range(idx_1+1, len(arms_idx)):
                ia2 = arms_idx[idx_2]
                arm2 = ns.arms[ia2]
                first_pair_a1 = arm1.base_pairs[1] # 
                first_pair_a2 = arm2.base_pairs[1]
                last_pair_a1 = arm1.base_pairs[dims_ls[0]] # tuples (s0[i],s1[i])
                last_pair_a2 = arm2.base_pairs[dims_ls[0]]
                base_ls = list(last_pair_a1)
                base_ls.extend(last_pair_a2)
                if len(set([base.strand_id for base in base_ls])) < 4:
                    is_sharing_strand = True
                else:
                    is_sharing_strand = False
                fp_a1_pos = np.average(np.array([base.position for base in first_pair_a1]),0)
                fp_a2_pos = np.average(np.array([base.position for base in first_pair_a2]),0)
                lp_a1_pos = np.average(np.array([base.position for base in last_pair_a1]),0) # assuming mass of nucs are the same.
                lp_a2_pos = np.average(np.array([base.position for base in last_pair_a2]),0)
                vec_1 = fp_a1_pos - junction_pos
                vec_2 = fp_a2_pos - junction_pos
                # vec_1 = lp_a1_pos - fp_a1_pos
                # vec_2 = lp_a2_pos - fp_a2_pos
                ang_cos = obtain_cos(vec_1, vec_2) # 0~180
                ang = ang_cos # if ang_cross >= 0 else (360 - ang_cos)
                angle_results_ls.append((ang, is_sharing_strand, (ia1, ia2)))
        print(angle_results_ls, len(angle_results_ls))
        p_angs_vtime_dic[t_stamp] = angle_results_ls
    return p_angs_vtime_dic, arms_idx # ns: last_conf

