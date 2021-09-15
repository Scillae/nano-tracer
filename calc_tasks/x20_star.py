from readers import Reader, NanoConstructor
import numpy as np

def x20_star(path_top = 'data/6-arm-nanostar-starlike-412kInit.top', path_traj = 'data/6-arm-nanostar-starlike-412kInit.conf', nm_input = None, sys_input = None, arm_n = 4):
    dims_ls = [20,2,7]
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_dic = reader.read_data(p=sys_input)
    else:
        strands_dic = reader.read_data(obj=sys_input)
    # savepoint loading: nano-stars
    nc = NanoConstructor(strands_dic, dims_ls, arm_n)
    if type(nm_input) is str or nm_input == None: 
        nm_tm = nc.construct(p=nm_input)
    else:
        nm_tm = nc.construct(obj=nm_input)
    # finish savepoint loading

    # patch angle
    p_angs_vtime_dic = {}
    for t_stamp, nm in nm_tm.time_capsule.items():
        p_angs_vnsid_dic = {}
        for ns_idx, (ns, arm_num) in nm.nodes.items():
            # Noella: Questioned_sci(dic or list) - dict
            # base is a tuple defined at line 21 in NanoConstructor
            # junction_pos = np.average(np.array([base.position for base in ns.center.values()]),-1)
            junction_pos = np.average(np.array([base.position for base in ns.center.values()]), 0) # assuming mass of nucs are the same.
            arms = ns.arms
            arms_idx = list(arms.keys())
            is_sharing_strand = True
            angle_results_ls = []   # [(ang_rad, is_sharing)]
            for idx_1, ia1 in enumerate(arms_idx):
                arm1 = arms[ia1]
                for idx_2 in range(idx_1, len(arms_idx)):
                    ia2 = arms_idx[idx_2]
                    if ia2 == ia1:
                        continue
                    arm2 = arms[ia2]
                    last_pair_a1 = list(arm1.base_pairs.values())[-1] # tuples (s0[i],s1[i])
                    last_pair_a2 = list(arm2.base_pairs.values())[-1]
                    base_ls = list(last_pair_a1)
                    base_ls.extend(last_pair_a2)
                    # Noella: Questioned_sci(dic or list); Other implementation is wrong: set()
                    # 1. other implementation for following 4 lines:
                    # is_sharing_strand = (len(base_ls) < 4)
                    # 2. base is a tuple
                    # if len(set([base.strand_id for base in base_ls])) < 4:
                    if len(set([base.strand_id for base in base_ls])) < 4:
                        is_sharing_strand = True
                    else:
                        is_sharing_strand = False
                    lp_a1_pos = np.average(np.array([base.position for base in last_pair_a1]),0) # assuming mass of nucs are the same.
                    lp_a2_pos = np.average(np.array([base.position for base in last_pair_a2]),0)
                    vec_1 = lp_a1_pos - junction_pos
                    vec_2 = lp_a2_pos - junction_pos
                    ang_rad = np.arccos(np.sum(vec_1*vec_2)/(np.sqrt(np.sum(np.square(vec_1))) * np.sqrt(np.sum(np.square(vec_2))))) # cos(x) = n1*n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs
                    angle_results_ls.append((ang_rad, is_sharing_strand))
            print(angle_results_ls, len(angle_results_ls))
            p_angs_vnsid_dic[ns_idx] = angle_results_ls
        p_angs_vtime_dic[t_stamp] = p_angs_vnsid_dic # {t:{nsid:[angles]}}
    return p_angs_vtime_dic, nm_tm, strands_dic