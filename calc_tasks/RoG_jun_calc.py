from readers import Reader, NanoConstructor
import numpy as np

def RoG_jun_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
    m = 1 # mass of a single nucleotide
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


    RoG_ls = []  # [(t_stamp, k2_val)]
    # for t_stamp, ns in ns_tm.time_capsule.items():
    for t_stamp, ns in ns_tm.time_capsule.items():
        RoG = 0
        CoM_pos = np.zeros(3)
        base_cnt = 0
        center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
        center_ls.extend([base for arm in ns.arms.values() for base in list(arm.base_pairs.values())[-1]])
        # CoM
        for base in center_ls:
            CoM_pos = np.add(CoM_pos, np.array(base.position))
            base_cnt += 1
        assert base_cnt == arm_num*(1+1+dims_ls[1])
        CoM_pos = np.divide(CoM_pos, base_cnt)
        for base in center_ls:
            v = base.position - CoM_pos
            RoG += np.dot(v,v)
        RoG = np.sqrt(RoG/base_cnt)
        RoG_ls.append((t_stamp, RoG))
    return RoG_ls