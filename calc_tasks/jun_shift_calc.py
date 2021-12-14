from readers import Reader, NanoConstructor
import numpy as np

def jun_shift_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
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


    js_ls = []  # [(t_stamp, k2_val)]
    # for t_stamp, ns in ns_tm.time_capsule.items():
    for t_stamp, ns in ns_tm.time_capsule.items():
        jun_shift = 0
        CoM_pos = np.zeros(3)
        CoM_cen_pos = np.zeros(3)
        base_cnt = 0
        if ns.center == None:
            center_ls = []
        else:
            center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
        center_ls.extend([base for arm in ns.arms.values() for base in arm.base_pairs[1]]) # list(arm.base_pairs.values())[-1]]
        # CoM
        for base in center_ls:
            CoM_cen_pos = np.add(CoM_cen_pos, np.array(base.position))
            base_cnt += 1
        assert base_cnt == arm_num*(1+1+dims_ls[1])
        CoM_cen_pos = np.divide(CoM_cen_pos, base_cnt)
        base_cnt = 0
        for strand in ns.strands.values():
            for base in strand.base_sequence.values():
                CoM_pos = np.add(CoM_pos, np.array(base.position))
                base_cnt += 1
        assert base_cnt == arm_num*(dims_ls[0]+dims_ls[0]+dims_ls[1]+dims_ls[2])
        CoM_pos = np.divide(CoM_pos, base_cnt)
        jun_shift = CoM_pos - CoM_cen_pos
        jun_shift = np.sqrt(np.dot(jun_shift, jun_shift))
        js_ls.append((t_stamp, jun_shift))
    return js_ls