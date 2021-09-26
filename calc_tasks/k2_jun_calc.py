from readers import Reader, NanoConstructor
import numpy as np

def k2_jun_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
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


    k2_ls = []  # [(t_stamp, k2_val)]
    # for t_stamp, ns in ns_tm.time_capsule.items():
    for t_stamp, ns in ns_tm.time_capsule.items():
        i_arr33 = np.zeros((3,3))
        CoM_pos = np.zeros(3)
        base_cnt = 0
        center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
        center_ls.extend([base for arm in ns.arms.values() for base in list(arm.base_pairs.values())[0]])
        # CoM
        for base in center_ls:
            CoM_pos = np.add(CoM_pos, np.array(base.position))
            base_cnt += 1
        assert base_cnt == arm_num*(1+1+dims_ls[1])
        CoM_pos = np.divide(CoM_pos, base_cnt)
        for base in center_ls:
            x, y, z = base.position - CoM_pos
            i_arr33[0][0] += (y**2 + z**2)*m
            i_arr33[1][1] += (x**2 + z**2)*m
            i_arr33[2][2] += (y**2 + x**2)*m
            i_arr33[0][1] += -x*y*m
            i_arr33[1][0] += -x*y*m
            i_arr33[0][2] += -x*z*m
            i_arr33[2][0] += -x*z*m
            i_arr33[1][2] += -z*y*m
            i_arr33[2][1] += -z*y*m
        e_vals, _ = np.linalg.eig(i_arr33)
        l1, l2, l3 = e_vals
        k2 = 1 - (27*l1*l2*l3)/((l1+l2+l3)**3)
        k2_ls.append((t_stamp, k2))
    return k2_ls