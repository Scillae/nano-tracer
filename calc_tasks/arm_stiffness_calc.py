from readers import Reader, NanoConstructor
import numpy as np
# Noella: Accepted_sci
# call this function in main.py
# if __name__ == '__main__':
def arm_stiffness_calc(path_top = 'data/6-arm-nanostar-starlike-412kInit.top', path_traj = 'data/6-arm-nanostar-starlike-412kInit.conf', ns_input = None, sys_input = None, arm_num = 6, dims_ls= [20,2,7], v_num = 2):
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

    arm_stf_vtime_dic = {} # {t_stamp: armstf_ls}
    for t_stamp, ns in ns_tm.time_capsule.items():
        arms = ns.arms
        armstf_ls = []
        v_len = dims_ls[0]//v_num # v_num: number of vectors used to evaluate the stiffness
        for arm in arms.values():
            v_ls = []
            for i in range(v_num):
                st_pair = list(arm.base_pairs.values())[v_len*i]
                ed_pair = list(arm.base_pairs.values())[v_len*(i+1)-1]
                st_pos = np.average(np.array([base.position for base in st_pair]),0)
                ed_pos = np.average(np.array([base.position for base in ed_pair]),0)
                vec = st_pos-ed_pos
                v_normed = vec / np.linalg.norm(vec)
                v_ls.append(v_normed)
            armstf = 0
            for i in range(len(v_ls)-2+1): # +1 inclusive
                v1 = v_ls[i]
                v2 = v_ls[i+1]
                ang_rad = np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2))))) # cos(x) = n1*n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs
                armstf += ang_rad / np.pi # armstf: add 0~1 per time
            armstf /= v_num-1
            armstf_ls.append((armstf, arm.arm_id)) # len = arm_num
        arm_stf_vtime_dic[t_stamp] = armstf_ls
    return arm_stf_vtime_dic, list(arms.keys())
