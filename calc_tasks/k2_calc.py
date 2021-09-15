from readers import Reader, NanoConstructor
import numpy as np

def k2_calc(path_top = 'data/6-arm-nanostar-starlike-412kInit.top', path_traj = 'data/6-arm-nanostar-starlike-412kInit.conf', sys_input = None, dims_ls = [20,2,7]):
    m = 1 # mass of a single nucleotide
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_tm = reader.read_data(p=sys_input)
    else:
        strands_tm = reader.read_data(obj=sys_input)

    # # Noella: Accepted_sci
    # nc = NanoConstructor(strands_tm, dims_ls)
    # ns_tm = nc.construct()
    k2_ls = []  # [(t_stamp, k2_val)]
    # for t_stamp, ns in ns_tm.time_capsule.items():
    for t_stamp, strands in strands_tm.time_capsule.items():
        i_arr33 = np.zeros((3,3))
        CoM_pos = np.zeros(3)
        base_cnt = 0
        # CoM
        for strand in strands.values():
            for base in strand.base_sequence.values():
                CoM_pos = np.add(CoM_pos, np.array(base.position))
                base_cnt += 1
        CoM_pos = np.divide(CoM_pos, base_cnt)
        for strand in strands.values():
            for base in strand.base_sequence.values():
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
        e_vals, e_vecs = np.linalg.eig(i_arr33)
        l1, l2, l3 = e_vals
        k2 = 1 - (27*l1*l2*l3)/((l1+l2+l3)**3)
        k2_ls.append((t_stamp, k2))
    return k2_ls
