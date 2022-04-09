from readers import Reader, NanoConstructor
import numpy as np
from collections import OrderedDict
# Noella: Accepted_sci
# call this function in main.py
# if __name__ == '__main__':

def obtain_cos(v1,v2):
    return np.degrees(np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))) # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs

def obtain_cross(v1,v2):
    return np.cross(v1,v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))) # sin(x) = (n1 x n2) / (|n1|*|n2|), angle only -90~90

def patch_angle_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
    deprecation = False # using the CoM_aligned method would make the length -1.
    if deprecation:
        return patch_angle_calc_CoM_aligned(path_top, path_traj, arm_num, dims_ls, ns_input, sys_input)
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_dic = reader.read_data(p=sys_input)
    else:
        strands_dic = reader.read_data(obj=sys_input)
    # savepoint loading: nano-stars
    nc = NanoConstructor(strands_dic, dims_ls, arm_num)
    if type(ns_input) is str or ns_input == None: 
        # box_dim hacking
        import re
        with open(path_traj,'r') as f:
            f.readline()
            ret=re.match('^b = ([0-9]+) ([0-9]+) ([0-9]+)\n',r)
        box_dim = np.array((ret.group(1),ret.group(2),ret.group(3)))
        ns_tm = nc.construct(p=ns_input, box_dim=box_dim)
    else:
        ns_tm = nc.construct(obj=ns_input)
    # finish savepoint loading

    p_angs_vtime_dic = OrderedDict() #{t_stamp: angle_results_ls}
    for t_stamp, ns in ns_tm.time_capsule.items():
        # Change of definition: origin~arm_end --> arm_start~arm_end
        # junction_pos = np.average(np.array([base.position for base in ns.center.values()]), 0) # assuming mass of nucs are the same.
        arms = ns.arms
        arms_idx = list(arms.keys())
        is_sharing_strand = True
        angle_results_ls = []   # [(ang_rad, is_sharing)]
        for idx_1, ia1 in enumerate(arms_idx):
            arm1 = arms[ia1]
            for idx_2 in range(idx_1+1, len(arms_idx)):
                ia2 = arms_idx[idx_2]
                arm2 = arms[ia2]
                first_pair_a1 = arm1.base_pairs[dims_ls[0]-2-10]
                first_pair_a2 = arm2.base_pairs[dims_ls[0]-2-10]
                last_pair_a1 = arm1.base_pairs[dims_ls[0]-2] # tuples (s0[i],s1[i])
                last_pair_a2 = arm2.base_pairs[dims_ls[0]-2]
                if dims_ls[0] < 13:
                    first_pair_a1 = arm1.base_pairs[dims_ls[0]-10] # b-form DNA, loop length == 10
                    first_pair_a2 = arm2.base_pairs[dims_ls[0]-10]
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
                # vec_1 = lp_a1_pos - junction_pos
                # vec_2 = lp_a2_pos - junction_pos
                vec_1 = lp_a1_pos - fp_a1_pos
                vec_2 = lp_a2_pos - fp_a2_pos
                ang_cos = obtain_cos(vec_1, vec_2) # 0~180
                # ang_cross = obtain_cross(vec_1, vec_2) # -90~90
                ang = ang_cos # if ang_cross >= 0 else (360 - ang_cos)
                angle_results_ls.append((ang, is_sharing_strand, (ia1, ia2), (vec_1, vec_2)))
        print(angle_results_ls, len(angle_results_ls))
        p_angs_vtime_dic[t_stamp] = angle_results_ls
    return p_angs_vtime_dic, arms_idx # ns: last_conf


def patch_angle_calc_CoM_aligned(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
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

    p_angs_vtime_dic = OrderedDict() #{t_stamp: angle_results_ls}
    z_arm_id = 0
    arm_normalize = True
    t_series = list(ns_tm.time_capsule.keys())
    t_series.sort()
    t_step = t_series[1] - t_series[0] # assume time series is continuous
    t_max = max(t_series)
    vec_dic = OrderedDict()
    vec_jun_dic = OrderedDict()
    arm_shift_log_dic = OrderedDict()
    arms_idx = list(list(ns_tm.time_capsule.values())[0].arms.keys())
    for t_stamp, ns in ns_tm.time_capsule.items():
        if t_stamp == t_max: # if-block: read next timestep
            vec_dic[t_stamp], vec_jun_dic[t_stamp], arm_shift_log_dic[t_stamp] = nanostar_vectorize(ns, dims_ls ,z_arm_id, arm_normalize, ns_tm.time_capsule[t_stamp])
        else:
            t_stamp_next = t_stamp + t_step
            if t_stamp_next not in ns_tm.time_capsule.keys():
                continue
            vec_dic[t_stamp], vec_jun_dic[t_stamp], arm_shift_log_dic[t_stamp] = nanostar_vectorize(ns, dims_ls ,z_arm_id, arm_normalize, ns_tm.time_capsule[t_stamp_next])

    # shielding...
    for t_stamp, vec_arms_ls in vec_dic.items():
        angle_results_ls = [] # [(ang, is_sharing_strand, (ia1, ia2))]
        if type(vec_arms_ls) == bool:
            continue
        vec_jun_ls = vec_jun_dic[t_stamp]
        for i, (ia2,(x,y,z),(shareia_2,shareia2_2), _) in enumerate(vec_arms_ls):
            if i == 0:
                continue
            _,(xj,yj,zj),_,_ = vec_jun_ls[i]
            vec_2 = np.array((x,y,z)) - np.array((xj,yj,zj))
            for j in range(i):
                ia1,(x,y,z),(shareia,shareia2),_ = vec_arms_ls[j]
                _,(xj,yj,zj),_,_ = vec_jun_ls[j]
                vec_1 = np.array((x,y,z)) - np.array((xj,yj,zj))
                ang_cos = obtain_cos(vec_1, vec_2) # 0~180
                # cross = obtain_cross(vec_1, vec_2) # cross product
                ang = ang_cos # if cross[2] >= 0 else (360 - ang_cos)
                is_sharing_strand = True if shareia_2 in [shareia,shareia2] or shareia2_2 in [shareia,shareia2] else False
                angle_results_ls.append((ang, is_sharing_strand, (ia1, ia2)))
        p_angs_vtime_dic[t_stamp] = angle_results_ls
    return p_angs_vtime_dic, arms_idx # ns: last_conf


#### Below: From report_plot.py ####

def CoM_calc(ns):
    # from jun_shift_calc.py
    CoM_pos = np.zeros(3)
    CoM_cen_pos = np.zeros(3)
    base_cnt = 0
    if ns.center == None:
        center_ls = []
    else:
        center_ls = list(ns.center.values()) # center_ls: [center_base], len == 4*arm_num
    center_ls.extend([base for arm in ns.arms.values() for base in arm.base_pairs[1]]) # list(arm.base_pairs.values())[-1]] # if they are basepairs or bases
    # CoM
    for base in center_ls:
        CoM_cen_pos = np.add(CoM_cen_pos, np.array(base.position))
        base_cnt += 1
    CoM_cen_pos = np.divide(CoM_cen_pos, base_cnt)
    base_cnt = 0
    for strand in ns.strands.values():
        for base in strand.base_sequence.values():
            CoM_pos = np.add(CoM_pos, np.array(base.position))
            base_cnt += 1
    CoM_pos = np.divide(CoM_pos, base_cnt)
    return CoM_pos, CoM_cen_pos

def coord_rotate(ns, dims_ls, z_arm_id, next_ns = None):
    '''
    Create a rotation matrix that transform the space so that: 
        the selected arm aligns with z-axis; 
        cross the 1st z-arm base's backbone-to-base versor with z-axis to obtain y axis.
    :ns: nanostar object
    :z_arm_id: the arm chosen to be z-axis.
    '''
    CoM_pos, CoM_cen_pos = CoM_calc(ns)
    next_CoM_pos, next_CoM_cen_pos = CoM_calc(next_ns)
    z_pair = [base for base in ns.arms[z_arm_id].base_pairs[1]] # start
    z_end_pair = [base for base in ns.arms[z_arm_id].base_pairs[dims_ls[0]]] # end
    z_arm = ns.arms[z_arm_id]
    z = (np.array(z_end_pair[1].position) - CoM_pos + np.array(z_end_pair[0].position) - CoM_pos)/2 - (np.array(z_pair[1].position) - CoM_pos + np.array(z_pair[0].position) - CoM_pos)/2 # z is Translation invariant though...
    z_unit = z/np.linalg.norm(z)
    y = np.cross(z_unit, np.array(z_pair[0].backbone))
    if np.dot(y,y) < 0.9 : print(f'Magnitude of Y: {np.dot(y,y)} < 0.90')
    y_unit = y/np.linalg.norm(y)
    x_unit = np.cross(y_unit,z_unit)
    if np.dot(x_unit,x_unit) < 0.9 : print(f'Magnitude of x-unit: {np.dot(x_unit,x_unit)} < 0.90')
    rot_mat = np.vstack((x_unit,y_unit,z_unit)).T
    rot_mat = np.linalg.inv(rot_mat)
    # arm displacement tracking
    arm_shift_log = OrderedDict() # shifts of arm-ending positions in **the moving CoM space**
    for ia, arm in ns.arms.items():
        arm_next = next_ns.arms[ia]
        current_pair = arm.base_pairs[dims_ls[0]]
        next_pair = arm_next.base_pairs[dims_ls[0]]
        current_pos = (np.array(current_pair[0].position) + np.array(current_pair[1].position))/2 - CoM_pos
        next_pos = (np.array(next_pair[0].position) + np.array(next_pair[1].position))/2 - next_CoM_pos
        arm_shift_log[ia] = next_pos - current_pos
        # print(f'For Patch Angle vTime use. Arm_{arm.arm_id} ~ strands: {arm.strand_id_0} & {arm.strand_id_1}')
    arm_shift_log['CoM'] = next_CoM_pos - CoM_pos
    # shielding last time_stamp
    if next_ns == ns or next_ns == None:
        return False, False   
    return rot_mat, arm_shift_log

def convert_spherical_to_rectangular(r,t,p):
    x = r*np.sin(t)*np.cos(p)
    y = r*np.sin(t)*np.sin(p)
    z = r*np.cos(t)
    return np.array((x,y,z))

def append_spherical_np(xyz):
    rtp = np.zeros(xyz.shape)
    xy = xyz[0,:]**2 + xyz[1,:]**2
    rtp[0,:] = np.sqrt(xy + xyz[2,:]**2)
    rtp[1,:] = np.arctan2(np.sqrt(xy), xyz[2,:]) # for elevation angle defined from Z-axis down
    rtp[2,:] = np.arctan2(xyz[1,:], xyz[0,:])
    return np.vstack((xyz,rtp))

def nanostar_vectorize(ns, dims_ls, sel_arm_id, arm_normalize, next_ns = None):
    '''
    Transform a nanostar into two sets of points. 
        CoM (of whole nanostar) ~ Origin, CoM to last pair of an arm ~ arm vector (vec_ns_ls), CoM to first pair of an arm ~ junction vector (vec_jun_ls)
    Vectors are NOT normalized yet.
    :sel_arm_id: the arm chosen to be z-axis.
    '''
    vec_ns_ls = []
    vec_jun_ls = []
    CoM, CoM_cen = CoM_calc(ns) # now arms are vectors pointing from CoM instead of 1st pair.
    rot_mat, arm_shift_log = coord_rotate(ns, dims_ls, sel_arm_id, next_ns) # _CoMs_align
    # shielding...
    if type(rot_mat) == bool: # bool: False
        return (False,False,False)
    for ia, arm in ns.arms.items():
        end_base_pair = [base for base in arm.base_pairs[dims_ls[0]]] # 
        start_base_pair = [base for base in arm.base_pairs[1]] # 
        end_pos = (np.array(end_base_pair[1].position) - CoM + np.array(end_base_pair[0].position) - CoM)/2 # now in CoM coord TODO: confirm
        start_pos = (np.array(start_base_pair[1].position) - CoM + np.array(start_base_pair[0].position) - CoM)/2 # now in CoM coord
        vec = end_pos
        vec_jun = start_pos
        vec = np.matmul(rot_mat,vec)
        vec_jun = np.matmul(rot_mat,vec_jun)
        if arm_normalize: 
            arm_len = np.sqrt(np.dot(vec,vec))
            vec = vec/arm_len
            vec_jun = vec_jun/arm_len
        strands_using = (arm.strand_id_0, arm.strand_id_1)
        is_sharing_sel = ns.arms[sel_arm_id].strand_id_0 in strands_using or ns.arms[sel_arm_id].strand_id_1 in strands_using
        if is_sharing_sel:
            is_sharing_sel = 1 if ns.arms[sel_arm_id].strand_id_0 in strands_using else 2
        # print(f'Arm ID: {ia} ; is_sharing: {is_sharing_sel}')
        vec_ns_ls.append((ia, tuple(vec), strands_using, is_sharing_sel))
        vec_jun_ls.append((ia, tuple(vec_jun), strands_using, is_sharing_sel))
    return vec_ns_ls, vec_jun_ls, arm_shift_log
