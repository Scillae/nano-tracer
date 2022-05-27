from readers import Reader, NanoConstructor
import numpy as np
from collections import OrderedDict

def stacking_local_identify_calc(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_tm = reader.read_data(p=sys_input)
    else:
        strands_tm = reader.read_data(obj=sys_input)
    # savepoint loading: nano-stars
    nc = NanoConstructor(strands_tm, dims_ls, arm_num)
    if type(ns_input) is str or ns_input == None: 
        # box_dim hacking
        import re
        with open(path_traj,'r') as f:
            f.readline()
            ret=re.match('^b = ([0-9]+) ([0-9]+) ([0-9]+)\n',f.readline())
        box_dim = np.array((int(ret.group(1)),int(ret.group(2)),int(ret.group(3))))
        ns_tm = nc.construct(p=ns_input, box_dim=box_dim)
    else:
        ns_tm = nc.construct(obj=ns_input)
    # finish savepoint loading

    from copy import deepcopy
    stacking_dict = OrderedDict() # OrderedDict{t_stamp: is_stacking, stack_ls}
    for t_stamp, ns in ns_tm.time_capsule.items():
        is_stacking = False
        stack_ls = [] # [(adj_bp, adj_bp2, arm_id, arm2_id)]
        # generate potential stk_arm_pairs
        ptt_aps = [] # [(arm1,arm2)]
        arms_idx = list(ns.arms.keys())
        for idx_0, ia0 in enumerate(arms_idx):
            arm0 = ns.arms[ia0]
            arm0_dir = get_arm_dir(arm0)
            for idx_1 in range(idx_0+1, len(arms_idx)):
                arm1 = ns.arms[arms_idx[idx_1]]
                # find if angle~ {arm_dir0:arm_dir1} < threshold
                arm1_dir = get_arm_dir(arm1)
                if obtain_cos(arm0_dir, arm1_dir*(-1)) < 30: # one arm reversed. 30deg --> ang==150deg
                    ptt_aps.append((arm0, arm1))
        # check if no ptt_aps
        if not ptt_aps: # bool(empty_dict) == false
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
        # identify & create jxn_bps
        jxn_bps = [(b0,b1) for b0 in ns.center.values() for b1 in ns.center.values() if check_if_bp(b0,b1)]
        for b0, b1 in jxn_bps:
            if (b1,b0) in jxn_bps:
                jxn_bps.remove((b1,b0))
        jxn_bps.extend([arm.base_pairs[1] for arm in ns.arms.values()])
        # check if potential stackings are real
        for arm0, arm1 in ptt_aps:
            armbp0 = arm0.base_pairs[1]
            armbp1 = arm1.base_pairs[1]
            incyl_cnt = 0
            for jxn_bp in jxn_bps:
                if jxn_bp is armbp0 or jxn_bp is armbp1:
                    continue
                ang0 = obtain_cos(v_normalize(CoM_bp(jxn_bp)-CoM_bp(armbp0)),v_normalize(CoM_bp(armbp1)-CoM_bp(armbp0))) # angle~ {jxn_bp:armbp0, armbp1:armbp0}
                ang1 = obtain_cos(v_normalize(CoM_bp(jxn_bp)-CoM_bp(armbp1)),v_normalize(CoM_bp(armbp0)-CoM_bp(armbp1))) # angle~ {jxn_bp:armbp0, armbp1:armbp0}
                if ang0 > 90 or ang1 > 90:
                    continue
                # check if dist_ax < threshold. dist_ax: dist of CoM_bp from line connecting the two armbps.
                dist_ax = np.sin(ang0) * np.linalg.norm(CoM_bp(jxn_bp)-CoM_bp(armbp0))
                if dist_ax < 1.5 * 1.3: # * np.sin(15*np.pi/180) * 0.5 * incyl_cnt: # normal value: ~1.3 SU. 2 is ~150%*1.3 [* typical_stk_ang_corr_dist(sin(15deg)*incyl_cnt*typical_dist_stk_bps(0.5))] 
                    incyl_cnt += 1 # TODO! Need to know jxn# & arm#? currently it is a cone!
            # check if dist_1bps < typical_dist_stk_bps(0.5) * cnt_bps_incyl * threshold(150%) * tolerance(200%, optional).
            if arc_approx_corr(armbp0,armbp1,np.linalg.norm(CoM_bp(armbp0) - CoM_bp(armbp1))) < 0.5 * (incyl_cnt + 1) * 1.5 * 2: # do we want this arc_approx?
                # stacking detected!
                stack_ls.append((None, None, arm0.arm_id,arm1.arm_id)) # stack_ls.append((adj_bp, adj_bp2, arm.arm_id,arm2.arm_id))
        if len(stack_ls) == 0:
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
        else:
            is_stacking = True
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
    return stacking_dict

def check_if_bp(b1,b2):
    is_oppos_dir = np.dot(np.array(b1.backbone),np.array(b2.backbone)) < 0
    is_normal_align = np.dot(np.array(b1.normal),np.array(b2.normal)) < 0 # anti-aligned
    disp = (np.array(b1.position)-np.array(b2.position))
    is_same_plane_b1 = np.abs(np.dot(disp,np.array(b1.normal))) < 0.1
    is_same_plane_b2 = np.abs(np.dot(disp,np.array(b2.normal))) < 0.1
    if np.linalg.norm(disp) == 0:
        return False
    is_disp_small = np.linalg.norm(disp) < 2 # normal value: ~1.3 SU. 2 is ~150%*1.3
    disp /= np.linalg.norm(disp)
    is_align_disp_b1 = np.abs(np.dot(disp, np.array(b1.backbone))) > 0.5
    is_align_disp_b2 = np.abs(np.dot(disp, np.array(b2.backbone))) > 0.5
    return is_oppos_dir and is_normal_align and is_disp_small and (is_align_disp_b1 or is_align_disp_b2) and (is_same_plane_b1 or is_same_plane_b2)

def arc_approx_corr(armbp0,armbp1,norm):
    ang = obtain_cos(np.array(armbp0[0].normal),np.array(armbp1[0].normal))
    central_ang = (180-ang)*2
    r_ang = (ang/2-(180-ang)/2)
    r = norm/np.sin(central_ang*(np.pi/180)) * np.sin(r_ang*(np.pi/180)) # law of sines
    return r

def get_arm_dir(arm):
    dir_bp1_bp2 = v_normalize(CoM_bp(arm.base_pairs[2])-CoM_bp(arm.base_pairs[1]))
    is_b1_rev = np.dot(np.array(arm.base_pairs[1][1].normal),dir_bp1_bp2) > 0 # check if b1 of bp1 points in the direction of arm
    # define arm_dir as the avg of normal of bp1
    arm_dir = (np.array(arm.base_pairs[1][0].normal) * (is_b1_rev*2-1) * (-1) + np.array(arm.base_pairs[1][1].normal) * (is_b1_rev*2-1))/2 # if b1 not rev, b0 is rev.
    return v_normalize(arm_dir)

def CoM_bp(bp):
    return (np.array(bp[0].position)+np.array(bp[1].position))/2

def v_normalize(v):
    return v/np.linalg.norm(v)

def obtain_cos(v1,v2):
    return np.degrees(np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))) # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs

def stacking_local_identify_calc_og(path_top, path_traj, arm_num, dims_ls, ns_input = None, sys_input = None):
    # savepoint loading: strands-sys
    reader = Reader(path_top, path_traj)
    if type(sys_input) is str or sys_input == None:
        strands_tm = reader.read_data(p=sys_input)
    else:
        strands_tm = reader.read_data(obj=sys_input)
    # savepoint loading: nano-stars
    nc = NanoConstructor(strands_tm, dims_ls, arm_num)
    if type(ns_input) is str or ns_input == None: 
        # box_dim hacking
        import re
        with open(path_traj,'r') as f:
            f.readline()
            ret=re.match('^b = ([0-9]+) ([0-9]+) ([0-9]+)\n',f.readline())
        box_dim = np.array((ret.group(1),ret.group(2),ret.group(3)))
        ns_tm = nc.construct(p=ns_input, box_dim=box_dim)
    else:
        ns_tm = nc.construct(obj=ns_input)
    # finish savepoint loading

    from copy import deepcopy
    stacking_dict = OrderedDict() # OrderedDict{t_stamp: is_stacking, stack_ls}
    for t_stamp, ns in ns_tm.time_capsule.items():
        is_stacking = False
        stack_ls = [] # [(adj_bp, adj_bp2, arm_id, arm2_id)]
        # obtain the paired 1st bps
        arm_ls = []
        for _,arm in ns.arms.items():
            f_bp = arm.base_pairs[1] # idx count from center, starting from 1. type(f_bp) == tuple (bp,bp)
            # check if broken
            is_oppos_dir = np.dot(np.array(f_bp[0].backbone),np.array(f_bp[1].backbone)) < 0
            is_normal_align = np.dot(np.array(f_bp[0].normal),np.array(f_bp[1].normal)) < 0 # anti-aligned
            disp = (np.array(f_bp[0].position)-np.array(f_bp[1].position))
            is_disp_small = np.linalg.norm(disp) < 2 # normal value: ~1.3 SU. 2 is ~150%*1.3
            disp /= np.linalg.norm(disp)
            is_align_disp_b0 = np.abs(np.dot(disp, np.array(f_bp[0].backbone))) > 0.5
            is_align_disp_b1 = np.abs(np.dot(disp, np.array(f_bp[1].backbone))) > 0.5
            if is_oppos_dir and is_normal_align and is_disp_small and (is_align_disp_b0 or is_align_disp_b1):
                arm_ls.append(arm)
        if len(arm_ls) <= 1:
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
        # find adjacent base, check if forming bp
        stack_cand_ls = [] # [(arm,adj_bp)]
        for arm in arm_ls:
            f_bp = arm.base_pairs[1]
            if f_bp[0].prev_id in ns.center.keys() or f_bp[0].next_id in ns.center.keys():
                adj_b0 = ns.center[f_bp[0].prev_id] if f_bp[0].prev_id in ns.center.keys() else ns.center[f_bp[0].next_id]
            else:
                assert 0 == 1 # safety check
            if f_bp[1].prev_id in ns.center.keys() or f_bp[1].next_id in ns.center.keys():
                adj_b1 = ns.center[f_bp[1].prev_id] if f_bp[1].prev_id in ns.center.keys() else ns.center[f_bp[1].next_id]
            else:
                assert 0 == 1
            # check if bbv aligning
            is_bbv_align_b0 = np.dot(np.array(adj_b0.backbone),np.array(f_bp[0].backbone)) > 0
            is_bbv_align_b1 = np.dot(np.array(adj_b1.backbone),np.array(f_bp[1].backbone)) > 0
            if not (is_bbv_align_b0 and is_bbv_align_b1):
                continue
            # check if adj_b0 and adj_b1 pair
            is_oppos_dir = np.dot(np.array(adj_b0.backbone),np.array(adj_b1.backbone)) < 0
            is_normal_align = np.dot(np.array(adj_b0.normal),np.array(adj_b1.normal)) < 0 # anti-aligned
            disp = (np.array(adj_b0.position)-np.array(adj_b1.position))
            is_disp_small = np.linalg.norm(disp) < 2 # normal value: ~1.3 SU. 2 is ~150%*1.3
            disp /= np.linalg.norm(disp)
            is_align_disp_b0 = np.abs(np.dot(disp, np.array(adj_b0.backbone))) > 0.5
            is_align_disp_b1 = np.abs(np.dot(disp, np.array(adj_b1.backbone))) > 0.5
            if is_oppos_dir and is_normal_align and is_disp_small and (is_align_disp_b0 or is_align_disp_b1):
                stack_cand_ls.append((arm,(adj_b0,adj_b1)))
        # find stacking
        if len(stack_cand_ls) < 2:
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
        cand_op_ls = deepcopy(stack_cand_ls)
        for arm, adj_bp in stack_cand_ls:
            if (arm, adj_bp) not in cand_op_ls:
                continue
            else:
                cand_op_ls.remove((arm, adj_bp))
            # find the other bulk stacking with the selected: f2_bp & adj2_bp vs f_bp & adj_bp
            for arm2, adj_bp2 in cand_op_ls:
                # check if normal align
                cross = np.cross(np.array(adj_bp[0].backbone),np.array(adj_bp[1].backbone))
                cross /= np.linalg.norm(cross)
                cross2 = np.cross(np.array(adj_bp2[0].backbone),np.array(adj_bp2[1].backbone))
                cross2 /= np.linalg.norm(cross2)
                if np.abs(np.dot(cross,cross2)) > 0.7:
                    # align detected!
                    stack_ls.append((adj_bp, adj_bp2, arm.arm_id,arm2.arm_id))
                    cand_op_ls.remove((arm2, adj_bp2))
                    break
        if len(stack_ls) == 0:
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
        else:
            is_stacking = True
            stacking_dict[t_stamp]=(is_stacking, stack_ls)
            continue
    return stacking_dict

