from models import Strand, Base, NanoStar, Arm, TimeMachine
from collections import OrderedDict
import numpy as np
import copy

def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))


class NanoMesh:
    def __init__(self, strands_dic, dims_ls):
        self.dim = dims_ls # dims_ls
        self.homogen_arms = True
        self.nodes = {} # {ns#: (ns_obj, arm#)}
        self.conn_map = None
        # construct all nanostars
        # for all single-ends, (energy or topology) determine if bonded; if bonded: find nearest base -> find strand -> find nanostar
        # iteratively fill up conn_map
        len_arm, len_cen, len_end = self.dim
        # construct nanostars
        strand_ids = list(strands_dic.keys())
        # pairing_pos = [(dims_ls[0]-dims_ls[2])//2 - (dims_ls[0]-dims_ls[2])//4, (dims_ls[0]-dims_ls[2])//2, (dims_ls[0]-dims_ls[2])//2 + (dims_ls[0]-dims_ls[2])//4]
        # pairing_extra_pos = [i for i in range((dims_ls[0]-dims_ls[2])//2 - (dims_ls[0]-dims_ls[2])//4, (dims_ls[0]-dims_ls[2])//2 + (dims_ls[0]-dims_ls[2])//4 + 1)]
        search_pts_num = 5 # amount of points to be used for paring.
        pairing_pos = [i*dims_ls[0]//search_pts_num for i in range(search_pts_num)]
        pairing_extra_pos = [i for i in range(0, dims_ls[0])]
        # start
        discard_cnt = 0 # debug ~ #ns instances discarded
        strand_groups = [] # [[strand_id]] of all nanostars
        emp_dist = 0 # temporary variable for Typical Distance between Paired Strands
        while len(strand_ids) > 0:
            strand_g = [] # [strand_id] that compose 1 nanostar
            s0 = strands_dic[strand_ids[0]] # s0 is strand, not list
            strand_g.append(s0.strand_id)
            # traverse all strands left
            while True:
                pairing_ls = []
                # Loose Pairing method.
                for loc in pairing_pos:
                    s0_b = list(s0.base_sequence.values())[loc]
                    pairing = (10000000000, None)
                    for i in strand_ids:
                        s1_b = list(strands_dic[i].base_sequence.values())[len_arm+len_arm+len_cen-1-loc]
                        dis = dist(s0_b.position, s1_b.position)
                        if dis < pairing[0]:
                            pairing = (dis, s1_b.strand_id)
                    pairing_ls.append(pairing[1])
                # safeguard check: consistency of pairing. 
                # Strict Pairing method.
                if len(set(pairing_ls)) != 1:
                    print('Pairing disturbed. Strict Pairing method.')
                    # 3 points not agreeing, search for a larger range.
                    pairing_ls = []
                    for loc in pairing_extra_pos:
                        s0_b = list(s0.base_sequence.values())[loc]
                        pairing = (10000000000, None)
                        for i in strand_ids:
                            s1_b = list(strands_dic[i].base_sequence.values())[len_arm+len_arm+len_cen-1-loc]
                            dis = dist(s0_b.position, s1_b.position)
                            if dis < pairing[0]:
                                if emp_dist == 0 or (dis < emp_dist*1.1 and dis > emp_dist*0.9):
                                    pairing = (dis, s1_b.strand_id)
                                else:
                                    print(f'abnormal distance {dis}, pairing: {s0_b.strand_id}~{s1_b.strand_id} is discarded.')
                        pairing_ls.append(pairing[1])
                elif emp_dist == 0: # elif
                    p_id = max(pairing_ls, key=pairing_ls.count)
                    emp_dist_ls = []
                    s0_b_ls = list(s0.base_sequence.values())
                    s1_b_ls = list(strands_dic[p_id].base_sequence.values())
                    for i in range(0, dims_ls[0]):
                        b0 = s0_b_ls[i]
                        b1 = s1_b_ls[len_arm+len_arm+len_cen-1-i]
                        emp_dist_ls.append(dist(b0.position, b1.position))
                    avg_dist = sum(emp_dist_ls)/(len(emp_dist_ls))
                    if (max(emp_dist_ls)-min(emp_dist_ls)) < avg_dist*0.15:
                        if avg_dist < 2.0:
                            emp_dist = avg_dist
                        else:
                            print(f'Large pairing distance: {avg_dist} of a successfully paired arm is detected! Discarded!')
                    # else:
                        # print(f'abnormal pairing! {s0_b.strand_id}~{s1_b.strand_id} is loose!')
                # safeguard ends.

                # sanity check: bad-formed loose nanostar. only reporting.
                if max(pairing_ls, key=pairing_ls.count) == None:
                # if None in pairing_ls:
                    print('Bad-formed loose nanostar detected! Majority of the pairing distances are abnormal!')
                    # debug ~ forcing 4 arms
                    print(f'#### DEBUG #### Malformed 4 arm nanostar! [Broken Arms] Group {strand_g} Discarded!')
                    del strand_ids[0]
                    break
                    # debug ends.
                p_ls = [i for i in pairing_ls if i != None]
                p_id = max(p_ls, key=p_ls.count)
                if p_id in strand_g:
                    # looped back. a completed nanostar.
                    strand_ids.remove(p_id)
                    # debug ~ forcing 4 arms
                    if len(strand_g) == 4:
                        strand_groups.append(strand_g)
                    else:
                        print(f'#### DEBUG #### Malformed 4 arm nanostar! [Less Arms] Group {strand_g} Discarded!')
                    # strand_groups.append(strand_g)
                    # debug ends.
                    break
                strand_g.append(p_id)
                # debug ~ forcing 4 arms
                if len(strand_g) == 5:
                    print(f'#### DEBUG #### Malformed 4 arm nanostar! [More Arms] Group {strand_g} Discarded!')
                    del strand_ids[0]
                    break
                # debug ends.
                strand_ids.remove(p_id)
                s0 = strands_dic[p_id]
        ns_idx = 0
        arm_num_ls = []
        for groups in strand_groups:
            strands_nano = {}
            for i in groups:
                strands_nano[i] = strands_dic[i]
            arm_num = len(strands_nano)
            ns = NanoStar(strands_nano, self.dim, arm_num)
            self.nodes[ns_idx] = (ns, arm_num)
            arm_num_ls.append(arm_num)
        if len(set(arm_num_ls)) != 1:
            self.homogen_arms == False
            print('Heterogenous Mesh detected!!')
            print(f'Mesh Arm# list: {arm_num_ls}')
        # construct conn_map TODO.


'''
class NanoMesh:
    def __init__(self, strands_dic, dims_ls):
        self.dim = dims_ls # dims_ls
        self.homogen_arms = True
        self.nodes = {} # {ns#: (ns_obj, arm#)}
        self.conn_map = None
        # construct all nanostars
        # for all single-ends, (energy or topology) determine if bonded; if bonded: find nearest base -> find strand -> find nanostar
        # iteratively fill up conn_map
        len_arm, len_cen, len_end = self.dim
        # construct nanostars
        strand_ids = list(strands_dic.keys())
        # pairing_pos = [(dims_ls[0]-dims_ls[2])//2 - (dims_ls[0]-dims_ls[2])//4, (dims_ls[0]-dims_ls[2])//2, (dims_ls[0]-dims_ls[2])//2 + (dims_ls[0]-dims_ls[2])//4]
        # pairing_extra_pos = [i for i in range((dims_ls[0]-dims_ls[2])//2 - (dims_ls[0]-dims_ls[2])//4, (dims_ls[0]-dims_ls[2])//2 + (dims_ls[0]-dims_ls[2])//4 + 1)]
        search_pts_num = 5 # amount of points to be used for paring.
        pairing_pos = [i*dims_ls[0]//search_pts_num for i in range(search_pts_num)]
        pairing_extra_pos = [i for i in range(0, dims_ls[0])]
        # start
        discard_cnt = 0 # debug ~ #ns instances discarded
        strand_groups = [] # [[strand_id]] of all nanostars
        emp_dist = 0 # temporary variable for Typical Distance between Paired Strands
        while len(strand_ids) > 0:
            strand_g = [] # [strand_id] that compose 1 nanostar
            s0 = strands_dic[strand_ids[0]] # s0 is strand, not list
            strand_g.append(s0.strand_id)
            # traverse all strands left
            while True:
                pairing_ls = []
                # Loose Pairing method.
                for loc in pairing_pos:
                    s0_b = list(s0.base_sequence.values())[loc]
                    pairing = (10000000000, None)
                    for i in strand_ids:
                        s1_b = list(strands_dic[i].base_sequence.values())[len_arm+len_arm+len_cen-1-loc]
                        dis = dist(s0_b.position, s1_b.position)
                        # temporarily fixed emp_dist & enabled emp_dist in loose paring.
                        # if dis < pairing[0]:
                        if dis < pairing[0] and (emp_dist == 0 or dis < emp_dist):
                            pairing = (dis, s1_b.strand_id)
                    pairing_ls.append(pairing[1])
                # # safeguard check: consistency of pairing. 
                # # Strict Pairing method.
                # if len(set(pairing_ls)) != 1:
                #     print('Pairing disturbed. Strict Pairing method.')
                #     # 3 points not agreeing, search for a larger range.
                #     pairing_ls = []
                #     for loc in pairing_extra_pos:
                #         s0_b = list(s0.base_sequence.values())[loc]
                #         pairing = (10000000000, None)
                #         for i in strand_ids:
                #             s1_b = list(strands_dic[i].base_sequence.values())[len_arm+len_arm+len_cen-1-loc]
                #             dis = dist(s0_b.position, s1_b.position)
                #             if dis < pairing[0]:
                #                 if emp_dist == 0 or (dis < emp_dist*1.1 and dis > emp_dist*0.9):
                #                     pairing = (dis, s1_b.strand_id)
                #                 else:
                #                     print(f'abnormal distance {dis}, pairing: {s0_b.strand_id}~{s1_b.strand_id} is discarded.')
                #         pairing_ls.append(pairing[1])
                # if emp_dist == 0: # elif
                #     p_id = max(pairing_ls, key=pairing_ls.count)
                #     emp_dist_ls = []
                #     s0_b_ls = list(s0.base_sequence.values())
                #     s1_b_ls = list(strands_dic[p_id].base_sequence.values())
                #     for i in range(0, dims_ls[0]):
                #         b0 = s0_b_ls[i]
                #         b1 = s1_b_ls[len_arm+len_arm+len_cen-1-i]
                #         emp_dist_ls.append(dist(b0.position, b1.position))
                #     avg_dist = sum(emp_dist_ls)/(len(emp_dist_ls))
                #     if (max(emp_dist_ls)-min(emp_dist_ls)) < avg_dist*0.15:
                #         if avg_dist < 2.0:
                #             emp_dist = avg_dist
                #         else:
                #             print(f'Large pairing distance: {avg_dist} of a successfully paired arm is detected! Discarded!')
                #     else:
                #         print(f'abnormal pairing! {s0_b.strand_id}~{s1_b.strand_id} is loose!')
                # # safeguard ends.

                # sanity check: bad-formed loose nanostar. only reporting.
                # if max(pairing_ls, key=pairing_ls.count) == None:
                if None in pairing_ls:
                    print('Bad-formed loose nanostar detected! Majority of the pairing distances are abnormal!')
                    # debug ~ forcing 4 arms
                    print(f'#### DEBUG #### Malformed 4 arm nanostar! [Broken Arms] Group {strand_g} Discarded!')
                    del strand_ids[0]
                    break
                    # debug ends.
                p_ls = [i for i in pairing_ls if i != None]
                p_id = max(p_ls, key=p_ls.count)
                if p_id in strand_g:
                    # looped back. a completed nanostar.
                    strand_ids.remove(p_id)
                    # debug ~ forcing 4 arms
                    if len(strand_g) == 4:
                        strand_groups.append(strand_g)
                    else:
                        print(f'#### DEBUG #### Malformed 4 arm nanostar! [Less Arms] Group {strand_g} Discarded!')
                    # strand_groups.append(strand_g)
                    # debug ends.
                    break
                strand_g.append(p_id)
                # debug ~ forcing 4 arms
                if len(strand_g) == 5:
                    print(f'#### DEBUG #### Malformed 4 arm nanostar! [More Arms] Group {strand_g} Discarded!')
                    del strand_ids[0]
                    break
                # debug ends.
                strand_ids.remove(p_id)
                s0 = strands_dic[p_id]
        ns_idx = 0
        arm_num_ls = []
        for groups in strand_groups:
            strands_nano = {}
            for i in groups:
                strands_nano[i] = strands_dic[i]
            arm_num = len(strands_nano)
            ns = NanoStar(strands_nano, self.dim, arm_num)
            self.nodes[ns_idx] = (ns, arm_num)
            arm_num_ls.append(arm_num)
        if len(set(arm_num_ls)) != 1:
            self.homogen_arms == False
            print('Heterogenous Mesh detected!!')
            print(f'Mesh Arm# list: {arm_num_ls}')
        # construct conn_map TODO.
        
'''
