from models import Strand, Base, Arm
from collections import OrderedDict
import numpy as np
import copy

# Noella: Accepted_sci
# dist function should be defined before using
def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))

class NanoStar:
    def __init__(self, strands_dic, dims_ls, arm_num):
        # revert all strands: done in TimeMachine.add_strands
        # strands_dic = OrderedDict()
        # sanity check of strand amount
        assert len(strands_dic) == arm_num
        strands = copy.deepcopy(strands_dic) # avoid modification to the ref var.
        '''
        for strand_id, strand in strands_dic.items():
            # Noella: Accepted_sci
            # `list.reverse()` return nothing, thus `base_seq_r_ls` will be assigned to None
            # base_seq_r_ls = list(strand.base_sequence.items()).reverse()
            base_seq_r_ls = list(strand.base_sequence.items())
            base_seq_r_ls.reverse()
            strands_dic[strand_id] = Strand(strand_id, base_seq_r_ls)
        '''
        self.strands = strands_dic
        self.center = self.center_gen(strands, dims_ls)
        self.arms = self.binding(strands, dims_ls)

    def center_gen(self, strands_dic, dims_ls):
        len_arm, len_cen, len_end = dims_ls
        cen_dic = OrderedDict()
        st = 0
        # Noella: Accepted_sci
        # for strand in strands_dic:
        for _, strand in strands_dic.items():
            for i in range(len_cen):
                # Noella: Questioned_sci(base sequence is a dic, and thus should not be accessed by index.)
                # `strand.base_sequence` is a list defined in line 15, `base_seq_r_ls` (a list) is passed to contructor as required parameter `base_seq`. All `strand.base_sequence.values()` should be replaced
                # around lines: 31 44 48 68 71 (90) [may not be exactly because of new added comments lines]
                # cen_dic[st*len_cen+i] = strand.base_sequence.values()[len_arm+i]
                # dic or list: dict
                cen_dic[st*len_cen+i] = list(strand.base_sequence.values())[len_arm+i]
            st += 1
        return cen_dic

    def binding(self, strands_dic, dims_ls):
        strand_id_idx = list(strands_dic.keys())
        # find 3' head of an arbitrary strand
        # len_arm, len_cen, len_end = dims_ls
        s0_i = min(list(strands_dic.keys()))
        s0 = strands_dic[s0_i]
        # Noella: Questioned_sci(0~7 after s1 being reversed is correct.)
        # s0_i should not be remove ? index 0 may pair index 6 ?
        # strand_id_idx.remove(s0_i)
        # Noella: Questioned_sci(dic or list) - dict
        # s0_head = s0.base_sequence[0]
        s0_head = list(s0.base_sequence.values())[0]
        # pool all bases with their locations
        pool_base_ls = []
        for idx in strand_id_idx:
            if idx == s0.strand_id:
                continue
            # Noella: Questioned_sci(dic or list) - dict
            # pool_base_ls.extend(list(strands_dic[idx].base_sequence))
            pool_base_ls.extend(list(strands_dic[idx].base_sequence.values()))
        # get the closest base
        min_dist = 1000000
        min_loc = 0
        for i in range(len(pool_base_ls)):
            # Noella: Accepted_sci
            # using duplicate name `dist` is not recommended
            #   dist = dist(pool_base_ls[i].position,s0_head.position)
            dis = dist(pool_base_ls[i].position,s0_head.position)
            if dis < min_dist:
                min_dist = dis
                min_loc = i
        bind_base = pool_base_ls[min_loc]
        # assert its relative location in the strand is len(strand) - len_end
        # construct arm
        s1 = strands_dic[bind_base.strand_id]
        arm_idx = 0
        # Noella: Accepted_sci
        # `arm_dic` is not defined
        arm_dic = OrderedDict()
        arm_dic[arm_idx] = Arm(arm_idx, dims_ls, s0.strand_id, {s0.strand_id:s0,s1.strand_id:s1})
        # loop out all other strands
        while len(strands_dic) != 0:
            s0 = copy.deepcopy(s1)
            # Noella: Questioned_sci(True, but may be further improved by first remove s0_i and then add it back at the last cycle. Circle topology must close at the starting pt.)
            # Maybe Scillae want to remove current strand id from `strand_id_idx` rather than `s0_i`
            # strand_id_idx.remove(s0_i)
            strand_id_idx.remove(s0.strand_id)
            # Noella: Questioned_sci(dic or list) - dict
            # s0_head = s0.base_sequence[0]
            s0_head = list(s0.base_sequence.values())[0]
            pool_base_ls = []
            # Noella: Questioned_sci(why not using len(strands_dic)?)
            # We may need to break the loop if strand_id_idx == []
            if strand_id_idx == []:
                break
            for idx in strand_id_idx:
                # Noella: Questioned_sci(dic or list) - dict
                # pool_base_ls.extend(list(strands_dic[idx].base_sequence))
                pool_base_ls.extend(list(strands_dic[idx].base_sequence.values()))
            min_dist = 1000000
            min_loc = 0
            for i in range(len(pool_base_ls)):
                # Noella: Accepted_sci
                # duplicate name `dist`
                dis = dist(pool_base_ls[i].position,s0_head.position)
                if dis < min_dist:
                    min_dist = dis
                    min_loc = i
            bind_base = pool_base_ls[min_loc]
            s1 = strands_dic[bind_base.strand_id]
            arm_idx += 1
            arm_dic[arm_idx] = Arm(arm_idx, dims_ls, s0.strand_id, {s0.strand_id:s0,s1.strand_id:s1})
            del strands_dic[s0.strand_id]
        return arm_dic
        '''
        len_arm, len_cen, len_end = dims_ls
        len_tot = 2*len_arm + len_cen + len_end
        binding_pts_2dls = []
        for strand in strands_dic.values():
            str_ls = list(strand.base_sequence)
            binding_pts_2dls.append(list(str_ls[0],str_ls[len_end],str_ls[len_tot-1-len_end],str_ls[len_tot-1]))
        # now find the min-dist
        #TODO: THE LOOP HERE MUST BE REDESIGNED!!!
        bind_tp_ls = []
        for i in range(len(binding_pts_2dls)):
            idx_arr = np.arange(len(binding_pts_2dls))
            idx_arr = idx_arr[idx_arr != i] # double-counting if '!=' ? Must Double Count! No Need to double count!
            min_dist_dic = {}
            for j in idx_arr: # still need to compare 'ALL' pairs of DNAs
                s0 = binding_pts_2dls[i]
                s1 = binding_pts_2dls[j]
                min_dist_pair_dic = {}
                for k in range(len(s0)):
                    for l in range(k+1, len(s1)):
                        min_dist_pair_dic[dist(s0[k], s1[l])] = [s0[k], s1[l]]
                        # min_dist_pair_ls.append((dist(s0[k], s1[l]), s0[k], s1[l]))
                # min_dist_ls.append(min(min_dist_pair_ls))
                m_d = min(list(min_dist_pair_dic.keys()))
                m_v = min_dist_pair_dic[m_d]
                min_dist_dic[j] = m_v.insert(0, m_d) # Find min in 16 pts. This generate tot_i * tot_j min. In each i, 2 comparable min out of tot_j should be found.
                # min_dist_dic: 5 items, each j: [s0_b, s1_b, dist]
            min_d_0 = min([m_v_ls[2] for m_v_ls in list(min_dist_dic.values())])
            min_item_0 = [item for item in min_dist_dic.items() if item[1][2] == min_d_0][0] # select tuples (j,[s0,s1,d]) if d==d_0
            min_dist_dic.pop(min_item_0[0])
            min_d_1 = min([m_v_ls[2] for m_v_ls in list(min_dist_dic.values())])
            min_item_1 = [item for item in min_dist_dic.items() if item[1][2] == min_d_0][0] # 2 contact points for each strand with all other ones
        '''
        '''
            min_dist_cp = min_dist_ls[:] # avoid potential reference/shallow copy issue
            min_2dists_ls = [tp for tp in min_dist_cp if not any(tp[0] > t[0] for t in min_dist_ls)] # unsure about this line
            min_dist_cp = [tp for tp in min_dist_ls if tp not in min_2dists_ls] # drop lowest
            min_2dists_ls.extend([tp for tp in min_dist_cp if not any(tp[0] < t[0] for t in min_dist_cp)]) # extend only 1 element.
        '''
        '''
                s0 = binding_pts_2dls[i]
                s1 = binding_pts_2dls[j]
                dist = [dist(s0[0].position, s1[2].position),
                        dist(s0[0].position, s1[1].position),
                        dist(s0[1].position, s1[0].position),
                        dist(s0[1].position, s1[3].position),
                        dist(s0[2].position, s1[3].position),
                        dist(s0[2].position, s1[0].position),
                        dist(s0[3].position, s1[1].position),
                        dist(s0[3].position, s1[2].position),
                ]
                min_dist_val = min(dist)
                min_dist_idx = dist.index(min_dist_val)
                # check if actually binding #
                dist.sort()
                if dist[-2] > dist[-1]*2:
                    continue # break out to the next
                # end check #
                # 2 binding, possiblilities: 0~0, 0~41, 41~0, 41~41
                # NOTE if strands are reversed here, they would be reversed again since each strand counts twice? No! No double count! But I cannot reverse it here...
                if s0[min_dist_idx].base_id == s0[0].base_id: #s0 leading normal
                    if s1[min_dist_idx].base_id == s1[1].base_id: #s1 normal
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (strands_dic[s0[min_dist_idx].strand_id], strands_dic[s1[min_dist_idx].strand_id]), (False, False))
                    else: # s1[3], s1 reversed
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (strands_dic[s0[min_dist_idx].strand_id], OrderedDict(reversed(list(strands_dic[s1[min_dist_idx].strand_id].items())))), (False, True))
                    
                elif s0[min_dist_idx].base_id == s0[3].base_id: #s0 leading reversed
                    if s1[min_dist_idx].base_id == s1[1].base_id: #s1 normal
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (OrderedDict(reversed(list(strands_dic[s0[min_dist_idx].strand_id].items()))), strands_dic[s1[min_dist_idx].strand_id]), (True, False))
                    else: # s1[3]
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (OrderedDict(reversed(list(strands_dic[s0[min_dist_idx].strand_id].items()))), OrderedDict(reversed(list(strands_dic[s1[min_dist_idx].strand_id].items())))), (True, True))

                elif s1[min_dist_idx].base_id == s0[2].base_id: #s0 ending normal

                    











        tp_base_ls = []
        bt_base_ls = []
        for strand_tp in strands_dic.items():
            (strand_id, strand) = strand_tp
            bt_base_ls.append(strand.base_sequence[max(list(strand.base_sequence.keys()))])
            tp_base_ls.append(strand.base_sequence[min(list(strand.base_sequence.keys()))])
        bind_ls = []
        for tp_base in tp_base_ls:
            pair_tp = (0,0,0)
            min_dist_dbl = 100000000000000000.00 #inf
            for bt_base in bt_base_ls:
                dist_dbl = np.sqrt(np.sum(np.square(np.array(tp_base.position) - np.array(bt_base.position))))
                assert min_dist_dbl != dist_dbl
                if dist_dbl < min_dist_dbl:
                    min_dist_dbl = dist_dbl
                    pair_tp = (tp_base,bt_base,dist_dbl)
            bind_ls.append(pair_tp)
        print(bind_ls)
        '''



'''
import numpy as np
from models import Arm
class NanoStar:
    def __init__(self, strands, dims_ls, bind_ls):
        self.arms = self.construct_arms(strands, dims_ls, bind_ls)
        self.single_ends = self.get_single_ends()

    # construct arms, main work is to determine strand pair,
    # then construct by Arm's ctor
    # TODO how to present centers: move Arm.pair() here, and construct Center class.
    def construct_arms(self, strands, dims_ls, bind_ls):
        arms_ls = []
        for i in range(len(bind_ls)):
            tp_base, bt_base, _ = bind_ls[i]
            a = Arm(i, dims_ls,[strands[tp_base.strand_id], strands[bt_base.strand_id]])
            arms_ls.append(a)
        return arms_ls


    # get single ends from self.arms
    def get_single_ends(self):
        pass
'''