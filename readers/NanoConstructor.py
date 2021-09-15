from models import Strand, Base, NanoStar, NanoMesh, Arm, TimeMachine
from collections import OrderedDict
import numpy as np
import copy
import pickle
import os.path

class NanoConstructor:
    def __init__(self, strands_tm, dims_ls, arm_num):
        self.time_machine = None
        self.dim = dims_ls
        self.strands = strands_tm
        self.arm_num = arm_num # nanomesh not using this.
    
    def save_load(self, p, obj):
        print(f'TM path: {p}')
        if p is None:
            print('TM skipping!')
            return obj
        if (obj is not None) and (not os.path.isfile(p)):
            if os.path.isdir(os.path.split(p)[0]) == False:
                os.makedirs(os.path.split(p)[0])
            print('TM saving!')
            pickle.dump(obj, open(p,"wb"))
            r_obj = obj
        elif (obj is None) and os.path.isfile(p):
            print('TM loading!')
            r_obj = pickle.load(open(p, "rb"))
        elif (obj is not None) and os.path.isfile(p):
            print('TM updating savepoint!')
            pickle.dump(obj, open(p,"wb"))
            r_obj = obj
        else:
            print('TM save_load both empty')
            r_obj = False
        return r_obj

    def construct(self, p = None, obj = None):
        self.time_machine = self.save_load(p, obj)
        if self.time_machine == False:
            self.time_machine = TimeMachine()
            # check if nanomesh
            if len(self.strands.time_capsule[self.strands.timeseries[0]]) > self.arm_num:
                for t_stamp in self.strands.timeseries:
                    nm = NanoMesh(self.strands.time_capsule[t_stamp],self.dim)
                    self.time_machine.add_strands(t_stamp, nm)
            else:
                for t_stamp in self.strands.timeseries:
                    ns = NanoStar(self.strands.time_capsule[t_stamp],self.dim, self.arm_num)
                    self.time_machine.add_strands(t_stamp, ns)
            self.time_machine = self.save_load(p, self.time_machine)
        return self.time_machine


