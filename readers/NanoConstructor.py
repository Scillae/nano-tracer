from models import Strand, Base, NanoStar, NanoMesh, Arm, TimeMachine
from collections import OrderedDict
from utils.tools import save_load as SL
import numpy as np
import copy
import pickle
import os.path

class NanoConstructor:
    def __init__(self, strands_tm, dims_ls, arm_num):
        self.time_machine = None
        self.dim = dims_ls
        self.box_dim = None # box_dim hacking
        self.strands = strands_tm
        self.arm_num = arm_num # nanomesh not using this.
    
    def save_load(self, p, obj = None):
        r_obj = SL(p, obj)
        return r_obj

    def construct(self, p = None, obj = None, box_dim = None): # # box_dim hacking
        self.time_machine = self.save_load(p, obj)
        if self.time_machine == False:
            if p == None:
                self.time_machine = TimeMachine()
            else:
                self.time_machine = TimeMachine(p)
            # check if nanomesh
            if len(self.strands.time_capsule[self.strands.timeseries[0]]) > self.arm_num:
                for t_stamp in self.strands.timeseries:
                    nm = NanoMesh(self.strands.time_capsule[t_stamp],self.dim)
                    self.time_machine.add_strands(t_stamp, nm)
            else:
                self.time_machine.box_dim = box_dim # box_dim. A hacking solution.
                for t_stamp in self.strands.timeseries:
                    ns = NanoStar(self.strands.time_capsule[t_stamp],self.dim, self.arm_num, box_dim=box_dim)
                    self.time_machine.add_strands(t_stamp, ns)
            self.time_machine = self.save_load(p, self.time_machine)
        return self.time_machine


