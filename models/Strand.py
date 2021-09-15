from collections import OrderedDict
class Strand:
    # Noella: Accepted_sci
    # 修个怪bug, 如果函数默认值在定义时构建的话, 其id会和类绑定, 变成类似于静态成员的东西 (特性 or BUG?)
    # def __init__(self, strand_id, base_seq = OrderDict()):
    def __init__(self, strand_id, base_seq = None, timestamp = None):
        self.strand_id = strand_id
        self.timestamp = timestamp
        # list or dic?
        self.base_sequence = OrderedDict({}) if base_seq is None else base_seq

    def add_base(self, base_entity):
        self.base_sequence[base_entity.base_id] = base_entity

    @staticmethod
    def parse(s):
        pass