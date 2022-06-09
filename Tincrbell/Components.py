from .Function import Function
from .F2n import F2n
from pysat.card import CardEnc, EncType

class __INVERT__(Function):
    def __init__(self):
        super().__init__(1, 1)

    def __call__(self, v: F2n):
        return ~v

    def get_ANF_propagation_model(self, input_key_mask, output_key_mask):
        return ((-1, 2), )


INVERT = __INVERT__()


class __XOR__(Function):
    def __init__(self):
        super().__init__(2, 1)
        self.ANF_prop_model = ((1, 2, -3), (-1, 3), (-2, 3), (-1, -2))

    def __call__(self, v: F2n):
        x1, x2 = v.split(1)
        return x1 ^ x2

    def get_ANF_propagation_model(self, input_key_mask, output_key_mask):
        if input_key_mask.hamming_weight() == output_key_mask.hamming_weight() == 0:
            return ((1, 2, -3), (-1, 3), (-2, 3), (-1, -2))
        else:
            return ((-1, 3), (-2, 3), (-1, -2))



XOR = __XOR__()

class XORn(Function):

    def __init__(self, n, encoding=EncType.pairwise):
        super().__init__(n, 1)
        self.atmost1 = tuple(tuple(x) for x in CardEnc.atmost(lits = list(range(1, n+1)), encoding = encoding).clauses)
        self.n_cnf_vars = max(map(lambda x: max(map(abs, x)), self.atmost1)) + 1
    
    def get_output_vars(self) -> tuple[int]:
        return (self.n_cnf_vars, )

    def __call__(self, v: F2n):
        return F2n(1, v.hamming_weight() % 2)

    def get_ANF_propagation_model(self, input_key_mask, output_key_mask):
        if input_key_mask.hamming_weight() == output_key_mask.hamming_weight() == 0:
            return self.atmost1 +  tuple((-v, self.n_cnf_vars) for v in self.get_input_vars()) + (self.get_input_vars() + (-self.n_cnf_vars,),)
        else:
            return self.atmost1 + tuple((-v, self.n_cnf_vars) for v in self.get_input_vars())


class __AND__(Function):
    def __init__(self):
        super().__init__(2, 1)
        self.ANF_prop_models = {(F2n(2, 0), F2n(1, 0)):((1, -3), (2, -1), (3, -2)),
                                (F2n(2, 0), F2n(1, 1)):((1, -2), (2, -1), (3, -2)),
                                (F2n(2, 1), F2n(1, 0)):((2, -1), (2, -3), (3, -2)),
                                (F2n(2, 1), F2n(1, 1)):((2, -1), (3, -2)),
                                (F2n(2, 2), F2n(1, 0)):((1, -2), (1, -3), (3, -1)),
                                (F2n(2, 2), F2n(1, 1)):((1, -2), (3, -1)),
                                (F2n(2, 3), F2n(1, 0)):((3, -1), (3, -2)),
                                (F2n(2, 3), F2n(1, 1)):((3, -1), (3, -2))}

    def __call__(self, v: F2n):
        x1, x2 = v.split(1)
        return x1 & x2


AND = __AND__()


class COPYn(Function):
    def __init__(self, n):
        self.n = n
        super().__init__(1, n)
        self.ANF_prop_model = tuple((1, -x-2) for x in range(self.n)) + (tuple(range(2, self.n+2)) + (-1, ), )

    def __call__(self, v: F2n):
        r = v.concatenate(v)
        for _ in range(self.n-2):
            r = r.concatenate(v)
        return r

    def get_ANF_propagation_model(self, input_key_mask, output_key_mask):
        if input_key_mask[0] == 1:
            return (tuple(range(2, self.n+2)) + (-1, ), )
        else:
            return tuple((1, -x-2) for x in range(self.n) if output_key_mask[x]==0) + (tuple(range(2, self.n+2)) + (-1, ), )


class SBox(Function):

    def __init__(self, input_size, output_size, lookup_table: list):
        self.lookup_table = lookup_table
        super().__init__(input_size, output_size)

    def __call__(self, v: F2n):
        return F2n(self.output_size, self.lookup_table[v])

class Sink(Function):
    def __init__(self, input_size: int):
        super().__init__(input_size, 0)

    def __call__(self, v: F2n) -> F2n:
        return F2n(0, 0)

    def get_ANF_propagation_model(self, input_key_mask, output_key_mask):
        return tuple((-x, ) for x in range(1, self.input_size+1))