from .Function import Function
from .F2n import F2n
from .CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from .Components import COPYn, XORn
from galois import FieldArray, GF2
import numpy as np

class LinearFunction(Function):
    def __init__(self, A: FieldArray):
        self.A = A
        super().__init__(A.shape[1], A.shape[0])

    def __call__(self, v: F2n):
        x = GF2([[x] for x in v])
        y = np.matmul(self.A, x)
        res = 0
        for i in range(self.output_size):
            res <<= 1
            if y[self.output_size-1-i]:
                res += 1
        return F2n(self.output_size, res)

def linear_to_compound(f: LinearFunction):
    Lc = CompoundFunction(f.input_size, f.output_size)
    cids=[]
    cnis=[]
    for i in range(f.input_size):
        ni = 0
        for j in range(f.output_size):
            if f.A[j, i] == 1:
                ni += 1
        cids.append(Lc.add_component(COPYn(ni)))
        Lc.connect_components(INPUT_ID, i, cids[-1], 0)
        cnis.append(0)
    for i in range(f.output_size):
        ni = 0
        for j in range(f.input_size):
            if f.A[i, j] == 1:
                ni += 1
        idl = Lc.add_component(XORn(ni))
        l = 0
        for j in range(f.input_size):
            if f.A[i, j] == 1:
                Lc.connect_components(cids[j], cnis[j], idl, l)
                cnis[j]+=1
                l+=1
        Lc.connect_components(idl, 0, OUTPUT_ID, i)
    return Lc