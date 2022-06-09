import sys
if not 'Tincrbell' in sys.modules:
    sys.path.append(".")
from Tincrbell.Function import Function
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from Tincrbell.Components import Sink, INVERT, XORn, AND
from Tincrbell.F2n import F2n
from Tincrbell.SearchAlgorithms import get_exact_propagation_given_output_monomials

class T(Function):
    def __init__(self):
        super().__init__(5, 1)
        
    def __call__(self, v: F2n):
        return F2n(1, (v[0] + v[1]*v[2] + v[3] + v[4]) % 2)
t = T()

z = XORn(6)

def connect_t_components(f, register, id0, inputs):
    for i in range(5):
        f.connect_components(*register[inputs[i]], id0, i)

def add_round(f, register):
    global t
    id0 = f.add_component(t)
    connect_t_components(f, register, id0, (65, 90, 91, 92, 170))
    id1 = f.add_component(t)
    connect_t_components(f, register, id1, (161, 174, 175, 176, 263))
    id2 = f.add_component(t)
    connect_t_components(f, register, id2, (242, 285, 286, 287, 68))
    register[92] = (id0,0)
    register[176] = (id1,0)
    register[287] = (id2,0)
    return register[287:] + register[:287]

def add_output(f, register):
    global z
    id0 = f.add_component(z)
    to_use = (65, 92, 161, 176, 242, 287)
    for i in range(6):
        f.connect_components(*register[to_use[i]], id0, i)
    f.connect_components(id0, 0, OUTPUT_ID, 0)
    id1 = f.add_component(Sink(288-6))
    j = 0
    for i in range(288):
        if i not in to_use:
            f.connect_components(*register[i], id1, j)
            j += 1

if __name__=="__main__":
    R = 841
    trivium = CompoundFunction(288, 1)
    register = [(INPUT_ID, x) for x in range(288)]
    # invert the last three inputs (since these are 1 and the tools only consider a zero input)
    for i in range(3):
        nid = trivium.add_component(INVERT)
        trivium.connect_components(*register[-i-1], nid, 0)
        register[-i-1] = (nid, 0)
    # build trivium
    for _ in range(R):
        register = add_round(trivium, register)
    add_output(trivium, register)
    # set fixed bit constraints
    invars = trivium.get_input_vars()
    for i in range(80, 93):
        trivium.add_correction((-invars[i],))
    for i in range(173, 288):
        trivium.add_correction((-invars[i],))
    # set cube
    I = F2n(80, 2**80 - 2**78 - 2**8 - 1)
    for i in range(80):
        trivium.add_correction(((2*I[i]-1)*invars[i+93],))

    # compute superpoly
    res = get_exact_propagation_given_output_monomials(trivium, F2n(1, 1))
    print(len(res), max(map(F2n.hamming_weight, res)))

