from .F2n import F2n
from galois import FieldArray, GF2
from math import log2
import sympy


def function_to_transition_table(f, imask: F2n, omask: F2n) -> FieldArray:
    res: FieldArray = GF2.Zeros((2**f.output_size, 2**f.input_size))
    for i in range(2**f.input_size):
        res[f(F2n(f.input_size, i)^imask)^omask, i] = 1
    return res


def mobius_transform(a: FieldArray):
    """
    Möbius transform implemented using the butterfly algorithm
    transforms the array in place
    """
    for i in range(int(log2(a.shape[1]))):
        for b in range(0, a.shape[1], 2**(i+1)):
            for x in range(b, b+2**i):
                a[:, x+2**i] += a[:, x]


def minterms_to_cnf_QM(minterms, nb_vars) -> tuple[tuple[int, ...], ...]:
    # this uses sympy's logic reduction using Quine-McCluskey. This is slow if the sum of inputs and outputs if above 12.
    # get cnf variables
    variables = sympy.symbols(
        ("x{} "*(nb_vars)).format(*list(range(1, nb_vars+1))))

    # build cnf
    str_clauses = str(sympy.POSform(variables, minterms)).replace("x", "").replace(
        " ", "").replace("~", "-").replace("(", "").replace(")", "").split("&")
    clauses: list[tuple[int]] = list()
    for clause in str_clauses:
        clauses.append(tuple(int(x) for x in clause.split("|")))
    return tuple(clauses)


def minterms_to_cnf(minterms, nb_vars) -> tuple[tuple[int, ...], ...]:
    return minterms_to_cnf_QM(minterms, nb_vars)

def maxterms_to_cnf(maxterms, nb_vars) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple((1-2*maxterm[i])*(i+1) for i in range(nb_vars)) for maxterm in maxterms)

def full_state_post_key_addition(a: FieldArray):
    for i in range(int(log2(a.shape[0]))):
        for b in range(0, a.shape[0], 2**(i+1)):
            for x in range(b, b+2**i):
                a[x+2**i, :] |= a[x, :]

def ANF_matrix(f):
    res = function_to_transition_table(f, F2n(f.input_size, 0), F2n(f.output_size, 0))
    mobius_transform(res[::-1, :].T)
    mobius_transform(res)
    return res

def compute_ANF_propagation_model(f, input_key_mask: F2n, output_key_mask: F2n) -> tuple[tuple[int, ...], ...]:
    """
    Computes the ANF propagation model for the given function. 
    """
    if input_key_mask == F2n(f.input_size, 2**f.input_size-1):
        if output_key_mask == F2n(f.output_size, 2**f.output_size-1):
            res = function_to_transition_table(f, F2n(f.input_size, 0), F2n(f.output_size, 0))
            mobius_transform(res[::-1, :].T)
            mobius_transform(res)
            full_state_post_key_addition(res[:, ::-1].T)
            full_state_post_key_addition(res)
        else:
            res: FieldArray = GF2.Zeros((2**f.output_size, 2**f.input_size))
            for j in output_key_mask.supp():
                xanf = function_to_transition_table(f, F2n(f.input_size, 0), j)
                mobius_transform(xanf[::-1, :].T)
                mobius_transform(xanf)
                res |= xanf
            full_state_post_key_addition(res[:, ::-1].T)
    elif output_key_mask == F2n(f.output_size, 2**f.output_size-1):
        res: FieldArray = GF2.Zeros((2**f.output_size, 2**f.input_size))
        for i in input_key_mask.supp():
            xanf = function_to_transition_table(f, i, F2n(f.output_size, 0))
            mobius_transform(xanf[::-1, :].T)
            mobius_transform(xanf)
            res |= xanf
        full_state_post_key_addition(res)
    else:
        res: FieldArray = GF2.Zeros((2**f.output_size, 2**f.input_size))
        for i in input_key_mask.supp():
            for j in output_key_mask.supp():
                xanf = function_to_transition_table(f, i, j)
                mobius_transform(xanf[::-1, :].T)
                mobius_transform(xanf)
                res |= xanf

    if f.input_size+f.output_size < 100:
        # extract all coefficients that are 1
        minterms: list[list[int]] = []
        for i in range(res.shape[0]):
            for j in range(res.shape[1]):
                if res[i, j] == 1:
                    minterms.append(
                        list(F2n(f.input_size, j).concatenate(F2n(f.output_size, i))))

        return minterms_to_cnf(minterms, f.input_size+f.output_size)
    else:
        # extract all coefficients that are 0
        maxterms: list[list[int]] = []
        for i in range(res.shape[0]):
            for j in range(res.shape[1]):
                if res[i, j] == 0:
                    maxterms.append(
                        list(F2n(f.input_size, j).concatenate(F2n(f.output_size, i))))
        return maxterms_to_cnf(maxterms, f.input_size+f.output_size)
