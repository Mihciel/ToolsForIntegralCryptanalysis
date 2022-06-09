from functools import reduce
from threading import Lock
from pysat.solvers import Glucose4 as Solver
from multiprocessing.managers import SyncManager, ValueProxy
from itertools import product
from pqdm.processes import pqdm
from time import time
from galois import GF2

from .Function import Function
from .F2n import F2n
from .CipherConstructor import CipherConstructor
from .SATSolverTools import enum_projected_models, xor_clauses

def contains_trail(f: Function, u: F2n, v: F2n) -> bool:
    # returns whether a trail through f from u to v exists
    assumptions = []
    for i in range(f.input_size):
        assumptions.append((-1+2*u[i])*f.get_input_vars()[i])
    for i in range(f.output_size):
        assumptions.append((-1+2*v[i])*f.get_output_vars()[i])
    with Solver(bootstrap_with=f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(f.output_size, 0))) as solver:
        res = solver.solve(assumptions=assumptions)
    return res

def contains_trail_bulk(f: Function, us: list[F2n], vs: list[F2n], nb_threads = 1):
    f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(f.output_size, 0))
    res = pqdm([(f,) + uv for uv in zip(us, vs)],
               contains_trail, n_jobs=nb_threads, argument_type='args', desc=f"Checking trails in bulk")
    return tuple(res)


def get_balanced_output_bits(f: Function, u: F2n, nb_threads: int = 1) -> tuple[int, ...]:
    # returns the index of all outputs bits that are balanced for input monomial u
    # takes as an optional parameter for the number of threads it can use and for the model solver to use
    f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(f.output_size, 0))
    res = pqdm([(f, u, F2n(f.output_size, 2**i)) for i in range(f.output_size)],
               contains_trail, n_jobs=nb_threads, argument_type='args', desc=f"Checking balanced bits for {u}")
    return tuple(x[0] for x in enumerate(res) if not x[1])

def get_all_properties_parity_set(f: Function, u: F2n, timeout: int =2**32,  nb_threads: int = 1):
    start = time()
    Vi = list()
    Vi.append(set(F2n(f.input_size, 2**i) for i in get_balanced_output_bits(f,u, nb_threads)))
    while len(Vi[-1]) > 0 and time() - start < timeout:
        Vi.append(set())
        to_test = set()
        for v1 in Vi[-2]:
            for v2 in Vi[0]:
                v = v1 | v2
                if v.hamming_weight() == len(Vi):
                    flag = True
                    for i, b in enumerate(v):
                        if b==1 and not (v^F2n(v.n, 2**i)) in Vi[-2]:
                            flag = False
                            break
                    if flag:
                        to_test.add(v)
        to_test = list(to_test)
        for v, trail in zip(to_test, contains_trail_bulk(f, [u]*len(to_test), to_test, nb_threads)):
            if not trail:
                Vi[-1].add(v)
    return Vi

def has_balanced_output_bit(f: Function, u: F2n) -> tuple[bool, int]:
    for i in range(f.output_size):
        r = contains_trail(f, u, F2n(f.output_size, 2**i))
        if not r:
            return True, i
    return False, -1

def lowest_data_input_properties(f: Function, U=[], nb_threads: int = 1):
    balanced_bits = {}
    if len(U) == 0:
        U = [set()]
        for i in range(f.input_size):
            u = F2n(f.input_size, 2**f.input_size-2**i-1)
            balanced_bits[u] = set(get_balanced_output_bits(f, u, nb_threads=nb_threads))
            if len(balanced_bits[u]) > 0:
                U[-1].add(u)
    else:
        for u in U[-1]:
            balanced_bits[u] = set(get_balanced_output_bits(f, u, nb_threads=nb_threads))
        
    while len(U[-1]) > 0:
        U.append(set())
        to_test = []
        for u1 in U[-2]:
            for u2 in U[0]:
                u = u1 & u2
                bits_to_check = set(range(f.input_size))
                if u.hamming_weight() == f.input_size - len(U):
                    for i, b in enumerate(u):
                        if b==0:
                            bits_to_check &= balanced_bits.get(u^F2n(u.n, 2**i), set())
                    if len(bits_to_check) > 0:
                        to_test += list((u, i) for i in bits_to_check)
                        balanced_bits[u] = set()
        for (u, i), trail in zip(to_test, contains_trail_bulk(f, list(map(lambda x: x[0], to_test)), list(map(lambda x: F2n(f.output_size, 2**x[1]), to_test)), nb_threads)):
            if not trail:
                balanced_bits[u].add(i)
                U[-1].add(u)
    return U[-2]

def lowest_data_input_properties_var(f: Function, iters: int, nb_threads: int = 1):
    startmon = [F2n(f.input_size, 2**f.input_size-2**i-1) for i in range(f.input_size)]
    U = [set(map(lambda y: y[0], filter(lambda x: x[1][0], zip(startmon, pqdm([(f, u) for u in startmon], has_balanced_output_bit, n_jobs=nb_threads, argument_type='args', desc=f"checking existence balanced bit")))))]
    it = 1
    while len(U[-1]) > 0 and it < iters:
        it += 1
        to_test = []
        for u1 in U[-1]:
            for u2 in U[0]:
                u = u1 & u2
                if u.hamming_weight() == f.input_size - len(U)-1:
                    flag = True
                    for i, b in enumerate(u):
                        if b==0 and not (u^F2n(u.n, 2**i)) in U[-1]:
                            flag = False
                            break
                    if flag:
                        to_test.append(u)
        U.append(set(map(lambda y: y[0], filter(lambda x: x[1][0], zip(to_test, pqdm([(f, u) for u in to_test], has_balanced_output_bit, n_jobs=nb_threads, argument_type='args', desc=f"checking existence balanced bit"))))))
    if len(U[-1]) == 0:
        return U[-2]
    else:
        return lowest_data_input_properties(f, U, nb_threads=nb_threads)

def __get_highest_round_given_input(constructor, u: F2n, global_highest_round: ValueProxy[int], global_lock: Lock):
    local_highest_round_out: int = -1
    local_highest_round = global_highest_round.get()
    res = has_balanced_output_bit(
        constructor(local_highest_round+1), u)
    while res[0]:
        global_lock.acquire(blocking=True)
        if global_highest_round.get() < local_highest_round+1:
            global_highest_round.set(local_highest_round+1)
        global_lock.release()
        local_highest_round_out = res[1]
        local_highest_round = global_highest_round.get()
        res = has_balanced_output_bit(constructor(local_highest_round+1), u)
    return local_highest_round, u, local_highest_round_out


def get_highest_round_balanced_bit(constructor: CipherConstructor, nb_threads: int = 1) -> tuple[int, F2n, int]:
    with SyncManager() as manager:
        global_highest_round = manager.Value(int, 0)
        global_lock = manager.Lock()
        res = pqdm([(constructor, F2n(constructor.get_input_size(), 2**constructor.get_input_size() - 2**i - 1), global_highest_round, global_lock) for i in range(
            constructor.get_input_size())], __get_highest_round_given_input, n_jobs=nb_threads, argument_type='args', desc="Searching for highest round with balanced bit")
        highest_achieved_round = global_highest_round.get()
    for local_highest_round, u, local_highest_round_out in res:
        if local_highest_round == highest_achieved_round and local_highest_round_out != -1:
            return local_highest_round, u, local_highest_round_out
    # this should not happen
    return 0, F2n(constructor.get_input_size(), 2**constructor.get_input_size()-1), -1


def get_exact_propagation_given_input_monomials(f: Function, u: F2n):
    res = []
    invars = f.get_input_vars()
    outvars = f.get_output_vars()
    constraints = f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(
        f.output_size, 0))+tuple(((-1+2*u[i])*invars[i], ) for i in range(f.input_size))
    for m in enum_projected_models(constraints, outvars):
        v = 0
        for ov in outvars[::-1]:
            v <<= 1
            if m[ov-1] > 0:
                v += 1
        v = F2n(f.output_size, v)
        with Solver(bootstrap_with=constraints+tuple(((-1+2*v[i])*outvars[i], ) for i in range(f.output_size))) as solver:
            if sum(1 for _ in solver.enum_models()) % 2 == 1:
                res.append(v)
    return res

def key_add(x, k, permutation = False):
    assert(x.n == k.n)
    r = set((x, ))
    for i in range(x.n):
        if x[i]==0 and k[i] == 1:
            new_r = set()
            for v in r:
                new_r.add(v|F2n(x.n, 2**i))
            r.update(new_r)
    if permutation:
        r.discard(F2n(x.n, x.n**2-1))
    return frozenset(r)

def key_add_set(s, k, permutation):
    return reduce(frozenset.symmetric_difference, map(lambda v: key_add(v, k, permutation), s))

def branching_paths(s: frozenset[F2n], recipe, permutation=False, nb_threads: int = 1):
    S = {s}
    i = 0
    for key_mask, f in recipe:
        i+=1
        # branch paths with key
        Snew = set()
        for j, s in enumerate(S):
            Snew.update(pqdm([(s, k, permutation) for k in (key_mask & (~reduce(F2n.__and__, s))).supp()], key_add_set, n_jobs=nb_threads, argument_type='args', desc=f'branching paths for round {i} and set {j+1} of {len(S)}'))
        # propagate through f
        if f is not None:
            f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(f.output_size, 0))
            sint = reduce(frozenset.intersection, Snew)
            base = reduce(lambda x, y: frozenset.symmetric_difference(x, frozenset(y)), pqdm([(f, v) for v in sint], get_exact_propagation_given_input_monomials, n_jobs=nb_threads, argument_type='args', desc=f'base propagation for round {i}'), frozenset())
            srest = list(reduce(frozenset.union, Snew) - sint)
            prop_elements = {u:frozenset(t) for u, t in zip(srest, pqdm([(f, v) for v in srest], get_exact_propagation_given_input_monomials, n_jobs=nb_threads, argument_type='args', desc=f'other propagation for round {i}'))}
            S = set(map(lambda s: reduce(frozenset.symmetric_difference, map(prop_elements.get, s - sint), base), Snew))
    return reduce(frozenset.union, S), reduce(frozenset.intersection, S)

def get_exact_propagation_given_output_monomials(f: Function, v: F2n):
    res = []
    invars = f.get_input_vars()
    outvars = f.get_output_vars()
    constraints = f.get_ANF_propagation_model(F2n(f.input_size, 0), F2n(
        f.output_size, 0))+tuple(((-1+2*v[i])*outvars[i], ) for i in range(f.output_size))
    for m in enum_projected_models(constraints, invars):
        u = 0
        for inv in invars[::-1]:
            u <<= 1
            if m[inv-1] > 0:
                u += 1
        u = F2n(f.input_size, u)
        with Solver(bootstrap_with=constraints+tuple(((-1+2*u[i])*invars[i], ) for i in range(f.input_size))) as solver:
            if sum(1 for _ in solver.enum_models()) % 2 == 1:
                res.append(u)
    return res


def get_linear_combination_CNF(constructor: CipherConstructor, rounds: int,  input_monomials: list[F2n], output_monomials: list[F2n], nb_threads: int = 1):
    # todo add output 0x000000000...0000 to computation, this will be able to find things that are 1 (for odd inputs)
    r = constructor.get_round(0)
    r.get_ANF_propagation_model(F2n(r.input_size, 0), F2n(r.output_size, 0))
    prop_monomials_in = list(map(set, pqdm([(r, mon) for mon in input_monomials],
                             get_exact_propagation_given_input_monomials, n_jobs=nb_threads, argument_type='args', desc=f"Propagating input monomials")))
    inputs_to_check = list(reduce(set.union, prop_monomials_in))

    r = constructor.get_round(rounds-1)
    r.get_ANF_propagation_model(F2n(r.input_size, 0), F2n(r.output_size, 0))
    prop_monomials_out = list(map(set, pqdm([(constructor.get_round(rounds-1), mon) for mon in output_monomials], get_exact_propagation_given_output_monomials,
                              n_jobs=nb_threads, argument_type='args', desc=f"Back propagating output monomials through round {rounds}")))
    outputs_to_check = list(reduce(set.union, prop_monomials_out))
    middle = constructor.build(1, rounds-1)
    middle.get_ANF_propagation_model(
        F2n(middle.input_size, 0), F2n(middle.output_size, 0))

    propin2propout = list(zip(pqdm([(middle, u, v) for u, v in product(inputs_to_check, outputs_to_check)], contains_trail, n_jobs=nb_threads,
                                   argument_type='args', desc=f"searching for properties of {rounds} rounds"), product(range(len(inputs_to_check)), range(len(outputs_to_check)))))

    vars_used = 0
    imonvars = list(range(vars_used+1, vars_used + len(input_monomials)+1))
    vars_used += len(input_monomials)
    isetvars = [0]*len(inputs_to_check)
    omonvars = list(range(vars_used+1, vars_used + len(output_monomials) +1))
    vars_used += len(output_monomials)
    osetvars = [0]*len(outputs_to_check)

    clauses = []

    # add summation constraints
    for i, prop_mon in enumerate(inputs_to_check):
        xor_vars = []
        for j in range(len(input_monomials)):
            if prop_mon in prop_monomials_in[j]:
                xor_vars.append(imonvars[j])
        if len(xor_vars) == 1:
            isetvars[i] = xor_vars[0]
        else:
            tvars = xor_vars[:1]
            xor_vars = xor_vars[1:]
            for j in range(len(xor_vars)):
                tvars.append(vars_used+1)
                vars_used += 1
                replacements = [
                    0, xor_vars[j], tvars[j], tvars[j+1], -tvars[j+1], -tvars[j], -xor_vars[j]]
                for clause in xor_clauses:
                    clauses.append(tuple(replacements[x] for x in clause))
            isetvars[i] = tvars[-1]
    for i, prop_mon in enumerate(outputs_to_check):
        xor_vars = []
        for j in range(len(output_monomials)):
            if prop_mon in prop_monomials_out[j]:
                xor_vars.append(omonvars[j])
        if len(xor_vars) == 1:
            osetvars[i] = xor_vars[0]
        else:
            tvars = xor_vars[:1]
            xor_vars = xor_vars[1:]
            for j in range(len(xor_vars)):
                tvars.append(vars_used+1)
                vars_used += 1
                replacements = [
                    0, xor_vars[j], tvars[j], tvars[j+1], -tvars[j+1], -tvars[j], -xor_vars[j]]
                for clause in xor_clauses:
                    clauses.append(tuple(replacements[x] for x in clause))
            osetvars[i] = tvars[-1]
    # add propagation constraints
    for c, (i, j) in propin2propout:
        if c:
            clauses.append((-isetvars[i], -osetvars[j]))
    # make sure the all zero solution can't exist
    clauses.append(tuple(imonvars))
    clauses.append(tuple(omonvars))

    return clauses, imonvars, omonvars


def find_highest_round_with_property(constructor: CipherConstructor, extra_inmons = tuple(), extra_outmons = tuple(), start_with=None, nb_threads: int = 1) -> tuple[int, frozenset[F2n], frozenset[F2n]]:
    if start_with is None:
        round_to_start, inmon, outbit = get_highest_round_balanced_bit(
            constructor, nb_threads)
        res = (round_to_start, frozenset((inmon, )), frozenset(
            (F2n(constructor.get_output_size(), 2**outbit), )))
    else:
        res = (start_with[0], frozenset(start_with[1:2]), frozenset(
            (F2n(constructor.get_output_size(), 2**start_with[2]), )))

    has_solution = True
    inmons = [F2n(constructor.get_input_size(), 2**constructor.get_input_size()-2**i-1)
              for i in range(constructor.get_input_size())] + list(extra_inmons)
    outmons = [F2n(constructor.get_output_size(), 2**i)
               for i in range(constructor.get_output_size())] + list(extra_outmons)

    while has_solution:
        clauses, imonvars, omonvars = get_linear_combination_CNF(
            constructor, res[0]+1, inmons, outmons, nb_threads=nb_threads)
        with Solver(bootstrap_with=clauses) as solver:
            has_solution = solver.solve()
            if has_solution:
                model = solver.get_model()
                res = (res[0]+1, frozenset(inmons[i] for i in range(len(inmons)) if model[imonvars[i]-1]>0),
                       frozenset(outmons[i] for i in range(len(outmons)) if model[omonvars[i]-1]>0))

    return res

def has_a_property(constructor: CipherConstructor, rounds: int, extra_inmons = tuple(), extra_outmons = tuple(), nb_threads: int = 1):
    inmons = [F2n(constructor.get_input_size(), 2**constructor.get_input_size()-2**i-1)
              for i in range(constructor.get_input_size())] + list(extra_inmons)
    outmons = [F2n(constructor.get_output_size(), 2**i)
               for i in range(constructor.get_output_size())] + list(extra_outmons)
    clauses, imonvars, omonvars = get_linear_combination_CNF(
            constructor, rounds, inmons, outmons, nb_threads=nb_threads)
    with Solver(bootstrap_with=clauses) as solver:
        has_solution = solver.solve()
        if has_solution:
            model = solver.get_model()
            return (True, frozenset(inmons[i] for i in range(len(inmons)) if model[imonvars[i]-1]>0),
                    frozenset(outmons[i] for i in range(len(outmons)) if model[omonvars[i]-1]>0))
        else:
            return (False, frozenset(), frozenset)
    

def compute_all_properties_from_input(constructor: CipherConstructor, rounds: int, u: F2n, extra_outmons = tuple(), nb_threads: int = 1):
    outmons = [F2n(constructor.get_output_size(), 2**i)
               for i in range(constructor.get_output_size())] + list(extra_outmons)

    r = constructor.get_round(rounds-1)
    r.get_ANF_propagation_model(F2n(r.input_size, 0), F2n(r.output_size, 0))
    prop_monomials_out = list(map(set, pqdm([(constructor.get_round(rounds-1), mon) for mon in outmons], get_exact_propagation_given_output_monomials,
                              n_jobs=nb_threads, argument_type='args', desc=f"Back propagating output monomials through round {rounds}")))
    outputs_to_check = list(reduce(set.union, prop_monomials_out))
    trails = list(map(lambda x: x[0], filter(lambda y: y[1], zip(outputs_to_check, contains_trail_bulk(constructor(rounds-1), [u]*len(outputs_to_check), outputs_to_check, nb_threads)))))
    M = GF2.Zeros((len(trails), len(outmons)))
    for i, y in enumerate(trails):
        for j, propmons in enumerate(prop_monomials_out):
            if y in propmons:
                M[i, j] = 1
    nullspace = M.null_space().T
    res = []
    for j in range(nullspace.shape[1]):
        res.append(list())
        for i in range(nullspace.shape[0]):
            if nullspace[i, j] == 1:
                res[-1].append(outmons[i])
        
    return tuple(map(frozenset, res))

def extended_search(constructor: CipherConstructor, rounds: int, extra_inmons = tuple(), extra_outmons = tuple(), nb_threads: int = 1):
    input_monomials = [F2n(constructor.get_input_size(), 2**constructor.get_input_size()-2**i-1)
              for i in range(constructor.get_input_size())] + list(extra_inmons)
    output_monomials = [F2n(constructor.get_output_size(), 2**i)
               for i in range(constructor.get_output_size())] + list(extra_outmons)
    r = constructor.get_round(0)
    r.get_ANF_propagation_model(F2n(r.input_size, 0), F2n(r.output_size, 0))
    prop_monomials_in = list(map(set, pqdm([(r, mon) for mon in input_monomials],
                             get_exact_propagation_given_input_monomials, n_jobs=nb_threads, argument_type='args', desc=f"Propagating input monomials")))
    inputs_to_check = list(reduce(set.union, prop_monomials_in))

    r = constructor.get_round(rounds-1)
    r.get_ANF_propagation_model(F2n(r.input_size, 0), F2n(r.output_size, 0))
    prop_monomials_out = list(map(set, pqdm([(constructor.get_round(rounds-1), mon) for mon in output_monomials], get_exact_propagation_given_output_monomials,
                              n_jobs=nb_threads, argument_type='args', desc=f"Back propagating output monomials through round {rounds}")))
    outputs_to_check = list(reduce(set.union, prop_monomials_out))
    middle = constructor.build(1, rounds-1)
    middle.get_ANF_propagation_model(
        F2n(middle.input_size, 0), F2n(middle.output_size, 0))

    trails = list(map(lambda y: y[1], filter(lambda x: x[0], zip(pqdm([(middle, u, v) for u, v in product(inputs_to_check, outputs_to_check)], contains_trail, n_jobs=nb_threads,
                                   argument_type='args', desc=f"searching for properties of {rounds} rounds"), product(inputs_to_check, outputs_to_check)))))
    pairs = list(product(enumerate(input_monomials), enumerate(output_monomials)))
    M = GF2.Zeros((len(trails), len(pairs)))
    for l, ((i, u), (j, v)) in enumerate(pairs):
        prod = set(product(prop_monomials_in[i], prop_monomials_out[j]))
        for k, p in enumerate(trails):
            if p in prod:
                M[k, l] = 1
    nullspace = M.null_space().T
    res = []
    for j in range(nullspace.shape[1]):
        res.append(list())
        for i in range(nullspace.shape[0]):
            if nullspace[i, j] == 1:
                res[-1].append((pairs[i][0][1], pairs[i][1][1]))
    return res
