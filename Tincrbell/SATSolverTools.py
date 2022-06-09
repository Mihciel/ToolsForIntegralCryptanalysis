from pysat.solvers import Glucose4 as Solver
from galois import GF2

xor_clauses = ((1, 2, -3), (1, -2,  3), (-1, 2, 3), (-1, -2, -3)) # x1+x2=x3

def gen_projected_models(clauses, proj_variables, assumptions = []):
    with Solver(bootstrap_with=clauses, incr=True) as solver:
        while solver.solve(assumptions=assumptions):
            model = solver.get_model()
            solver.add_clause(tuple(-model[v-1] for v in proj_variables))
            yield model

def enum_projected_models(clauses, proj_variables, assumptions = []):
    res = []
    with Solver(bootstrap_with=clauses, incr=True) as solver:
        while solver.solve(assumptions=assumptions):
            model = solver.get_model()
            solver.add_clause(tuple(-model[v-1] for v in proj_variables))
            res.append(model)
    return res

def enum_projected_linearly_independend_models(clauses, proj_variables, assumptions=[]):
    vars_used: int = max(map(lambda x: max(map(abs, x)), clauses))
    M = GF2.Zeros((len(proj_variables), len(proj_variables)))

    nbsols = 0

    has_solution = False
    
    with Solver(bootstrap_with=clauses) as solver:
        solver.add_clause(proj_variables)
        if solver.solve(assumptions=assumptions):
            model = solver.get_model()
            for i, v in enumerate(proj_variables):
                if model[v-1] > 0:
                    M[nbsols, i] = 1
            nbsols += 1
            has_solution = True

    while has_solution:
        has_solution = False
        local_vars_used = vars_used
        # add constraints such that the result is contained in the right null space of M
        N = M.null_space()
        null_space_vars = tuple(range(local_vars_used+1, local_vars_used+N.shape[0]+1))
        local_vars_used += N.shape[0]
        with Solver(bootstrap_with=clauses) as solver:
            solver.add_clause(proj_variables)
            for j, v in enumerate(proj_variables):
                xor_vars = [v,]
                for i in range(N.shape[0]):
                    if N[i, j] == 1:
                        xor_vars.append(null_space_vars[i])
                tvars = xor_vars[:1]
                xor_vars = xor_vars[1:]
                if len(xor_vars) > 0:
                    for i in range(len(xor_vars)):
                        tvars.append(local_vars_used+1)
                        local_vars_used += 1
                        replacements = [
                            0, xor_vars[i], tvars[i], tvars[i+1], -tvars[i+1], -tvars[i], -xor_vars[i]]
                        for clause in xor_clauses:
                            solver.add_clause(tuple(replacements[x] for x in clause))
                solver.add_clause((-tvars[-1], ))
            if solver.solve(assumptions=assumptions):
                model = solver.get_model()
                for i, v in enumerate(proj_variables):
                    if model[v-1] > 0:
                        M[nbsols, i] = 1
                        nbsols += 1
                has_solution = True
    
    reducedM = M.row_space()
    # convert matrix back to solutions
    res = []
    for i in range(reducedM.shape[0]):
        res.append([None]*vars_used)
        for j, v in enumerate(proj_variables):
            if reducedM[i, j] == 1:
                res[-1][v-1] = v
            else:
                res[-1][v-1] = -v
    return res
    