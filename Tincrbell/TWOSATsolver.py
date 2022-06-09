from typing import Collection
from networkx import Digraph

def solve2SAT(clauses: Collection[tuple[int, int]], assumptions: Collection[int]) -> tuple[bool, list[int]]:
    nvars = max(map(lambda x: max(abs(x), abs(y)), clauses))
    res = [None]*(nvars+1)
    assumptions = set(assumptions)
    for a in assumptions:
        if -a in assumptions:
            return False,  []

    # propagate assumptions
    while len(assumptions):
        new_clauses = []
        a = assumptions.pop()
        if a > 0:
            res[a] = True
        else:
            res[-a] = False
        for clause in clauses:
            if clause[0] == -a:
                if -clause[1] in assumptions:
                    return False, []
                else:
                    assumptions.add(clause[1])
            elif clause[1] == -a:
                if -clause[0] in assumptions:
                    return False, []
                else:
                    assumptions.add(clause[0])
            elif a not in clause:
                new_clauses.append(a)
        clauses = new_clauses

    # construct implication graph
    edges = set()
    for clause in clauses:
        edges.add((-clause[0], clause[1]))
        edges.add((-clause[1], clause[0]))
    G = Digraph(edges)

    # if graph has cycles, no solutions exist
    if not G.is_directed_acyclic_graph():
        return False, []

    # get solution
    inv_top_sort = {x: i for i, x in enumerate(G.topological_sort())}
    for i in range(1, nvars+1):
        if res[i] is None:
            if inv_top_sort[i] < inv_top_sort[-i]:
                res[i] = False
            else:
                res[i] = True
    return True, res