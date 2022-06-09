from .Function import Function
from .F2n import F2n
from .Components import COPYn
from random import Random

INPUT_ID: int = 0
OUTPUT_ID: int = 2**64
ERROR_ID: int = -1


class Null(Function):
    def __init__(self, input_size: int, output_size: int):
        super().__init__(input_size, output_size)
        self.n_cnf_vars = 0

    def __call__(self, v: F2n) -> F2n:
        return F2n(self.output_size, 0)

    def compute_ANF_propagation_model(self):
        return tuple()

class component_record(object):
    def __init__(self, f: Function):
        self.f = f
        self.input_connections: list[tuple[int, int]] = [(ERROR_ID, ERROR_ID)]*f.input_size
        self.output_connections: list[list[tuple[int, int]]] = [list() for _ in range(f.output_size)]
        self.input_vars: list[int] = [0]*f.input_size
        self.output_vars: list[int] = [0]*f.output_size
        self.copy_output_vars: list[list[int]] = [list() for _ in range(f.output_size)]
        self.input_key_mask = F2n(f.input_size, 0)
        self.output_key_mask = F2n(f.output_size, 0)


class CompoundFunction(Function):
    def __init__(self, input_size: int, output_size: int):
        super().__init__(input_size, output_size)
        # component list contains component, inputs, outputs, input_vars, output_vars
        self.components: dict[int, component_record] =\
            {INPUT_ID: component_record(Null(0, self.input_size)) ,
             OUTPUT_ID: component_record(Null(self.output_size, 0))}
        self.used_vars: int = 0
        self.n_components: int = 0
        self.assigned_vars: bool = False
        self.n_cnf_vars = 0
        self.seed = 0

    def add_component(self, component: Function) -> int:
        """
        Add component to Compound Function.
        Components should be added in order of execution.
        The id of the component is returned
        """
        self.n_components += 1
        self.components[self.n_components] = component_record(component)
        self.n_cnf_vars += component.get_number_cnf_vars() - \
            component.input_size - component.output_size
        return self.n_components

    def connect_components(self, from_component: int, from_wire: int, to_component: int, to_wire: int) -> None:
        """
        connect one component to another
        """
        # check existence of components
        assert(-1 < from_component <= self.n_components)
        assert(-1 < to_component <= self.n_components or to_component == OUTPUT_ID)
        # check existence of wires
        assert(from_wire < self.components[from_component].f.output_size)
        assert(to_wire < self.components[to_component].f.input_size)
        # enforce order of computation
        assert(from_component < to_component)
        # update information
        self.components[from_component].output_connections[from_wire].append((to_component, to_wire))
        assert(self.components[to_component].input_connections[to_wire] == (ERROR_ID, ERROR_ID))
        self.components[to_component].input_connections[to_wire] = (from_component, from_wire)

    def assign_vars(self) -> None:
        # check if everything is connected
        for record in self.components.values():
            assert((ERROR_ID, ERROR_ID) not in record.input_connections)
            assert(0 not in map(len, record.output_connections))
        # set variables
        for record in self.components.values():
            if not isinstance(record.f, Null):
                assert((ERROR_ID, ERROR_ID) not in record.input_connections)
                assert(0 not in map(len, record.output_connections))
        for id in range(self.n_components+1):
            if id == INPUT_ID or not isinstance(self.components[id].f, Null):
                for w in range(self.components[id].f.output_size):
                    self.n_cnf_vars += 1
                    self.used_vars += 1
                    self.components[id].output_vars[w] = self.used_vars
                    if len(self.components[id].output_connections[w]) > 1:
                        for c, wi in self.components[id].output_connections[w]:
                            self.n_cnf_vars += 1
                            self.used_vars += 1
                            self.components[id].copy_output_vars[w].append(self.used_vars)
                            self.components[c].input_vars[wi] = self.used_vars
                    else:
                        c, wi = self.components[id].output_connections[w][0]
                        self.components[c].input_vars[wi] = self.used_vars
        for record in self.components.values():
            assert(0 not in record.input_vars or  isinstance(record.f, Null))
            assert(0 not in record.output_vars or  isinstance(record.f, Null))
        # fix output masks
        for id in range(1, self.n_components+1):
            for w in range(self.components[id].f.output_size):
                if all(map(lambda x: self.components[x[0]].input_key_mask[x[1]] == 1, self.components[id].output_connections[w])):
                    self.components[id].output_key_mask ^= F2n(self.components[id].f.output_size, 2**w)
        # disallow implicit copy between component and output
        for id in range(self.n_components+1):
            for w in range(self.components[id].f.output_size):
                assert(len(self.components[id].output_connections[w]) == 1 or OUTPUT_ID not in map(lambda x: x[0], self.components[id].output_connections[w]))
        self.assigned_vars = True


    def get_input_vars(self) -> tuple[int]:
        if not self.assigned_vars:
            self.assign_vars()
        return tuple(self.components[INPUT_ID].output_vars)

    def get_output_vars(self) -> tuple[int]:
        if not self.assigned_vars:
            self.assign_vars()
        return tuple(self.components[OUTPUT_ID].input_vars)

    def get_number_cnf_vars(self) -> int:
        if not self.assigned_vars:
            self.assign_vars()
        return self.n_cnf_vars

    def __call__(self, v: F2n) -> F2n:
        # evaluates with random key based on self.seed
        r = Random(self.seed)
        outputs: list[F2n] = [v]
        for i in range(1, self.n_components+1):
            # construct input
            x = 0
            for wire in self.components[i].input_connections[::-1]:
                x <<= 1
                x += outputs[wire[0]][wire[1]]
            # get output
            key = F2n(self.components[i].f.input_size, r.randrange(2**self.components[i].f.input_size)) & self.components[i].input_key_mask
            outputs.append(self.components.get(i, component_record(Null(0, 0))).f(F2n(self.components[i].f.input_size, x)^key))
        # contstruct output
        x = 0
        for wire in self.components[OUTPUT_ID].input_connections[::-1]:
            x <<= 1
            x += outputs[wire[0]][wire[1]]
        key = F2n(self.output_size, r.randrange(2**self.output_size)) & self.components[OUTPUT_ID].input_key_mask
        return F2n(self.output_size, x) ^ key

    def set_input_key_mask(self, component: int, mask: F2n):
        # add key addition to input wire of component
        assert(mask.n == self.components[component].f.input_size)
        self.components[component].input_key_mask = mask

    def compute_ANF_propagation_model(self, input_key_mask, output_key_mask) -> tuple[tuple[int, ...], ...]:
        if not self.assigned_vars:
            self.assign_vars()
        res: list[tuple[int]] = []
        # build copy functions
        for component in self.components.values():
            for w in range(len(component.output_vars)):
                if len(component.output_connections[w]) > 1:
                    cnf_replacement_vars = [0, component.output_vars[w]] + component.copy_output_vars[w]
                    cnf_replacement_vars += [-x for x in cnf_replacement_vars[1:][::-1]]
                    for clause in COPYn(len(component.copy_output_vars[w])).get_ANF_propagation_model(F2n(1, 0), F2n(len(component.copy_output_vars[w]), 0)):
                        res.append(tuple(cnf_replacement_vars[x] for x in clause))

        # build other functions
        for component in self.components.values():
            f = component.f
            if not isinstance(f, Null):
                input_vars = component.input_vars
                output_vars = component.output_vars
                cnf_replacement_vars = [0]*(f.get_number_cnf_vars()+1) # first element is never used
                input_vars_to_replace = f.get_input_vars()
                output_vars_to_replace = f.get_output_vars()
                for j in range(f.input_size):
                    cnf_replacement_vars[input_vars_to_replace[j]] = input_vars[j]
                for j in range(f.output_size):
                    cnf_replacement_vars[output_vars_to_replace[j]] = output_vars[j]
                for j in range(1, len(cnf_replacement_vars)):
                    if cnf_replacement_vars[j] == 0:
                        self.used_vars += 1
                        cnf_replacement_vars[j] = self.used_vars
                cnf_replacement_vars += [-x for x in cnf_replacement_vars[1:][::-1]] # for negative vars
                # modify component input/output key mask based on global input/output key mask
                cinput_key_mask = component.input_key_mask
                for i in range(f.input_size):
                    c, w = component.input_connections[i]
                    if c == INPUT_ID:
                        cinput_key_mask |= F2n(f.input_size, input_key_mask[w]*(2**i))
                coutput_key_mask = component.output_key_mask
                for i in range(f.output_size):
                    c, w = component.output_connections[i][0] # should only have one component
                    if c == OUTPUT_ID:
                        coutput_key_mask |= F2n(f.output_size, output_key_mask[w]*(2**i))
                for clause in f.get_ANF_propagation_model(input_key_mask=cinput_key_mask, output_key_mask=coutput_key_mask):
                    res.append(tuple(cnf_replacement_vars[x] for x in clause))
        return tuple(res)
