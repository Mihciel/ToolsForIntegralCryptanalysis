from .Function import Function
from .CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from .F2n import F2n

class CipherConstructor():
    # constructs a block cipher by using the given round functions
    # the key masks are seperated from the round functions
    # when more rounds are constructed than there are different round_functions give,
    # it will reuse the first round functions

    def __init__(self, round_functions: Function | tuple[Function], key_masks: F2n | tuple[F2n]) -> None:
        if isinstance(round_functions, Function):
            self.round_functions: tuple[Function] = (round_functions,)
        else:
            self.round_functions: tuple[Function] = tuple(round_functions)
        if isinstance(key_masks, F2n):
            self.key_masks = (key_masks,)
        else:
            self.key_masks = tuple(key_masks)
        self.input_size = self.round_functions[0].input_size
        self.output_size = self.round_functions[0].output_size
        
    def __call__(self, n: int) -> Function:
        assert(n > 0)
        f = CompoundFunction(self.input_size, self.input_size)
        idc = f.add_component(self.get_round(0))
        f.set_input_key_mask(idc, self.get_key_mask(0))
        for i in range(self.input_size):
            f.connect_components(INPUT_ID, i, idc, i)
        for i in range(1, n):
            idp = idc
            idc = f.add_component(self.get_round(i))
            f.set_input_key_mask(idc, self.get_key_mask(i))
            for j in range(self.input_size):
                f.connect_components(idp, j, idc, j)
        for i in range(self.input_size):
            f.connect_components(idc, i, OUTPUT_ID, i)
        f.set_input_key_mask(OUTPUT_ID, self.get_key_mask(n))
        return f

    def build(self, startr, stopr):
        # build cipher from round startr to (but not including) round stopr
        assert(stopr>startr)
        assert(startr>=0)
        f = CompoundFunction(self.input_size, self.input_size)
        idc = f.add_component(self.get_round(startr))
        f.set_input_key_mask(idc, self.get_key_mask(startr))
        for i in range(self.input_size):
            f.connect_components(INPUT_ID, i, idc, i)
        for i in range(startr+1, stopr):
            idp = idc
            idc = f.add_component(self.get_round(i))
            f.set_input_key_mask(idc, self.get_key_mask(i))
            for j in range(self.input_size):
                f.connect_components(idp, j, idc, j)
        for i in range(self.input_size):
            f.connect_components(idc, i, OUTPUT_ID, i)
        f.set_input_key_mask(OUTPUT_ID, self.get_key_mask(stopr))
        return f

    def get_input_size(self) -> int:
        return self.input_size

    def get_output_size(self) -> int:
        return self.output_size

    def get_round(self, r: int) -> Function:
        return self.round_functions[r % len(self.round_functions)]

    def get_key_mask(self, r: int) -> F2n:
        return self.key_masks[r % len(self.key_masks)]