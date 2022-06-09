from abc import ABC, abstractmethod
from .F2n import F2n
from ._ANF_prop import compute_ANF_propagation_model


class Function(ABC):
    """
    Abstraction of a function.
    Can be called.
    Can generated or give the model needed for ANF propagation.
    """

    def __init__(self, input_size: int, output_size: int):
        self.input_size: int = input_size
        self.n_cnf_vars: int = input_size + output_size
        self.output_size: int = output_size
        self.ANF_prop_models: dict[tuple[F2n, F2n], tuple[tuple[int, ...], ...]] = {}
        self.corrections: list[tuple[int, ...]] = []

    def get_input_vars(self) -> tuple[int]:
        return tuple(range(1, self.input_size+1))

    def get_output_vars(self) -> tuple[int]:
        return tuple(range(self.input_size+1, self.output_size + self.input_size+1))

    def get_number_cnf_vars(self):
        return self.n_cnf_vars

    @abstractmethod
    def __call__(self, v: F2n):
        return F2n(0, 0)

    def compute_ANF_propagation_model(self, input_key_mask, output_key_mask) -> tuple[tuple[int, ...], ...]:
        return compute_ANF_propagation_model(self, input_key_mask, output_key_mask)

    def get_ANF_propagation_model(self, input_key_mask: F2n, output_key_mask: F2n) -> tuple[tuple[int, ...], ...]:
        assert(input_key_mask.n == self.input_size)
        assert(output_key_mask.n == self.output_size)
        if (input_key_mask, output_key_mask) not in self.ANF_prop_models:
            self.ANF_prop_models[(input_key_mask, output_key_mask)] = self.compute_ANF_propagation_model(input_key_mask, output_key_mask)
        return self.ANF_prop_models[(input_key_mask, output_key_mask)] + tuple(self.corrections)

    def add_correction(self, correction: tuple[int]):
        self.corrections.append(correction)
