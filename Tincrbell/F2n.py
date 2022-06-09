from typing import Generator

class F2n(object):
    """
    Boolean vectors of size n
    """

    def __init__(self, n: int, v: int) -> None:
        assert(isinstance(v, int))
        self.v: int = v
        self.n: int = n

    def __iter__(self) -> 'F2n':
        self.iter_v: int = self.v << 1
        self.iter_n: int = self.n
        return self

    def __next__(self) -> int:
        if self.iter_n == 0:
            raise StopIteration
        self.iter_n -= 1
        self.iter_v >>= 1
        return (self.iter_v & 1)

    def __invert__(self) -> 'F2n':
        return F2n(self.n, ~self.v)

    def __and__(self, other: 'F2n') -> 'F2n':
        assert(self.n == other.n)
        return F2n(self.n, self.v & other.v)

    def __or__(self,  other: 'F2n') -> 'F2n':
        assert(self.n == other.n)
        return F2n(self.n, self.v | other.v)

    def __xor__(self,  other: 'F2n') -> 'F2n':
        assert(self.n == other.n)
        return F2n(self.n, self.v ^ other.v)

    def __pow__(self, exp: 'F2n') -> int:
        assert(self.n == exp.n)
        if (self.v | ~exp.v) == -1:
            return 1
        return 0

    def __eq__(self, other: object) -> bool:
        return isinstance(other, type(self)) and self.n == other.n and self.v == other.v

    def __le__(self,  other: 'F2n') -> bool:
        assert(self.n == other.n)
        # self <= other
        # equivalent with other**self == 1
        return (other.v | (~self.v)) == -1

    def __ge__(self,  other: 'F2n') -> bool:
        # self >= other
        # equivalent with self**other == 1
        return (self.v | (~other.v)) == -1

    def __index__(self) -> int:
        return self.v & (2**self.n-1)

    def __getitem__(self, i: int) -> int:
        return (self.v >> i) & 1

    def __repr__(self) -> str:
        return f"F2n({self.n}, 0x{self.v & (2**self.n-1):0{self.n//4}x})"

    def __str__(self) -> str:
        return f"F2n({self.n}, 0x{self.v & (2**self.n-1):0{self.n//4}x})"

    def __hash__(self) -> int:
        self.v &= 2**self.n-1
        return (self.n.__hash__() + self.v.__hash__()).__hash__()

    def change_bit(self, i: int, b: int) -> 'F2n':
        res: 'F2n' = F2n(self.n, self.v)
        b &= 1
        b <<= i
        c = 1 << i
        res.v &= ~c
        res.v |= b
        return res

    def concatenate(self, other: 'F2n') -> 'F2n':
        # return other||self
        r: int = other.v << self.n
        r += self.v & (2**self.n-1)
        return F2n(self.n+other.n, r)

    def split(self, n: int) -> tuple:
        return F2n(n, self.v & (2**self.n-1)), F2n(self.n - n, self.v >> n)

    def hamming_weight(self) -> int:
        n = self.v & (2**self.n-1)
        return n.bit_count()

    def reverse(self) -> 'F2n':
        r = 0
        for i in self:
            r <<=1
            r += i
        return F2n(self.n, r)

    def supp(self) -> Generator['F2n', None, None]:
        n = self.hamming_weight()
        positions1 = []
        for i in range(self.n):
            if self[i] == 1:
                positions1.append(i)
        for i in range(2**n):
            res = F2n(self.n, 0)
            j = F2n(n, i)
            for k in range(n):
                res ^= F2n(self.n, j[k]*(2**positions1[k]))
            yield res

    def get_slice(self, position: int, size: int):
        return F2n(size, (self.v >> position) & (2**size - 1))