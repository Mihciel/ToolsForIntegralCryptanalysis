import sys
if not 'Tincrbell' in sys.modules:
    sys.path.append(".")
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from Tincrbell.F2n import F2n
from Tincrbell.Components import SBox
from Tincrbell.CipherConstructor import CipherConstructor
from Tincrbell.SearchAlgorithms import get_highest_round_balanced_bit, find_highest_round_with_property, get_balanced_output_bits, get_all_properties_parity_set, has_balanced_output_bit
from Tincrbell.Linear import LinearFunction, linear_to_compound
from Tincrbell.Function import Function
from galois import GF, GF2
from numpy import matmul
from time import time
from functools import reduce


RIJNDAEL_sbox = SBox(8, 8, 
        [
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
        ])
RIJNDAEL128_ShiftRows = [8*(((i//8)-4*((i//8) % 4)) % 16)+(i % 8) for i in range(128)]
RIJNDAEL128_MixColumns_M256 = GF(2**8)([[2, 3, 1, 1], [1, 2, 3, 1], [1, 1, 2, 3], [3, 1, 1, 2]])

# https://crypto.stackexchange.com/questions/33085/how-to-perform-aes-mixcolumns-as-matrix-multiplication-in-gf2-boolean-values
RIJNDAEL128_MixColumns_M2 = GF2.Zeros((32, 32))
for i in range(4):
    for j in range(8):
        x = GF(2**8).Zeros((4, 1))
        x[i, 0] = GF(2**8)(1<<j)
        v = matmul(RIJNDAEL128_MixColumns_M256, x)
        for a in range(4):
            for b in range(8):
                RIJNDAEL128_MixColumns_M2[8*a+b, 8*i+j] = (int(v[a])>>b)&1
RIJNDAEL128_MixColumns = linear_to_compound(LinearFunction(RIJNDAEL128_MixColumns_M2))

RIJNDAEL128_SubBytes = CompoundFunction(128, 128)
for i in range(128//8):
    ids = RIJNDAEL128_SubBytes.add_component(RIJNDAEL_sbox)
    for j in range(8):
        RIJNDAEL128_SubBytes.connect_components(INPUT_ID, 8*i+j, ids, j)
        RIJNDAEL128_SubBytes.connect_components(ids, j, OUTPUT_ID, 8*i+j)

RIJNDAEL128_roundfunction = CompoundFunction(128, 128)
sbox_ids = []
for i in range(128//8):
    sbox_ids.append(RIJNDAEL128_roundfunction.add_component(RIJNDAEL_sbox))
    for j in range(8):
        RIJNDAEL128_roundfunction.connect_components(INPUT_ID, i*8+j, sbox_ids[-1], j)
lbox_ids = []
for i in range(128//32):
    lbox_ids.append(RIJNDAEL128_roundfunction.add_component(RIJNDAEL128_MixColumns))
    RIJNDAEL128_roundfunction.set_input_key_mask(lbox_ids[-1], F2n(32, 2**32-1))
    for j in range(32):
        RIJNDAEL128_roundfunction.connect_components(lbox_ids[-1], j, OUTPUT_ID, i*32 + j)
for i in range(128):
    RIJNDAEL128_roundfunction.connect_components(sbox_ids[i//8], i%8, lbox_ids[RIJNDAEL128_ShiftRows[i]//32], RIJNDAEL128_ShiftRows[i] % 32)

class Rijndael128Constructor():
    # constructs a block cipher by using the given round functions
    # the key masks are seperated from the round functions
    # when more rounds are constructed than there are different round_functions give,
    # it will reuse the first round functions

    def __call__(self, n: int) -> Function:
        assert(n > 0)
        f = CompoundFunction(128, 128)
        idc = f.add_component(RIJNDAEL128_roundfunction)
        f.set_input_key_mask(idc, F2n(128, 2**128-1))
        for i in range(128):
            f.connect_components(INPUT_ID, i, idc, i)
        for i in range(1, n):
            idp = idc
            if i == n-1:
                idc = f.add_component(RIJNDAEL128_SubBytes)
            else:
                idc = f.add_component(RIJNDAEL128_roundfunction)
            f.set_input_key_mask(idc, F2n(128, 2**128-1))
            for j in range(128):
                f.connect_components(idp, j, idc, j)
        for i in range(128):
            f.connect_components(idc, i, OUTPUT_ID, i)
        f.set_input_key_mask(OUTPUT_ID, F2n(128, 2**128-1))
        return f

    def get_input_size(self) -> int:
        return 128

    def get_output_size(self) -> int:
        return 128

RIJNDAEL128_builder = Rijndael128Constructor()

if __name__=="__main__":

    if len(sys.argv) > 1:
        nb_threads = int(sys.argv[1])
    else:
        nb_threads = int(input("Please input the number of threads to use: "))

    start = time()
    m = RIJNDAEL_sbox.get_ANF_propagation_model(F2n(8, 0xff), F2n(8, 0xff))
    print(f"{len(m)} constraints in {time()-start}s")

    # verify literature
    # start = time()
    # res = get_balanced_output_bits(RIJNDAEL128_builder(3), F2n(128, 0xff), nb_threads=nb_threads)
    # print(f"Daemen and Rijmen 3 round 2^8 data (1998): {len(res) == 128} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(RIJNDAEL128_builder(4), F2n(128, 0xff00000000ff00000000ff00000000ff), nb_threads=nb_threads)
    # print(f"Daemen and Rijmen 4 round 2^32 data (1998): {len(res) == 128} in {time()-start}s")

    # improving literature
    # start = time()
    # res = get_all_properties_parity_set(RIJNDAEL128_builder(3), F2n(128, 0xff), timeout=120, nb_threads=nb_threads)
    # print(f"Daemen and Rijmen 3 round 2^8 data (1998): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(RIJNDAEL128_builder(4), F2n(128, 0xff00000000ff00000000ff00000000ff), timeout=120, nb_threads=nb_threads)
    # print(f"Daemen and Rijmen 4 round 2^32 data (1998): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")

    # res = get_highest_round_balanced_bit(RIJNDAEL128_builder, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]}")

    # trying to find lower data distinguisher on 3 round AES
    inputs = [F2n(128, 0xff)^F2n(128, 1<<i) for i in range(8)]
    f = RIJNDAEL128_builder(3)
    f.get_ANF_propagation_model(F2n(128, 0), F2n(128, 0))
    for u in inputs:
        print(u, has_balanced_output_bit(f, u))