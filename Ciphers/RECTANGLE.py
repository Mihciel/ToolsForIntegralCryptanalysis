import sys
if not 'Tincrbell' in sys.modules:
    sys.path.append(".")
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from Tincrbell.F2n import F2n
from Tincrbell.Components import SBox
from Tincrbell.CipherConstructor import CipherConstructor
from Tincrbell.SearchAlgorithms import find_highest_round_with_property
from Tincrbell.SearchAlgorithms import get_highest_round_balanced_bit, find_highest_round_with_property

RECTANGLE_sbox = SBox(4, 4, [6, 5, 0xC, 0xA, 1, 0xE, 7, 9, 0xB, 0, 3, 0xD, 8, 0xF, 4, 2])
RECTANGLE_bitpermutation = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 44, 45, 46, 47, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 61, 62, 63, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]
RECTANGLE_SubColumn = CompoundFunction(64, 64)
for i in range(16):
    sid =  RECTANGLE_SubColumn.add_component(RECTANGLE_sbox)
    for j in range(4):
        RECTANGLE_SubColumn.connect_components(INPUT_ID, 16*j+i, sid, j)
        RECTANGLE_SubColumn.connect_components(sid, j, OUTPUT_ID, 16*j+i)

RECTANGLE_roundfunction = CompoundFunction(64, 64)
subid = RECTANGLE_roundfunction.add_component(RECTANGLE_SubColumn)
for i in range(64):
    RECTANGLE_roundfunction.connect_components(INPUT_ID, i, subid, i)
    RECTANGLE_roundfunction.connect_components(subid, i, OUTPUT_ID, RECTANGLE_bitpermutation[i])
RECTANGLE_builder = CipherConstructor(RECTANGLE_roundfunction, F2n(64, 2**64-1))

if __name__=="__main__":

    if len(sys.argv) > 1:
        nb_threads = int(sys.argv[1])
    else:
        nb_threads = int(input("Please input the number of threads to use: "))

    res = get_highest_round_balanced_bit(RECTANGLE_builder, nb_threads=nb_threads)
    print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]}")
    # rounds: 9
    # input: F2n(64, 0xfffffffffffffffe)
    # output bit: 0

    res = find_highest_round_with_property(RECTANGLE_builder, start_with=res, nb_threads=nb_threads)
    print(f"rounds: {res[0]}\ninput set: {res[1]}\noutput set: {res[2]}")
    # rounds: 10
    # input set: frozenset({F2n(64, 0xfffffdffffffffff), F2n(64, 0xfffffffffdffffff)})
    # output set: frozenset({F2n(64, 0x0000000000004000)})