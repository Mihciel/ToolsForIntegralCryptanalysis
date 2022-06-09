import sys
if not 'Tincrbell' in sys.modules or not "Ciphers" in sys.modules:
    sys.path.append(".")
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID, ERROR_ID
from Tincrbell.F2n import F2n
from Tincrbell.Components import SBox
from Tincrbell.CipherConstructor import CipherConstructor
from Tincrbell.SearchAlgorithms import get_highest_round_balanced_bit, find_highest_round_with_property, get_balanced_output_bits, get_all_properties_parity_set, compute_all_properties_from_input, lowest_data_input_properties_var
from Ciphers.GIFT64 import GIFT_sbox, GIFT_modsbox, GIFT_round_constants
from time import time
from functools import reduce

GIFT128_bitpermutation = [4*(i//16) + 32*((3 * ((i%16)//4) + (i % 4))%4) + (i % 4) for i in range(128)]
GIFT128_round_functions = []
for c in GIFT_round_constants:
    GIFT128_round_functions.append(CompoundFunction(128, 128))
    for i in range(32):
        sboxid = ERROR_ID
        if (GIFT128_bitpermutation[4*i+3]//4 < 6 and c[GIFT128_bitpermutation[4*i+3]//4] == 1) or GIFT128_bitpermutation[4*i+3]//4 == 31:
            sboxid = GIFT128_round_functions[-1].add_component(GIFT_modsbox)
        else:
            sboxid = GIFT128_round_functions[-1].add_component(GIFT_sbox)
        for j in range(4):
            GIFT128_round_functions[-1].connect_components(INPUT_ID, 4*i+j, sboxid, j)
            GIFT128_round_functions[-1].connect_components(sboxid, j, OUTPUT_ID, GIFT128_bitpermutation[4*i+j])
# don't add the S and P layer of the first round, only the key addition
# since these can trivially be inverted
GIFT128_builder = CipherConstructor(tuple(GIFT128_round_functions[1:]), F2n(128, 0x33333333333333333333333333333333))

if __name__=="__main__":
    if len(sys.argv) > 1:
        nb_threads = int(sys.argv[1])
    else:
        nb_threads = int(input("Please input the number of threads to use: "))
    
    # check highest rounds
    # start = time()
    # res = get_highest_round_balanced_bit(GIFT128_builder, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]} in {time() - start}s")
    # rounds: 11
    # input: F2n(128, 0xfffffffffffffffffffffffffffffffb)
    # output bit: 0

    # res = find_highest_round_with_property(GIFT128_builder, start_with=res, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]} in {time()-start}s")
    # rounds: 11
    # input: frozenset({F2n(128, 0xfffffffffffffffffffffffffffffffb)})
    # output bit: frozenset({F2n(128, 0x00000000000000000000000000000001)})

    # confirm literature
    # start = time()
    # res = get_balanced_output_bits(GIFT128_builder(11), F2n(128, 0xffffffffffffffffffffffffffffffbf), nb_threads=nb_threads)
    # print(f"Banik et al. 12 round 2^127 data (2017): {len(res) >= 32} in {time()-start}s")

    # further analysis of literature results:
    # start = time()
    # res1 = get_all_properties_parity_set(GIFT128_builder(11), F2n(128, 0xffffffffffffffffffffffffffffffbf), nb_threads=nb_threads)
    # res2 = compute_all_properties_from_input(GIFT128_builder, 11,  F2n(128, 0xffffffffffffffffffffffffffffffbf), nb_threads=nb_threads)
    # print(f"Improving 12-round 2^127 Banik et al. (2017): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")

    # searching for highest number of properties on 12 rounds
    start = time()
    extra_inmons = tuple(F2n(4*i, 2**(4*i)-1).concatenate(F2n(4, 0x3)).concatenate(F2n(4*(31-i), 2**(4*(31-i))-1)) for i in range(32))
    extra_outmons = tuple(~x for x in extra_inmons)
    res = []
    f = GIFT128_builder(11)
    f.get_ANF_propagation_model(F2n(128, 0), F2n(128, 0))
    inputs = [F2n(128, 2**128-2**i-1) for i in range(128)]
    for u in inputs:
        res2 = compute_all_properties_from_input(GIFT128_builder, 11, u, extra_outmons=extra_outmons, nb_threads=nb_threads)
        if len(res2) > 0:
            res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
        else:
            res1 = [{}]
        res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    print(f"max number of properties on 12 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # searching lowest property 12 rounds
    # start = time()
    # res = []
    # f = GIFT128_builder(11)
    # f.get_ANF_propagation_model(F2n(128, 0), F2n(128, 0))
    # inputs = lowest_data_input_properties_var(f, 12, nb_threads=nb_threads)
    # print(inputs)
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(GIFT128_builder, 11, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"lowest data, highest number properties on 12 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")
