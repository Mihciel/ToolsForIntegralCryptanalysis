import sys
if not 'Tincrbell' in sys.modules:
    sys.path.append(".")
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID, ERROR_ID
from Tincrbell.F2n import F2n
from Tincrbell.Components import SBox
from Tincrbell.CipherConstructor import CipherConstructor
from Tincrbell.SearchAlgorithms import get_highest_round_balanced_bit, find_highest_round_with_property, get_balanced_output_bits, get_all_properties_parity_set, compute_all_properties_from_input, lowest_data_input_properties_var, has_balanced_output_bit, extended_search
from time import time
from itertools import product, combinations
from functools import reduce

GIFT_sbox = SBox(4, 4, [1, 0xa, 4, 0xc, 6, 0xf, 3, 9, 2, 0xd, 0xb, 7, 5, 0, 8, 0xe])
# there is only one way the round constant is merged into the sbox:
GIFT_modsbox = SBox(4, 4, [x ^ 8 for x in [1, 0xa, 4, 0xc, 6, 0xf, 3, 9, 2, 0xd, 0xb, 7, 5, 0, 8, 0xe]])
GIFT64_bitpermutation = [4*(i//16) + 16*((3 * ((i%16)//4) + (i % 4))%4) + (i % 4) for i in range(64)]
GIFT_round_constants = [F2n(6, x) for x in [0x01,0x03,0x07,0x0F,0x1F,0x3E,0x3D,0x3B,0x37,0x2F,0x1E,0x3C,0x39,0x33,0x27,0x0E,0x1D,0x3A,0x35,0x2B,0x16,0x2C,0x18,0x30,0x21,0x02,0x05,0x0B,0x17,0x2E,0x1C,0x38,0x31,0x23,0x06,0x0D,0x1B,0x36,0x2D,0x1A,0x34,0x29,0x12,0x24,0x08,0x11,0x22,0x04]]
GIFT64_round_functions = []
for c in GIFT_round_constants:
    GIFT64_round_functions.append(CompoundFunction(64, 64))
    for i in range(16):
        sboxid = ERROR_ID
        if (GIFT64_bitpermutation[4*i+3]//4 < 6 and c[GIFT64_bitpermutation[4*i+3]//4] == 1) or GIFT64_bitpermutation[4*i+3]//4 == 15:
            sboxid = GIFT64_round_functions[-1].add_component(GIFT_modsbox)
        else:
            sboxid = GIFT64_round_functions[-1].add_component(GIFT_sbox)
        for j in range(4):
            GIFT64_round_functions[-1].connect_components(INPUT_ID, 4*i+j, sboxid, j)
            GIFT64_round_functions[-1].connect_components(sboxid, j, OUTPUT_ID, GIFT64_bitpermutation[4*i+j])
# don't add the S and P layer of the first round, only the key addition
# since these can trivially be inverted
GIFT64_builder = CipherConstructor(tuple(GIFT64_round_functions[1:]), F2n(64, 0x3333333333333333))


if __name__=="__main__":
    if len(sys.argv) > 1:
        nb_threads = int(sys.argv[1])
    else:
        nb_threads = int(input("Please input the number of threads to use: "))
    
    # check highest possible round

    # start = time()
    # res = get_highest_round_balanced_bit(GIFT64_builder, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]} in {time()-start}s")
    # rounds: 9
    # input: F2n(64, 0xfffffffffffffffe)
    # output bit: 4

    # start = time()
    # extra_inmons = tuple(F2n(4*i, 2**(4*i)-1).concatenate(F2n(4, 0x3)).concatenate(F2n(4*(15-i), 2**(4*(15-i))-1)) for i in range(16))
    # extra_outmons = tuple(~x for x in extra_inmons)
    # res = find_highest_round_with_property(GIFT64_builder, extra_inmons=extra_inmons, extra_outmons=extra_outmons, start_with=res, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput set: {res[1]}\noutput set: {res[2]} in {time()-start}s")
    # rounds: 9
    # input set: frozenset({F2n(64, 0xfffffffffffffff)})
    # output set: frozenset({F2n(64, 0x0000000000000010)})
    
    # validate litaratue
    # start = time()
    # res = get_balanced_output_bits(GIFT64_builder(9), F2n(64, 0xfffffffffffffffb), nb_threads=nb_threads)
    # print(f"Banik et al. 10 round 2^63 data (2017): {len(res) >= 32} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(GIFT64_builder(9), F2n(64, 0xfffffffffffffff8), nb_threads=nb_threads)
    # print(f"Eskandari et al. 10 round 2^61 data (2017): {len(res) >= 5} in {time()-start}s")

    # improving literature results
    # start = time()
    # res1 = get_all_properties_parity_set(GIFT64_builder(9), F2n(64, 0xfffffffffffffffb), timeout=120, nb_threads=nb_threads)
    # res2 = compute_all_properties_from_input(GIFT64_builder, 9,  F2n(64, 0xfffffffffffffffb), nb_threads=nb_threads)
    # print(f"Improving 10-round 2^63 Banik et al. (2017): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res1 = get_all_properties_parity_set(GIFT64_builder(9), F2n(64, 0xfffffffffffffff8), timeout=120, nb_threads=nb_threads)
    # res2 = compute_all_properties_from_input(GIFT64_builder, 9,  F2n(64, 0xfffffffffffffff8), nb_threads=nb_threads)
    # print(f"Improving 10-round 2^61 Eskandari et al. (2018): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")

    # searching for highest number of properties on 10 rounds
    # start = time()
    # extra_inmons = tuple(F2n(4*i, 2**(4*i)-1).concatenate(F2n(4, 0x3)).concatenate(F2n(4*(15-i), 2**(4*(15-i))-1)) for i in range(16))
    # extra_outmons = tuple(~x for x in extra_inmons)
    # res = []
    # f = GIFT64_builder(9)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # inputs = [F2n(64, 2**64-2**i-1) for i in range(64)]
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(GIFT64_builder, 9, u, extra_outmons=extra_outmons, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, timeout=120, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"max number of properties on 9 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # searching lowest property 10 rounds
    # start = time()
    # res = []
    # f = GIFT64_builder(9)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # inputs = lowest_data_input_properties_var(f, 9, nb_threads=nb_threads)
    # print(inputs)
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(GIFT64_builder, 9, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"lowest data, highest number properties on 9 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")
    # lowest data, highest number properties on 9 rounds is (F2n(64, 0xfffffff1ffffffff), (11, [11, 0], 11))

    # extended search
    # start = time()
    # extra_inmons = tuple(F2n(4*i, 2**(4*i)-1).concatenate(F2n(4, 0x3)).concatenate(F2n(4*(15-i), 2**(4*(15-i))-1)) for i in range(16))
    # extra_outmons = tuple(~x for x in extra_inmons)
    # res = extended_search(GIFT64_builder, 10, extra_inmons=extra_inmons, extra_outmons=extra_outmons, nb_threads=nb_threads)
    # print(res)
    # print(f"extended search: found {len(res)} linearly independent results in {time()-start}s")