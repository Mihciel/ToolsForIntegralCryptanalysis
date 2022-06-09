import sys
if not 'Tincrbell' in sys.modules:
    sys.path.append(".")
from Tincrbell.CompoundFunction import CompoundFunction, INPUT_ID, OUTPUT_ID
from Tincrbell.F2n import F2n
from Tincrbell.Components import SBox
from Tincrbell.CipherConstructor import CipherConstructor
from Tincrbell.SearchAlgorithms import get_balanced_output_bits, get_highest_round_balanced_bit, find_highest_round_with_property, get_all_properties_parity_set, compute_all_properties_from_input, has_a_property, lowest_data_input_properties, extended_search, lowest_data_input_properties_var, has_balanced_output_bit
from time import time
from functools import reduce
import random
from itertools import combinations

PRESENT_sbox = SBox(4, 4, [0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2])
PRESENT_bitpermutation = [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63]
PRESENT_roundfunction = CompoundFunction(64, 64)
ids = []
# add sboxes
for _ in range(16):
    ids.append(PRESENT_roundfunction.add_component(PRESENT_sbox))
# connect input to components
for i in range(64):
    PRESENT_roundfunction.connect_components(INPUT_ID, i, ids[i//4], i%4)
# connect components to output
for i in range(64):
    PRESENT_roundfunction.connect_components(ids[i//4], i%4, OUTPUT_ID, PRESENT_bitpermutation[i])
PRESENT_builder = CipherConstructor(PRESENT_roundfunction, F2n(64, 2**64-1))

if __name__=="__main__":

    if len(sys.argv) > 1:
        nb_threads = int(sys.argv[1])
    else:
        nb_threads = int(input("Please input the number of threads to use: "))

    # search max round

    # start = time()
    # res = get_highest_round_balanced_bit(PRESENT_builder, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]} in {time()-start}s")
    # # rounds: 9
    # # input: F2n(64, 0xfffffffffffffffe)
    # # output bit: 0

    # start = time()
    # res = find_highest_round_with_property(PRESENT_builder, start_with=res, nb_threads=nb_threads)
    # print(f"rounds: {res[0]}\ninput: {res[1]}\noutput bit: {res[2]} in {time()-start}s")
    # # rounds: 9
    # # input set: frozenset({F2n(64, 0xfffffffffffffffe)})
    # # output set: frozenset({F2n(64, 0x0000000000000001)})

    # Check literature results
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(3), F2n(64, 0x1111), nb_threads=nb_threads)
    # print(f"Verifying Z'aba et al. (2008): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(6), F2n(64, 0x0fffffffffffffff), nb_threads=nb_threads)
    # print(f"Verifying Todo (2015): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(9), F2n(64, 0xfffffffffffffff0), nb_threads=nb_threads)
    # print(f"Verifying Xiang et al. (2016): {len(res) == 1 and 0 in res} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(4), F2n(64, 0xf), nb_threads=nb_threads)
    # print(f"Verifying 4-round Boura and Canteaut (2016): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(5), F2n(64, 0xfff0), nb_threads=nb_threads)
    # print(f"Verifying 5-round Boura and Canteaut (2016): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(6), F2n(64, 0xffffffff), nb_threads=nb_threads)
    # print(f"Verifying 6-round Boura and Canteaut (2016): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(7), F2n(64, 0xfffffffffffff000), nb_threads=nb_threads)
    # print(f"Verifying 7-round Boura and Canteaut (2016): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(8), F2n(64, 0xfffffffffffffffe), nb_threads=nb_threads)
    # print(f"Verifying 8-round Boura and Canteaut (2016): {len(res) == 64} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(6), F2n(64, 0xfff0), nb_threads=nb_threads)
    # print(f"Verifying 6-round reduced advantage Boura and Canteaut (2016): {len(res) == 16} in {time()-start}s")
    # start = time()
    # res = get_balanced_output_bits(PRESENT_builder(9), F2n(64, 0xfffffffffffffffe), nb_threads=nb_threads)
    # print(f"Verifying Hu et al. (2018): {len(res) == 28} in {time()-start}s")


    # improve literature results, that do not require linear combinations
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(3), F2n(64, 0x1111), timeout=120, nb_threads=nb_threads)
    # print(f"Improving Z'aba et al. (2008): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(6), F2n(64, 0x0fffffffffffffff), timeout=120, nb_threads=nb_threads)
    # print(f"Improving Todo (2015): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(4), F2n(64, 0xf), timeout=120, nb_threads=nb_threads)
    # print(f"Improving 4-round Boura and Canteaut (2016): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(5), F2n(64, 0xfff0), timeout=120, nb_threads=nb_threads)
    # print(f"Improving 5-round Boura and Canteaut (2016): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(6), F2n(64, 0xffffffff), timeout=120, nb_threads=nb_threads)
    # print(f"Improving 6-round Boura and Canteaut (2016): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(7), F2n(64, 0xfffffffffffff000), timeout=120, nb_threads=nb_threads)
    # print(f"Improving 7-round Boura and Canteaut (2016): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")
    # start = time()
    # res = get_all_properties_parity_set(PRESENT_builder(8), F2n(64, 0xfffffffffffffffe), timeout=120, nb_threads=nb_threads)
    # print(f"Improving 8-round Boura and Canteaut (2016): {len(reduce(set.union, res, set())), list(map(len, res))} in {time()-start}s")

    # improve literature results, requires linear combination
    # start = time()
    # res1 = get_all_properties_parity_set(PRESENT_builder(6), F2n(64, 0xfff0), timeout=120, nb_threads=nb_threads)
    # res2 = compute_all_properties_from_input(PRESENT_builder, 6,  F2n(64, 0xfff0), nb_threads=nb_threads)
    # print(f"Improving 6-round reduced advantage Boura and Canteaut (2016): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 9,  F2n(64, 0xfffffffffffffff0), nb_threads=nb_threads)
    # print(f"Improving Xiang et al. (2016): {1, [1], len(res2)} in {time()-start}s")
    # start = time()
    # res1 = get_all_properties_parity_set(PRESENT_builder(9), F2n(64, 0xfffffffffffffffe), timeout=120, nb_threads=nb_threads)
    # res2 = compute_all_properties_from_input(PRESENT_builder, 9,  F2n(64, 0xfffffffffffffffe), nb_threads=nb_threads)
    # print(f"Improving Hu et al. (2018): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")

    # improve literature results by extending with one round
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 4, F2n(64, 0x1111), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(4), F2n(64, 0x1111), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending Z'aba et al. (2008): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 7, F2n(64, 0x0fffffffffffffff), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(7), F2n(64, 0x0fffffffffffffff), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending Todo (2015): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 5, F2n(64, 0xf), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(5), F2n(64, 0xf), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending 4-round Boura and Canteaut (2016): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 7, F2n(64, 0xffffffff), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(7), F2n(64, 0xffffffff), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending 6-round Boura and Canteaut (2016): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 8, F2n(64, 0xfffffffffffff000), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(8), F2n(64, 0xfffffffffffff000), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending 7-round Boura and Canteaut (2016): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 7, F2n(64, 0xfff0), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(7), F2n(64, 0xfff0), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"Extending 6-round reduced advantage Boura and Canteaut (2016): {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")

    # highest number of properties on 9 rounds
    # start = time()
    # res = []
    # f = PRESENT_builder(9)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # inputs = [F2n(64, 2**64-2**i-1) for i in range(64)]
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(PRESENT_builder, 9, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, timeout=120, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"max number of properties on 9 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # lowest data on 9 rounds
    # start = time()
    # res = []
    # f = PRESENT_builder(9)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # inputs = lowest_data_input_properties(f, nb_threads=nb_threads)
    # print(inputs)
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(PRESENT_builder, 9, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"lowest data, highest number properties on 9 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # lowest data on 8 rounds
    # start = time()
    # res = []
    # f = PRESENT_builder(8)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # inputs = lowest_data_input_properties_var(f, 12, nb_threads=nb_threads)
    # print(inputs)
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(PRESENT_builder, 8, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"lowest data, highest number properties on 8 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # heuristic lower data search on 7 rounds
    # start = time()
    # inputs_subset = [F2n(64, 0xf<<(4*i)) for i in range(16)]
    # nb_saturated_sboxes = 12
    # flag=True
    # f = PRESENT_builder(7)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # while flag:
    #     flag=False
    #     nb_saturated_sboxes-=1
    #     for l in combinations(inputs_subset, nb_saturated_sboxes):
    #         u = reduce(F2n.__xor__, l)
    #         if has_balanced_output_bit(f, u)[0]:
    #             flag = True
    #             break
    # nb_saturated_sboxes+=1
    # inputs = [u for u in map(lambda l: reduce(F2n.__xor__, l), combinations(inputs_subset, nb_saturated_sboxes)) if has_balanced_output_bit(f, u)]
    # res = []
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(PRESENT_builder, 7, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"low data, highest number properties on 7 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")

    # heuristic lower data search on 6 rounds
    # start = time()
    # inputs_subset = [F2n(64, 0xf<<(4*i)) for i in range(16)]
    # nb_saturated_sboxes = 4
    # flag=True
    # f = PRESENT_builder(6)
    # f.get_ANF_propagation_model(F2n(64, 0), F2n(64, 0))
    # while flag:
    #     flag=False
    #     nb_saturated_sboxes-=1
    #     for l in combinations(inputs_subset, nb_saturated_sboxes):
    #         u = reduce(F2n.__xor__, l)
    #         if has_balanced_output_bit(f, u)[0]:
    #             flag = True
    #             break
    # nb_saturated_sboxes+=1
    # inputs = [u for u in map(lambda l: reduce(F2n.__xor__, l), combinations(inputs_subset, nb_saturated_sboxes)) if has_balanced_output_bit(f, u)]
    # res = []
    # for u in inputs:
    #     res2 = compute_all_properties_from_input(PRESENT_builder, 6, u, nb_threads=nb_threads)
    #     if len(res2) > 0:
    #         res1 = get_all_properties_parity_set(f, u, nb_threads=nb_threads)
    #     else:
    #         res1 = [{}]
    #     res.append((len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)))
    # print(f"low data, highest number properties on 6 rounds is {max(zip(inputs, res), key=lambda x: x[1][0] + x[1][2] - x[1][1][0])} in {time()-start}s")
    
    # search for 4 round distinguisher with 2^4 data and highest nb properties
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 4, F2n(64, 0x8000000000000011), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(4), F2n(64, 0x8000000000000011), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"low data, highest number properties on 4 rounds is {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")

    # search for 3 round distinguisher with 2^4 data and highest nb properties
    # start = time()
    # res2 = compute_all_properties_from_input(PRESENT_builder, 3, F2n(64, 0x3), nb_threads=nb_threads)
    # if len(res2) > 0:
    #     res1 = get_all_properties_parity_set(PRESENT_builder(3), F2n(64, 0x3), timeout=120, nb_threads=nb_threads)
    # else:
    #     res1 = [{}]
    # print(f"low data, highest number properties on 3 rounds is {len(reduce(set.union, res1, set())), list(map(len, res1)), len(res2)} in {time()-start}s")
    
    # extended search
    # start = time()
    # res = extended_search(PRESENT_builder, 10, nb_threads=nb_threads)
    # print(res)
    # print(f"extended search: found {len(res)} linearly independent results in {time()-start}s")

    # test the results from 6 round 2**12 inputs
    # f = PRESENT_builder(6)
    # random.seed(31415)
    # u = F2n(64, 0xfff0)
    # res1 = reduce(set.union, get_all_properties_parity_set(PRESENT_builder(6), u, timeout=120, nb_threads=nb_threads))
    # res2 = compute_all_properties_from_input(PRESENT_builder, 6, u, nb_threads=nb_threads)
    # for _ in range(2**10):
    #     f.seed = random.randrange(2**128)
    #     xs = [f(x) for x in u.supp()]
    #     for v in res1:
    #         assert(sum(map(lambda x: x**v, xs))%2 == 0)
    #     for vs in res2:
    #         assert(sum(map(lambda x: sum(map(lambda v: x**v, vs)), xs))%2==0)

    # test the results from 7 round 2**18 inputs
    f = PRESENT_builder(7)
    random.seed(31415)
    u = F2n(64, 0xffff)
    res1 = reduce(set.union, get_all_properties_parity_set(PRESENT_builder(7), u, timeout=120, nb_threads=nb_threads))
    res2 = compute_all_properties_from_input(PRESENT_builder, 7, u, nb_threads=nb_threads)
    for _ in range(2**6):
        f.seed = random.randrange(2**128)
        xs = [f(x) for x in u.supp()]
        for v in res1:
            assert(sum(map(lambda x: x**v, xs))%2 == 0)
        for vs in res2:
            assert(sum(map(lambda x: sum(map(lambda v: x**v, vs)), xs))%2==0)