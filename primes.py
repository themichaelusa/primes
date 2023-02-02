import math
import pathlib
from collections import defaultdict

MAX_INT = (2**63) - 1

def is_prime(n, primes):
    for prime in primes:
        if n % prime == 0:
            return False
    return True

def primes_to_nth_prime(n):
    primes = [2]
    current_index = 0
    while(len(primes) < n):
        if current_index == 0:
            start_value = 3
            end_value = 2 * 2
        else:
            start_value = primes[current_index - 1] * primes[current_index - 1] + 1
            end_value = primes[current_index] * primes[current_index]
        for i in range(start_value, end_value):
            if is_prime(i, primes):
                primes.append(i)
        current_index += 1
    return primes[:n]

def prime_frequency_map(fpath):
    primes = primes_to_nth_prime(256)
    primes_map = {k:v for k,v in enumerate(primes)}
    freq_map = defaultdict(int)
   
    # construct prime map
    byte_array = []
    with open(fpath, "rb") as f:
        data = f.read(8192)
        while data:
            data = f.read(8192)
            for byte in data:
                freq_map[primes_map[byte]] += 1
                byte_array.append(byte)

    # get occurence map of valid sequences e.g prime multiples
    occurence_map_largest = defaultdict(int)
    occurence_map_mults_cnts = defaultdict(int)
    occurence_map_seqs = defaultdict(list)
    bytes_idx = 1
    accum = byte_array[0]
    largest_term = []
    while bytes_idx < len(byte_array):
        prime_byte = primes_map[byte_array[bytes_idx]]
        prime_prod = accum * prime_byte
        if prime_prod > MAX_INT:
            occurence_map_mults_cnts[accum] += 1
            occurence_map_largest[accum] = max(largest_term)
            if len(occurence_map_seqs[accum]) == 0:
                occurence_map_seqs[accum].append([0])
            last_seqs_term = occurence_map_seqs[accum][-1]
            if last_seqs_term != largest_term:
                occurence_map_seqs[accum].append(largest_term)
            accum = 1
            largest_term = []
        accum *= prime_byte
        largest_term.append(prime_byte)
        bytes_idx += 1


    occ_seqs_cnts = list(occurence_map_seqs.items())
    occ_seqs_cnts.sort(key = lambda l: len(l[1]))
    total_unique_cnts = [(occ[0], len(occ[1])-1) for occ in occ_seqs_cnts]

    #occ_seqs_cnts.reverse()
    #total_unique_cnts = sum([len(occ[1])-1 for occ in occ_seqs_cnts])
    #print(total_unique_cnts)

    """
    occ_cnts = list(occurence_map_mults_cnts.items())
    occ_cnts.sort(key = lambda l: l[1])
    occ_cnts.reverse()
    total_unique_cnts = sum([occ[1] for occ in occ_cnts])
    """

    #print("occurence_map_largest", occurence_map_largest)
    #print("occurence_map_mults_cnts", occurence_map_mults_cnts)

    # get primes as keys
    #fq = list(freq_map.items())
    #fq.sort(key = lambda l: l[1])
    #return {"freq" : fq, "product_freq" : None}


if __name__ == '__main__':
    print(prime_frequency_map('./test.csv'))

