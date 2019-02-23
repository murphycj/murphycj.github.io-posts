"""
Given a set of parameters, simulates two rounds of PCR with a binomial
selection in between them.
"""

import random
import argparse
import time
import multiprocessing as mp

import numpy as np

MUTATOR = {
    '0': ['1', '2', '3'],
    '1': ['0', '2', '3'],
    '2': ['0', '1', '3'],
    '3': ['0', '1', '2']
}


def pcr_simulator(population, polymerase_error, pcr_efficiency, nbp, cycles):

    population_size = sum(population.values())
    start = time.time()

    range_bp = range(0, nbp)

    # pre draw large sample from binomial distribution for how many
    # mutations a sequence gets after PCR amplification
    # this step to save on some computational efficiency

    num_mutation_draws = 10000000
    num_mutations = np.random.binomial(
        nbp,
        polymerase_error,
        num_mutation_draws
    )
    num_mutations_counter = 0

    for cycle in range(0, cycles):

        population_t = population.items()
        pcr_sampling = np.random.binomial(
            [i[1] for i in population_t], pcr_efficiency)
        population_t = [
            [i[0], i[1], j] for i, j in zip(population_t, pcr_sampling)]

        for seq, count, pcr_sample in population_t:

            if pcr_sample > num_mutation_draws:
                num_mutation_draws = pcr_sample
                num_mutations = np.random.binomial(
                    nbp,
                    polymerase_error,
                    pcr_sample
                )
                num_mutations_counter = 0
            if (num_mutations_counter+pcr_sample) > num_mutation_draws:
                np.random.shuffle(num_mutations)
                num_mutations_counter = 0

            # sample how many mutations occur on the sequences chosen for
            # PCR amplification

            mutations_nn = num_mutations[
                num_mutations_counter:(num_mutations_counter+pcr_sample)]
            num_mutations_counter += pcr_sample

            # increment the current popluation by number of sequences
            # that do not get mutated
            population[seq] += np.count_nonzero(mutations_nn == 0)

            # randomly select base pair locations for mutation

            mutations_nn = mutations_nn[mutations_nn != 0]
            mutations_select = [
                random.sample(range_bp, i) for i in mutations_nn]

            for i in range(0, len(mutations_nn)):
                for mm in mutations_select[i]:
                    seq_new = seq[0:mm] + random.sample(
                        MUTATOR[seq[mm]], 1)[0] + seq[(mm+1):]
                if seq_new in population:
                    population[seq_new] += 1
                else:
                    population[seq_new] = 1
        # print 'Cycle: ' + str(cycle+1) + ' - ' + str(time.time() - start), str(sum(population.values()))

    return population, sum(population.values())


def simulation(x, args):
    population = {'0'*args.nbp: args.initpop}

    # REPLI-g simulation

    population, population_size = pcr_simulator(
        population=population,
        polymerase_error=args.polymerase_error_replig,
        pcr_efficiency=args.pcr_efficiency_replig,
        nbp=args.nbp,
        cycles=args.cycles_replig
    )

    # sample from the REPLI-g preduct for TruSeq PCR

    for seq, count in population.items():
        population[seq] = np.random.binomial(count, args.replig_sampling)
    for seq, count in population.items():
        if count == 0:
            population.pop(seq, None)

    population_size = sum(population.values())

    # TruSeq simulation

    population, population_size = pcr_simulator(
        population=population,
        polymerase_error=args.polymerase_error_truseq,
        pcr_efficiency=args.pcr_efficiency_truseq,
        nbp=args.nbp,
        cycles=args.cycles_truseq
    )

    # format the output data

    data = []
    for i in range(0, args.nbp):
        n_wt = 0
        n_mut_1 = 0
        n_mut_2 = 0
        n_mut_3 = 0
        for seq, count in population.items():
            if seq[i] == '0':
                n_wt += count
            elif seq[i] == '1':
                n_mut_1 += count
            elif seq[i] == '2':
                n_mut_2 += count
            elif seq[i] == '3':
                n_mut_3 += count
        total = n_wt + n_mut_1 + n_mut_2 + n_mut_3
        data += [n_mut_1, n_mut_2, n_mut_3]

    data = map(lambda x: str(x), data)
    print('Done with iteration: {}'.format(x))

    return str(total), data


def main(args):

    if args.p > 1:
        pool = mp.Pool(processes=args.p)
        results = [pool.apply_async(
            simulation, args=(x, args,)) for x in range(0, args.n)]
        results = [p.get() for p in results]
    else:
        results = [simulation(1, args)]

    totals = [i[0] for i in results]
    data = [i[1] for i in results]
    fout = open(args.prefix + '.csv', 'w')
    fout.write(','.join(map(lambda x: str(x), range(1, args.n+1))) + '\n')
    fout.write(','.join(totals) + '\n')
    for i in range(0, len(data[0])):
        fout.write(
            ','.join(map(lambda x: x[i], data)) + '\n'
        )
    fout.close()


parser = argparse.ArgumentParser()
parser.add_argument(
    '--initpop',
    type=int,
    required=True
)
parser.add_argument(
    '--polymerase_error_truseq',
    type=float,
    required=True
)
parser.add_argument(
    '--polymerase_error_replig',
    type=float,
    required=True
)
parser.add_argument(
    '--cycles_truseq',
    type=int,
    required=True
)
parser.add_argument(
    '--cycles_replig',
    type=int,
    required=True
)
parser.add_argument(
    '--pcr_efficiency_truseq',
    type=float,
    required=True
)
parser.add_argument(
    '--pcr_efficiency_replig',
    type=float,
    required=True
)
parser.add_argument(
    '--replig_sampling',
    type=float,
    required=True
)
parser.add_argument(
    '--nbp',
    type=int,
    required=False,
    default=3172
)
parser.add_argument(
    '--n',
    type=int,
    required=False,
    default=1,
    help='Number replicates.'
)
parser.add_argument(
    '--p',
    type=int,
    required=False,
    default=1,
    help='Number cores.'
)
parser.add_argument(
    '--prefix',
    type=str,
    required=True,
    help='Output prefix'
)
args = parser.parse_args()

main(args)
