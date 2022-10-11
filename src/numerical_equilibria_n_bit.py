"""This script contains the code to check if a given strategy of n memory
is Nash in the space of all memory-n strategies."""

import sys

import itertools
import numpy as np

import tqdm

import dask

import pandas as pd

from importlib.machinery import SourceFileLoader

main = SourceFileLoader("main", "src/main.py").load_module()

eq = SourceFileLoader("eq", "src/numerical_equilibria_n_bit_vs_n_bit.py").load_module()

from main import invariant_distribution

def match_payoff(player, coplayer, Sx):
    M = eq.calculate_M_memory_two(player, coplayer)
    ss = invariant_distribution(M)

    return ss @ Sx


def task(i, strategy, coplayers, labels, filename, Sx, R, P):

    sx = match_payoff(strategy, strategy, Sx)
    data = []

    for label, coplayer in zip(labels, coplayers):

        sy = match_payoff(coplayer, strategy, Sx)
        A = np.isclose(sx, sy, atol=10 ** -4) or sx > sy
        B = np.isclose(sy, R, atol=10 ** -4) or sy < R

        data_point = [
            i,
            *strategy,
            label,
            sx,
            sy,
            A,
            B,
            R,
            P,
        ]
        data.append(data_point)

    df = pd.DataFrame(data)
    df.to_csv(filename, header=False)


if __name__ == "__main__":
    max_simulation_number = 1000
    dimensions = int(sys.argv[1])
    # b = 2
    # c = 1
    n = 5
    R = 0.6
    P = 0.1
    seed = 0
    folder = "prisoners_dilemma_two_bit_reactive" # "two_bit_against_memory_two"
    lbound = 100
    ubound = 200

    deterministic_strategies = list(
        itertools.product([0, 1], repeat=2 ** (2 * dimensions))
    )

    labels = [f"N{i}" for i, _ in enumerate(deterministic_strategies)]
    Sx = eq.payoffs(R, P, dim=4)
    #np.random.seed(seed)
    # steps = np.arange(100, max_simulation_number, 100)
    
    # for lbound, ubound in zip(steps[:-1], steps[1:]):
    jobs = []
    for i in tqdm.tqdm(range(lbound, ubound)):
        np.random.seed(i)
        filename = f"{folder}/dimensions_{dimensions}_iter_{i}_number_of_trials_{max_simulation_number}.csv"
        p1, p2, p3, p4 = np.random.random((1, 4)).round(5)[0]
        p1 = 1
        strategy = [
                p1,
                p2,
                p1,
                p2,
                p3,
                p4,
                p3,
                p4,
                p1,
                p2,
                p1,
                p2,
                p3,
                p4,
                p3,
                p4,
            ]
        jobs.append(
                    dask.delayed(task)(
                        i,
                        strategy,
                        deterministic_strategies,
                        labels,
                        filename,
                        Sx,
                        R,
                        P,
                    )
                )
    dask.compute(*jobs, nworkers=n)
