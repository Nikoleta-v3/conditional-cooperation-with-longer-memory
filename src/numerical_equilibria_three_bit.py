"""This script contains the code to check if a given strategy of n memory
is Nash in the space of all memory-n strategies."""

import itertools
import sys
from importlib.machinery import SourceFileLoader

import numpy as np
import pandas as pd
import sympy as sym
import tqdm

main = SourceFileLoader("main", "src/main.py").load_module()

eq = SourceFileLoader(
    "eq", "src/numerical_equilibria_n_bit_vs_n_bit.py"
).load_module()

from main import invariant_distribution


def calculate_M_memory_three(player, coplayer, analytical=False):

    row_iterations = [range(16), range(16, 32), range(32, 48), range(48, 64)]

    column_iterations = np.linspace(0, 63, 64).reshape(16, 4)

    if analytical == False:
        M = np.zeros((64, 64))

    elif analytical == True:
        M = sym.zeros(64, 64)

    player_probabilities = np.linspace(0, 63, 64).reshape(16, 4)
    coplayer_probabilities = np.linspace(0, 63, 64).reshape(16, 4)

    # adjusting co-player

    coplayer_probabilities[[1, 2]] = coplayer_probabilities[[2, 1]]

    coplayer_probabilities[[4, 8]] = coplayer_probabilities[[8, 4]]
    coplayer_probabilities[[5, 10]] = coplayer_probabilities[[10, 5]]
    coplayer_probabilities[[6, 9]] = coplayer_probabilities[[9, 6]]
    coplayer_probabilities[[7, 11]] = coplayer_probabilities[[11, 7]]

    coplayer_probabilities[[13, 14]] = coplayer_probabilities[[14, 13]]

    coplayer_probabilities[:, [1, 2]] = coplayer_probabilities[:, [2, 1]]

    for row in row_iterations:
        for i, irow in enumerate(row):
            pi = int(player_probabilities.flatten()[irow])
            qi = int(coplayer_probabilities.flatten()[irow])

            p = player[pi]
            q = coplayer[qi]

            combos = [p * q, p * (1 - q), (1 - p) * q, (1 - p) * (1 - q)]

            for j, column in enumerate(column_iterations[i, :]):
                M[irow, int(column)] = combos[j]
    return M


def match_payoff(player, coplayer, Sx):
    M = calculate_M_memory_three(player, coplayer)
    ss = invariant_distribution(M)

    return ss @ Sx


def task(i, strategy, coplayers, labels, filename, Sx, R, P, trans):

    sx = match_payoff(strategy, strategy, Sx)
    data = []

    for label, coplayer in zip(labels, coplayers):

        sy = match_payoff(coplayer, strategy, Sx)
        A = np.isclose(sx, sy, atol=10 ** -4) or sx > sy
        B = np.isclose(sy, R, atol=10 ** -4) or sy < R

        data_point = [
            i,
            *trans,
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
    max_simulation_number = 1
    R = .6
    P = .1
    folder = "prisoners_dilemma_n_three"

    deterministic_transitions = list(itertools.product([0, 1], repeat=8))

    deterministic_strategies = []
    for transition in deterministic_transitions:

        idx = [
            0,
            0,
            1,
            1,
            0,
            0,
            1,
            1,
            2,
            2,
            3,
            3,
            2,
            2,
            3,
            3,
            4,
            4,
            5,
            5,
            4,
            4,
            5,
            5,
            6,
            6,
            7,
            7,
            6,
            6,
            7,
            7,
        ]

        self_reactive = [transition[i] for i in idx] * 2
        deterministic_strategies.append(self_reactive)

    labels = [f"N{i}" for i, _ in enumerate(deterministic_strategies)]
    Sx = np.array([R, 0, 1, P] * 16)#eq.payoffs_donation(b, c, dim=16)

    for i in tqdm.tqdm(range(max_simulation_number)):
        np.random.seed(i)
        filename = f"{folder}/three_bit_iteration_{i}.csv"
        p1, p2, p3, p4, p5, p6, p7, p8 = np.random.random((1, 8)).round(5)[0]
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
            p5,
            p6,
            p5,
            p6,
            p7,
            p8,
            p7,
            p8,
            p5,
            p6,
            p5,
            p6,
            p7,
            p8,
            p7,
            p8,
        ] * 2
        task(
            i,
            strategy,
            deterministic_strategies,
            labels,
            filename,
            Sx,
            R,
            P,
            (p1, p2, p3, p4, p5, p6, p7, p8),
        )
