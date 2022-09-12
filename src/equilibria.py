from email import header
import sys

import itertools
import glob
import numpy as np

import tqdm

import dask

import pandas as pd

import sympy as sym

from importlib.machinery import SourceFileLoader

main = SourceFileLoader("main", "src/main.py").load_module()

from main import invariant_distribution


def payoffs(R, P, dim=4):
    return np.array([R, 0, 1, P] * dim)


def coplayer_payoffs(R, P, dim=4):
    return np.array([R, 1, 0, P] * dim)


def payoffs_donation(b, c, dim=4):
    return np.array([b - c, -c, b, 0] * dim)


def coplayer_payoffs_donation(b, c, dim=4):
    return np.array([b - c, b, -c, 0] * dim)


def calculate_M_memory_two(player, coplayer):

    (
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        p7,
        p8,
        p9,
        p10,
        p11,
        p12,
        p13,
        p14,
        p15,
        p16,
    ) = player
    (
        q1,
        q2,
        q3,
        q4,
        q5,
        q6,
        q7,
        q8,
        q9,
        q10,
        q11,
        q12,
        q13,
        q14,
        q15,
        q16,
    ) = coplayer

    M = np.zeros((16, 16))

    col, row = 0, 0

    for p, q in [[p1, q1], [p2, q3], [p3, q2], [p4, q4]]:
        for i, combo in enumerate(
            [(p * q), ((1 - q) * p), ((1 - p) * q), ((1 - p) * (1 - q))]
        ):

            M[row, col + i] = combo

        col += 4
        row += 1

    col = 0
    for p, q in [[p5, q9], [p6, q11], [p7, q10], [p8, q12]]:
        for i, combo in enumerate(
            [(p * q), ((1 - q) * p), ((1 - p) * q), ((1 - p) * (1 - q))]
        ):

            M[row, col + i] = combo

        col += 4
        row += 1

    col = 0
    for p, q in [[p9, q5], [p10, q7], [p11, q6], [p12, q8]]:
        for i, combo in enumerate(
            [(p * q), ((1 - q) * p), ((1 - p) * q), ((1 - p) * (1 - q))]
        ):

            M[row, col + i] = combo

        col += 4
        row += 1

    col = 0
    for p, q in [[p13, q13], [p14, q15], [p15, q14], [p16, q16]]:
        for i, combo in enumerate(
            [(p * q), ((1 - q) * p), ((1 - p) * q), ((1 - p) * (1 - q))]
        ):

            M[row, col + i] = combo

        col += 4
        row += 1

    return M


def calculate_M(player, opponent):
    """
    Returns a Markov transition matrix for a game of memory one strategies.
    """
    return np.array(
        [
            [
                player[0] * opponent[0],
                player[0] * (1 - opponent[0]),
                (1 - player[0]) * opponent[0],
                (1 - player[0]) * (1 - opponent[0]),
            ],
            [
                player[1] * opponent[2],
                player[1] * (1 - opponent[2]),
                (1 - player[1]) * opponent[2],
                (1 - player[1]) * (1 - opponent[2]),
            ],
            [
                player[2] * opponent[1],
                player[2] * (1 - opponent[1]),
                (1 - player[2]) * opponent[1],
                (1 - player[2]) * (1 - opponent[1]),
            ],
            [
                player[3] * opponent[3],
                player[3] * (1 - opponent[3]),
                (1 - player[3]) * opponent[3],
                (1 - player[3]) * (1 - opponent[3]),
            ],
        ]
    )


def match_payoff(player, coplayer, Sx):
    M = calculate_M_memory_two(player, coplayer)
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
            *coplayer,
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
    max_simulation_number = 10 ** 5
    dimensions = int(sys.argv[1])
    b = 2
    c = 1
    # R = 0.6
    # P = 0.1
    seed = 0

    deterministic_strategies = list(
        itertools.product([0, 1], repeat=dimensions)
    )
    combos = list(itertools.product([0, 1], repeat=4))
    deterministic_strategies = []
    for (q1, q2, q3, q4) in combos:
        deterministic_strategies.append(
            [q1, q2, q1, q2, q3, q4, q3, q4, q1, q2, q1, q2, q3, q4, q3, q4]
        )

    labels = [f"N{i}" for i, _ in enumerate(deterministic_strategies)]
    Sx = payoffs(R, P, dim=4)  # payoffs(R, P, dim=4)

    np.random.seed(seed)
    jobs = []
    for i in tqdm.tqdm(range(max_simulation_number)):
        filename = f"two_bit_reactive_pd/dimensions_{dimensions}_iter_{i}_number_of_trials_{max_simulation_number}.csv"
        # strategy = np.random.random((1, dimensions)).round(5)[0]
        # strategy[0] = 1
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
                b,
                c,
            )
        )
    dask.compute(*jobs, nworkers=2)

    columns = (
        ["", "ID"]
        + [f"p{i+1}" for i in range(16)]
        + [f"q{i+1}" for i in range(16)]
        + ["label", "Sp", "Sq", "condition A", "condition B", "c", "b"]
    )

    files = glob.glob("../two_bit_reactive_pd/*csv")

    dfs = [pd.read_csv(file, index_col=0, names=columns) for file in files]

    df = pd.concat(dfs)

    df.to_csv("../two_bit_reactive_nash_pd.csv")
