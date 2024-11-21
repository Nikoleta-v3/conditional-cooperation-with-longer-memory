import sympy as sym

import numpy as np

import repeated_play


import itertools

from tqdm import tqdm

import pandas as pd


def trnsf_transition_m_memory_two(player, analytical=True):
    if analytical == True:
        M = sym.Matrix(
            [
                [player[0], (1 - player[0]), 0, 0],
                [0, 0, player[1], (1 - player[1])],
                [player[2], (1 - player[2]), 0, 0],
                [0, 0, player[3], (1 - player[3])],
            ]
        )
    else:
        M = np.array(
            [
                [player[0], (1 - player[0]), 0, 0],
                [0, 0, player[1], (1 - player[1])],
                [player[2], (1 - player[2]), 0, 0],
                [0, 0, player[3], (1 - player[3])],
            ]
        )

    return M


def trnsf_transition_m_memory_three(player, analytical=True):
    if analytical == False:
        return np.array(
            [
                [player[0], (1 - player[0]), 0, 0, 0, 0, 0, 0],
                [0, 0, player[1], (1 - player[1]), 0, 0, 0, 0],
                [0, 0, 0, 0, player[2], (1 - player[2]), 0, 0],
                [0, 0, 0, 0, 0, 0, player[3], (1 - player[3])],
                [player[4], (1 - player[4]), 0, 0, 0, 0, 0, 0],
                [0, 0, player[5], (1 - player[5]), 0, 0, 0, 0],
                [0, 0, 0, 0, player[6], (1 - player[6]), 0, 0],
                [0, 0, 0, 0, 0, 0, player[7], (1 - player[7])],
            ]
        )

    if analytical == True:
        return sym.Matrix(
            [
                [player[0], (1 - player[0]), 0, 0, 0, 0, 0, 0],
                [0, 0, player[1], (1 - player[1]), 0, 0, 0, 0],
                [0, 0, 0, 0, player[2], (1 - player[2]), 0, 0],
                [0, 0, 0, 0, 0, 0, player[3], (1 - player[3])],
                [player[4], (1 - player[4]), 0, 0, 0, 0, 0, 0],
                [0, 0, player[5], (1 - player[5]), 0, 0, 0, 0],
                [0, 0, 0, 0, player[6], (1 - player[6]), 0, 0],
                [0, 0, 0, 0, 0, 0, player[7], (1 - player[7])],
            ]
        )


def match_payoffs_efficiently(player, coplayer, b, c, memory):
    if memory == "two":
        M = trnsf_transition_m_memory_two(coplayer, analytical=False)
    elif memory == "three":
        M = trnsf_transition_m_memory_three(coplayer, analytical=False)

    ss = repeated_play.stationary_distribution(M)[0]

    if memory == "two":
        rho_q = ss[0] + ss[1]

    elif memory == "three":
        rho_q = ss[0] + ss[1] + ss[4] + ss[5]

    rho_p = sum([ss[i] * p for i, p in enumerate(player)])

    payoff_q = rho_p * b - c * rho_q

    return payoff_q


def match_payoffs_markov_approach(player, coplayer, b, c, memory):
    assert memory == "one"

    M = repeated_play.transition_matrix_repeated_game(
        player, coplayer, memory=memory
    )

    ss = repeated_play.stationary_distribution(M)[0]

    return ss @ np.array([b - c, -c, b, 0])


def match_payoffs(player, coplayer, b, c, memory):
    if memory == "one":
        assert len(player) == 4
        assert len(coplayer) == 4
        return match_payoffs_markov_approach(player, coplayer, b, c, memory)

    if memory == "two":
        assert len(player) == 4
        assert len(coplayer) == 4
        return match_payoffs_efficiently(player, coplayer, b, c, memory)

    if memory == "three":
        assert len(player) == 8
        assert len(coplayer) == 8
        return match_payoffs_efficiently(player, coplayer, b, c, memory)
