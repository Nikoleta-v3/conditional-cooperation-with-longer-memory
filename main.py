import numpy as np

import sympy as sym

import tqdm

import sys

import itertools


def transition_matrix_three_bits(p, q):

    p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8 = p
    q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8 = q

    M = np.zeros((64, 64))

    coop_combos = (
        list(itertools.product([q_1, q_2], [p_1, p_2]))
        + list(itertools.product([q_1, q_2], [p_3, p_4]))
        + list(itertools.product([q_3, q_4], [p_1, p_2]))
        + list(itertools.product([q_3, q_4], [p_3, p_4]))
    )
    entries = []
    for i, pair in enumerate(coop_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4), i] = entry

    cd_combos = (
        list(itertools.product([q_1, q_2], [p_5, p_6]))
        + list(itertools.product([q_1, q_2], [p_7, p_8]))
        + list(itertools.product([q_3, q_4], [p_5, p_6]))
        + list(itertools.product([q_3, q_4], [p_7, p_8]))
    )
    entries = []
    for i, pair in enumerate(cd_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 16, i] = entry

    dc_combos = (
        list(itertools.product([q_5, q_6], [p_1, p_2]))
        + list(itertools.product([q_5, q_6], [p_3, p_4]))
        + list(itertools.product([q_7, q_8], [p_1, p_2]))
        + list(itertools.product([q_7, q_8], [p_3, p_4]))
    )
    entries = []
    for i, pair in enumerate(dc_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 32, i] = entry

    def_combos = (
        list(itertools.product([q_5, q_6], [p_5, p_6]))
        + list(itertools.product([q_5, q_6], [p_7, p_8]))
        + list(itertools.product([q_7, q_8], [p_5, p_6]))
        + list(itertools.product([q_7, q_8], [p_7, p_8]))
    )
    entries = []
    for i, pair in enumerate(def_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 48, i] = entry

    return M


def transition_matrix_two_bit(p, q):
    p_1, p_2, p_3, p_4 = p
    q_1, q_2, q_3, q_4 = q

    M = np.zeros((16, 16))

    for i, expr in enumerate(
        [p_1 * q_1, p_1 * (1 - q_1), (1 - p_1) * q_1, (1 - p_1) * (1 - q_1)]
    ):
        M[0, i] = expr

    for i, expr in enumerate(
        [p_2 * q_1, p_2 * (1 - q_1), (1 - p_2) * q_1, (1 - p_2) * (1 - q_1)]
    ):
        M[1, 4 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_2, p_1 * (1 - q_2), (1 - p_1) * q_2, (1 - p_1) * (1 - q_2)]
    ):
        M[2, 8 + i] = expr

    for i, expr in enumerate(
        [p_2 * q_2, p_2 * (1 - q_2), (1 - p_2) * q_2, (1 - p_2) * (1 - q_2)]
    ):
        M[3, 12 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_1, p_3 * (1 - q_1), (1 - p_3) * q_1, (1 - p_3) * (1 - q_1)]
    ):
        M[4, i] = expr

    for i, expr in enumerate(
        [p_4 * q_1, p_4 * (1 - q_1), (1 - p_4) * q_1, (1 - p_4) * (1 - q_1)]
    ):
        M[5, 4 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_2, p_3 * (1 - q_2), (1 - p_3) * q_2, (1 - p_3) * (1 - q_2)]
    ):
        M[6, 8 + i] = expr

    for i, expr in enumerate(
        [p_4 * q_2, p_4 * (1 - q_2), (1 - p_4) * q_2, (1 - p_4) * (1 - q_2)]
    ):
        M[7, 12 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_3, p_1 * (1 - q_3), (1 - p_1) * q_3, (1 - p_1) * (1 - q_3)]
    ):
        M[8, i] = expr

    for i, expr in enumerate(
        [p_2 * q_3, p_2 * (1 - q_3), (1 - p_2) * q_3, (1 - p_2) * (1 - q_3)]
    ):
        M[9, 4 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_4, p_1 * (1 - q_4), (1 - p_1) * q_4, (1 - p_1) * (1 - q_4)]
    ):
        M[10, 8 + i] = expr

    for i, expr in enumerate(
        [p_2 * q_4, p_2 * (1 - q_4), (1 - p_2) * q_4, (1 - p_2) * (1 - q_4)]
    ):
        M[11, 12 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_3, p_3 * (1 - q_3), (1 - p_3) * q_3, (1 - p_3) * (1 - q_3)]
    ):
        M[12, i] = expr

    for i, expr in enumerate(
        [p_4 * q_3, p_4 * (1 - q_3), (1 - p_4) * q_3, (1 - p_4) * (1 - q_3)]
    ):
        M[13, 4 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_4, p_3 * (1 - q_4), (1 - p_3) * q_4, (1 - p_3) * (1 - q_4)]
    ):
        M[14, 8 + i] = expr

    for i, expr in enumerate(
        [p_4 * q_4, p_4 * (1 - q_4), (1 - p_4) * q_4, (1 - p_4) * (1 - q_4)]
    ):
        M[15, 12 + i] = expr

    return M


def transition_matrix_one_bit(p, q):
    """
    Returns a Markov transition matrix for a game of reactive strategies.
    """
    p_1, p_2 = p
    q_1, q_2 = q
    return np.array(
        [
            [
                p_1 * q_1,
                p_1 * (1 - q_1),
                q_1 * (1 - p_1),
                (1 - p_1) * (1 - q_1),
            ],
            [
                q_1 * p_2,
                p_2 * (1 - q_1),
                q_1 * (1 - p_2),
                (1 - q_1) * (1 - p_2),
            ],
            [
                p_1 * q_2,
                p_1 * (1 - q_2),
                q_2 * (1 - p_1),
                (1 - p_1) * (1 - q_2),
            ],
            [
                p_2 * q_2,
                p_2 * (1 - q_2),
                q_2 * (1 - p_2),
                (1 - p_2) * (1 - q_2),
            ],
        ],
    )


def invariant_distribution(M):

    eigenvalues, eigenvectors = np.linalg.eig(M.T)
    eigenvectors_one = eigenvectors[:, np.argmax(eigenvalues)]

    stationary = eigenvectors_one / eigenvectors_one.sum()

    return stationary.real


def invariant_distribution_analytically(M):
    size = M.shape[1]
    pi = sym.symbols(f"b_1:{size + 1}")
    ss = sym.solve(
        [sum(pi) - 1]
        + [a - b for a, b in zip(M.transpose() * sym.Matrix(pi), pi)],
        pi,
    )

    v_vector = sym.Matrix(
        [
            [ss[p] for p in pi],
        ]
    )

    return v_vector


def payoffs_vector(c, b, dim=4):
    return np.array([b - c, -c, b, 0] * dim)


def strategies_set(vector_size):
    """Set of determinist strategies for a given vector size."""
    actions = [0, 1]
    set_of_strategies = list(itertools.product(actions, repeat=vector_size))

    return np.array(set_of_strategies)


def transition_matrix_two_bit_analytical(p, q):
    p_1, p_2, p_3, p_4 = p
    q_1, q_2, q_3, q_4 = q

    M = sym.zeros(16, 16)

    for i, expr in enumerate(
        [p_1 * q_1, p_1 * (1 - q_1), (1 - p_1) * q_1, (1 - p_1) * (1 - q_1)]
    ):
        M[0, i] = expr

    for i, expr in enumerate(
        [p_2 * q_1, p_2 * (1 - q_1), (1 - p_2) * q_1, (1 - p_2) * (1 - q_1)]
    ):
        M[1, 4 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_2, p_1 * (1 - q_2), (1 - p_1) * q_2, (1 - p_1) * (1 - q_2)]
    ):
        M[2, 8 + i] = expr

    for i, expr in enumerate(
        [p_2 * q_2, p_2 * (1 - q_2), (1 - p_2) * q_2, (1 - p_2) * (1 - q_2)]
    ):
        M[3, 12 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_1, p_3 * (1 - q_1), (1 - p_3) * q_1, (1 - p_3) * (1 - q_1)]
    ):
        M[4, i] = expr

    for i, expr in enumerate(
        [p_4 * q_1, p_4 * (1 - q_1), (1 - p_4) * q_1, (1 - p_4) * (1 - q_1)]
    ):
        M[5, 4 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_2, p_3 * (1 - q_2), (1 - p_3) * q_2, (1 - p_3) * (1 - q_2)]
    ):
        M[6, 8 + i] = expr

    for i, expr in enumerate(
        [p_4 * q_2, p_4 * (1 - q_2), (1 - p_4) * q_2, (1 - p_4) * (1 - q_2)]
    ):
        M[7, 12 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_3, p_1 * (1 - q_3), (1 - p_1) * q_3, (1 - p_1) * (1 - q_3)]
    ):
        M[8, i] = expr

    for i, expr in enumerate(
        [p_2 * q_3, p_2 * (1 - q_3), (1 - p_2) * q_3, (1 - p_2) * (1 - q_3)]
    ):
        M[9, 4 + i] = expr

    for i, expr in enumerate(
        [p_1 * q_4, p_1 * (1 - q_4), (1 - p_1) * q_4, (1 - p_1) * (1 - q_4)]
    ):
        M[10, 8 + i] = expr

    for i, expr in enumerate(
        [p_2 * q_4, p_2 * (1 - q_4), (1 - p_2) * q_4, (1 - p_2) * (1 - q_4)]
    ):
        M[11, 12 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_3, p_3 * (1 - q_3), (1 - p_3) * q_3, (1 - p_3) * (1 - q_3)]
    ):
        M[12, i] = expr

    for i, expr in enumerate(
        [p_4 * q_3, p_4 * (1 - q_3), (1 - p_4) * q_3, (1 - p_4) * (1 - q_3)]
    ):
        M[13, 4 + i] = expr

    for i, expr in enumerate(
        [p_3 * q_4, p_3 * (1 - q_4), (1 - p_3) * q_4, (1 - p_3) * (1 - q_4)]
    ):
        M[14, 8 + i] = expr

    for i, expr in enumerate(
        [p_4 * q_4, p_4 * (1 - q_4), (1 - p_4) * q_4, (1 - p_4) * (1 - q_4)]
    ):
        M[15, 12 + i] = expr

    return M


def transition_matrix_one_bit_analytical(p, q):
    """
    Returns a Markov transition matrix for a game of reactive strategies.
    """
    p_1, p_2 = p
    q_1, q_2 = q
    return sym.Matrix(
        [
            [
                p_1 * q_1,
                p_1 * (1 - q_1),
                q_1 * (1 - p_1),
                (1 - p_1) * (1 - q_1),
            ],
            [
                q_1 * p_2,
                p_2 * (1 - q_1),
                q_1 * (1 - p_2),
                (1 - q_1) * (1 - p_2),
            ],
            [
                p_1 * q_2,
                p_1 * (1 - q_2),
                q_2 * (1 - p_1),
                (1 - p_1) * (1 - q_2),
            ],
            [
                p_2 * q_2,
                p_2 * (1 - q_2),
                q_2 * (1 - p_2),
                (1 - p_2) * (1 - q_2),
            ],
        ],
    )


def transition_matrix_three_bits_analytical(p, q):

    p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8 = p
    q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8 = q

    M = sym.zeros(64, 64)

    coop_combos = (
        list(itertools.product([q_1, q_2], [p_1, p_2]))
        + list(itertools.product([q_1, q_2], [p_3, p_4]))
        + list(itertools.product([q_3, q_4], [p_1, p_2]))
        + list(itertools.product([q_3, q_4], [p_3, p_4]))
    )
    entries = []
    for i, pair in enumerate(coop_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4), i] = entry

    cd_combos = (
        list(itertools.product([q_1, q_2], [p_5, p_6]))
        + list(itertools.product([q_1, q_2], [p_7, p_8]))
        + list(itertools.product([q_3, q_4], [p_5, p_6]))
        + list(itertools.product([q_3, q_4], [p_7, p_8]))
    )
    entries = []
    for i, pair in enumerate(cd_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 16, i] = entry

    dc_combos = (
        list(itertools.product([q_5, q_6], [p_1, p_2]))
        + list(itertools.product([q_5, q_6], [p_3, p_4]))
        + list(itertools.product([q_7, q_8], [p_1, p_2]))
        + list(itertools.product([q_7, q_8], [p_3, p_4]))
    )
    entries = []
    for i, pair in enumerate(dc_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 32, i] = entry

    def_combos = (
        list(itertools.product([q_5, q_6], [p_5, p_6]))
        + list(itertools.product([q_5, q_6], [p_7, p_8]))
        + list(itertools.product([q_7, q_8], [p_5, p_6]))
        + list(itertools.product([q_7, q_8], [p_7, p_8]))
    )
    entries = []
    for i, pair in enumerate(def_combos):
        actions = [
            (pair[1] * pair[0]),
            (pair[1] * (1 - pair[0])),
            ((1 - pair[1]) * pair[0]),
            ((1 - pair[1]) * (1 - pair[0])),
        ]
        for _, action in enumerate(actions):
            entries.append(action)
    for i, entry in enumerate(entries):
        M[int(i / 4) + 48, i] = entry

    return M


def simulate_process(N, c, b, beta, max_steps):
    N = N
    c = c
    b = b
    beta = beta
    max_steps = max_steps

    residents = [[0, 0, 0, 0, 0, 0]]

    for t in tqdm.tqdm(range(1, max_steps)):
        current_resident = np.array(residents[-1][:4])

        p, q = np.random.random(2)
        mutant = [p, q, p, q]

        steady_states = []

        for p, q in itertools.product([mutant, current_resident], repeat=2):

            M = transition_matrix_two_bit(p, q)
            steady_states.append(invariant_distribution(M))

        payoff_MM, payoff_MR, payoff_RM, payoff_RR = [
            state @ payoffs_vector(c, b) for state in steady_states
        ]

        lminus, lplus = [], []

        for k in range(1, N):
            expected_payoff_mutant = ((k - 1) / (N - 1) * payoff_MM) + (
                (N - k) / (N - 1)
            ) * payoff_MR
            expected_payoff_resident = (k / (N - 1) * payoff_RM) + (
                (N - k - 1) / (N - 1)
            ) * payoff_RR

            lplus.append(
                1
                / (
                    1
                    + np.exp(
                        float(
                            -beta
                            * (
                                expected_payoff_mutant
                                - expected_payoff_resident
                            )
                        )
                    )
                )
            )
            lminus.append(
                1
                / (
                    1
                    + np.exp(
                        float(
                            -beta
                            * (
                                expected_payoff_resident
                                - expected_payoff_mutant
                            )
                        )
                    )
                )
            )
        gammas = np.array(lminus) / np.array(lplus)

        if np.random.random() < 1 / (1 + np.sum(np.cumprod(gammas))):
            cooperation_rate = steady_states[0][0] + steady_states[0][1]
            residents.append(list(mutant) + [t] + [cooperation_rate])

    return np.array(residents)


if __name__ == "__main__":
    N = 100
    c = float(sys.argv[1])
    b = 1
    beta = 1
    max_steps = 10 ** 5

    residents = simulate_process(N, c, b, beta, max_steps)

    np.savetxt(f"data/population_one_bit_c_{c}.csv", residents, delimiter=",")
