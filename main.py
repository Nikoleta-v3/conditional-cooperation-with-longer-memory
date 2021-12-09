import numpy as np

import sympy as sym

import tqdm

import sys

import itertools


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


def payoffs_vector(c, b):
    return np.array([b - c, -c, b, 0] * 4)


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


def simulate_process(N, c, b, beta, max_steps):
    N = N
    c = c
    b = b
    beta = beta
    max_steps = max_steps

    residents = [[1, 1, 1, 1, 0, 1]]

    for t in tqdm.tqdm(range(1, max_steps)):
        current_resident = np.array(residents[-1][:4])

        mutant = np.random.random(4)

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
            residents.append(
                list(mutant) + [t] + [cooperation_rate]
            )

    return np.array(residents)


if __name__ == "__main__":
    N = 100
    c = float(sys.argv[1])
    b = 1
    beta = 1
    max_steps = 10 ** 5

    residents = simulate_process(N, c, b, beta, max_steps)

    np.savetxt(f"data/population_two_bits_c_{c}.csv", residents, delimiter=",")
