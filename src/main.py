import numpy as np

import sympy as sym

import tqdm

import sys

import itertools


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


def transition_matrix(p, q, bits, analytical=False):

    shape = bits ** 2

    ones = [
        p[j - int(bits / 2) : j]
        for _, j in itertools.product([int(bits / 2), bits], repeat=2)
    ]

    twos = [
        q[i - int(bits / 2) : i]
        for i, _ in itertools.product([int(bits / 2), bits], repeat=2)
    ]

    twos = [[two[i : i + 2] for i in np.arange(0, len(two), 2)] for two in twos]

    ones = [[one[i : i + 2] for i in np.arange(0, len(one), 2)] for one in ones]

    twos = [[two[i : i + 2] for i in np.arange(0, len(two), 2)] for two in twos]

    ones = [[one[i : i + 2] for i in np.arange(0, len(one), 2)] for one in ones]

    if analytical == True:
        M = sym.zeros(shape, shape)
    else:
        M = np.zeros((shape, shape))

    row = 0
    for two, one in zip(twos, ones):
        for pi, qj in itertools.product(two, one, repeat=1):
            for pj, qi in itertools.product(pi, qj, repeat=1):
                for i, combo in enumerate(itertools.product(pj, qi, repeat=1)):
                    for j, expr in enumerate(
                        [
                            combo[0] * combo[1],
                            (1 - combo[0]) * combo[1],
                            (1 - combo[1]) * combo[0],
                            (1 - combo[1]) * (1 - combo[0]),
                        ]
                    ):

                        M[row, ((row * 4) % shape) + j] = expr

                    row += 1
    return M


def cooperation_rate(ss, size):
    return sum(
        [s for i, s in enumerate(ss) if i in np.arange(0, size, 4)]
        + [s for i, s in enumerate(ss) if i in np.arange(1, size, 4)]
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


def payoffs_vector_coplayer(c, b, dim=4):
    return np.array([b - c, b, -c, 0] * dim)


def payoffs_vector_coplayer_prime(c, b, dim=4):
    return np.array([b - c] * dim + [b] * dim +  [-c] * dim + [0] * dim)


def strategies_set(vector_size):
    """Set of determinist strategies for a given vector size."""
    actions = [0, 1]
    set_of_strategies = list(itertools.product(actions, repeat=vector_size))

    return np.array(set_of_strategies)


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
