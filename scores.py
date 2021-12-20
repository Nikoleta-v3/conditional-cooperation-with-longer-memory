import axelrod as axl
import numpy as np

import sys

if __name__ == "__main__":

    for repetitions in [1, 100, 1000]:
        Delta = 1

        players = [axl.CURE(Delta), axl.CURE(Delta)]

        match = axl.Tournament(
            players=players, turns=10 ** 6, repetitions=repetitions, noise=0.1
        )
        results = match.play(progress_bar=True)

        score = np.mean(results.normalised_scores[0])

        with open(f"reps_{repetitions}_score.txt", "w") as f:
            f.write(f"{score}")
