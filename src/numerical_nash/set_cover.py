import pandas as pd

import tqdm

import itertools

import random

import copy

class DictList(dict):
    def __setitem__(self, key, value):
        try:
            # Assumes there is a list on the key
            self[key].append(value)
        except KeyError: # If it fails, because there is no key
            super(DictList, self).__setitem__(key, value)
        except AttributeError: # If it fails because it is not a list
            super(DictList, self).__setitem__(key, [self[key], value])


def update_sets(sets_, covered):
    sets_ = [s.difference(covered) for s in sets_]
    return sets_

def run(parameters, original_sets_):
    for starting_size, starting_index, seed in parameters:

        ### Initialisation ###
        sets_ = copy.deepcopy(original_sets_)
        sizes = DictList()

        for i, s in enumerate(sets_):
            sizes[len(s)] = i

        random.seed(seed)
        universe = {i for i in range(94)}
        covered = set()
        indices = []

        ### First Iteration ###

        if starting_index:
            index = sizes[starting_size][starting_index]
            covered |= sets_[index]
            sets_.pop(index)
            indices.append(index)
        else:
            index = random.choice(sizes[starting_size])
            covered |= sets_[index]
            sets_.pop(index)
            indices.append(index)

        sets_ = update_sets(sets_, covered)

        ### Main body ###
        while covered != universe:

            sizes = DictList()

            for i, s in enumerate(sets_):
                sizes[len(s)] = i

            max_coverage = max(sizes.keys())

            if isinstance(sizes[max_coverage], list):
                index = random.choice(sizes[max_coverage])
            else:
                index = sizes[max_coverage]

            subset = sets_[index]

            covered |= subset

            indices.append(index)

            sets_.pop(index)

            sets_ = update_sets(sets_, covered)
            
        with open(f'set_cover_output/output_starting_size_{starting_size}_starting_index_{starting_index}_seed_{seed}.txt', 'w') as f:
            data = [starting_index, len(indices), *indices]
                    
            f.write(",".join([str(d) for d in data]))

search_parameters = [(51, 36),
 (50, 97),
 (49, 35),
 (48, 37),
 (47, 25),
 (46, 42),
 (45, 32),
 (44, 4),
 (43, 6),
 (42, 83),
 (41, 54),
 (40, 60),
 (39, 24),
 (38, 77),
 (37, 34),
 (36, 51),
 (35, 83),
 (34, 79),
 (33, 72),
 (32, 68),
 (31, 221),
 (30, 232),
 (29, 415),
 (28, 200),
 (27, 210),
 (26, 189),
 (25, 214),
 (24, 259),
 (23, 246),
 (22, 282)]

if __name__ == "__main__":

    df = pd.read_csv("sets_per_strategy.csv")

    original_sets_ = []
    for i, row in df.iterrows():
        original_sets_.append(set(row.dropna().values))

    seeds = range(10)

    for size, max_index in  tqdm.tqdm(search_parameters[:1]):
        parameters = (itertools.product([size], range(max_index), seeds))
        run(parameters, original_sets_)