## Numerical Exploration of Nash

The project focuses on characterizing Nash Equilibria within the space of
reactive strategies. Even though we have analytical expressions for the
strategies that are cooperative Nash Equilibria, we also verified our results by
running numerical experiments.

The numerical evaluation can be summarized as follows:

1. For a given $n$, create the space of pure memory-$n$ strategies.

2. For the same $n$, get a random reactive strategy within the space of
   feasible strategies.

3. Calculate the payoff that the random strategy received against itself.

4. Calculate all the payoffs of pure strategies against the random strategy.

If the payoff of step 3 is higher or equal to all the payoffs of step 4, then
the strategy is Nash.

To run the scripts in this folder, you need to install one package that was
created while working on the project (see `requirements.txt`). The package is
called `repeated-play`, and it can be installed by simply running:

```
pip install repeated_play
```