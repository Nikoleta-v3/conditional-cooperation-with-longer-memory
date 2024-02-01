% Script to run the evolutionary process for a range of seeds
function evolRunSeeds();

sdim = 4;
starting_resident = zeros(1, sdim);
N = 100;
b=1;
numberIterations = 10 ^ 7;

beta = 1;
c = .5;
seeds = linspace(1, 10, 10);

parfor (i = 1:10)
    seed = seeds(i);
    [xDat]=evolSimulationCounting(N, c, b, beta, numberIterations, starting_resident, seed, sdim);
end
end
