% Script to run the evilutionary process for a given set of parameters
function evolRun();

sdim = 4;
starting_resident = zeros(1, sdim);
N = 100;
b=1;
numberIterations = 10 ^ 7;
seed = 1;

beta = 1;
c = .5;
seeds = linspace(0, 1, 11);

parfor (i = 1:11)
    seed = seeds(i);
    [xDat]=evolSimulationCounting(N, c, b, beta, numberIterations, starting_resident, seed, sdim);
end
end