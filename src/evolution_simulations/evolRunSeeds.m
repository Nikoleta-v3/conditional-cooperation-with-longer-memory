% Script to run the evolutionary process for a range of seeds
function evolRunSeeds();

sdim = 4;
starting_resident = zeros(1, sdim);
N = 100;
b = 1;
numberIterations = 10 ^ 7;

beta = 1;
c = .5;
seeds = linspace(1, 10, 10);
errorprobability = 0;
% seed =1;
% 
% filename = append("test_with_mutants_take_one_", num2str(seed), ".csv");
% [xDat]=evolSimulationGrid(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);

parfor (i = 1:10)
    seed = seeds(i);
    filename = append("Grid/test_with_mutants_take_one_", num2str(seed));
    [xDat]=evolSimulationGrid(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
end