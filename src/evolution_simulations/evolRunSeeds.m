% Script to run the evolutionary process for a range of seeds
function evolRunSeeds(seed);

sdim = 4;
starting_resident = zeros(1, sdim);
N = 100;
b = 1;
numberIterations = 10 ^ 5;

beta = 1;
c = .5;
% seeds = linspace(1, 10, 10);
errorprobability = 0;

filename = append("test_with_mutants_take_one_", num2str(seed), ".csv");
[xDat]=evolSimulationGrid(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);

% for (i = 1:2)
%     seed = seeds(i);
%     filename = append("test_with_mutants_take_one", num2str(seed), ".csv");
%     [xDat]=evolSimulationGrid(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
% %        [xDat]=evolSimulationCounting(N, c, b, beta, numberIterations, starting_resident, seed, sdim,  errorprobability, filename);
% end