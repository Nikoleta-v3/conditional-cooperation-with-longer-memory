% Script to run the evolutionary process for a given set of parameters
function evolRun();

N = 100;
b = 1;
numberIterations = 10 ^ 7;
seed = 1;
memory = 1;

beta = 1; % to change
costs = [0:0.1:1];
errorprobability = 0;

for index=1:11;
    c = costs(index);
    filename = append("data/memory_one_beta_1_cost_", num2str(c));
    [xDat]=evolSimulationMemoryN(N, c, b, beta, numberIterations, seed, memory, errorprobability, filename);
end
end
