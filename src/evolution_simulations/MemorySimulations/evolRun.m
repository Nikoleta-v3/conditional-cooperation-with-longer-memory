% Script to run the evolutionary process for a given set of parameters
function evolRun(beta);

N = 100;
b = 1;
numberIterations = 10 ^ 7;
seed = 1;
memory = 2;
counting = true;

costs = [0:0.1:1];
errorprobability = 0;

parfor index=1:11
    c = costs(index);
    filename = append("data/memory_two_counting_beta_", num2str(beta) ,"_cost_", num2str(c));
    [xDat]=evolSimulationMemoryN(N, c, b, beta, numberIterations, seed, memory, errorprobability, filename, counting);
end
end
