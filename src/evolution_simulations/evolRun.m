% Script to run the evolutionary process for a given set of parameters
function evolRun();

sdims = [2; 4];

N = 100;
b = 1;
numberIterations = 10 ^ 3;
seed = 1;
errorprobabilities = [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4];

beta = 1;
c = 0.3;

filenames = ["ErrorSim/reactive_one_error_"; "ErrorSim/reactive_two_error_"];

for i=1:2
    for j=1:7 
    sdim = sdims(i);
    errorprobability = errorprobabilities(j);
    starting_resident = zeros(1, sdim);
    filename = append(filenames(i), num2str(errorprobability));
    [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
    end
end
end
