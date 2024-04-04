function ErrorGrid(sdim);

error_probabilities = -5:0.5:-1;
num_error_values = 9;

cost_values = 0.1:0.1:1;

num_cost_values = 10;

N = 100;
b = 1;
numberIterations = 10 ^ 7;
seed = 1;
beta = 1;

folder = "ErrorGrid/dimension_";

for i=1:num_cost_values
    for j=1:num_error_values
        c = cost_values(i);
        starting_resident = zeros(1, sdim);
        errorprobability = 10 ^ error_probabilities(j);
        filename = append(folder, num2str(sdim), "_error_", num2str(errorprobability), "_cost_", num2str(c));
        [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
    end
end
end
