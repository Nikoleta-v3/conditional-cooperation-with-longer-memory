function ErrorGrid(sdim);

error_probabilities = -5:0.5:-1; cost_values = 0.1:0.1:1;
num_error_values = 9; num_cost_values = 10;

% parameters = zeros(num_error_values * num_cost_values, 2);
% 
% for i=1:num_error_values
%     for j=1:num_cost_values
%         index = (i - 1) * 10 + j;
%         parameters(index, 1) = 10 ^ error_probabilities(i);
%         parameters(index, 2) = cost_values(j);
%     end
% end
errorprobability = 10 ^ error_probabilities(1);
N = 100;
b = 1;
numberIterations = 10 ^ 3;
beta = 1;
starting_resident = zeros(1, sdim);

folder = "ErrorGrid/dimension_";

parfor i=1:num_cost_values
    seed = i;
    c = cost_values(i);
    filename = append(folder, num2str(sdim), "_error_", num2str(errorprobability), "_cost_", num2str(c));
    [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
end
end
