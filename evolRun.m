function evolRun();

starting_resident = [0, 0, 0, 0];
N = 100;
beta = 0.1;
b=1;
numberIterations = 10 ^ 7;

cs = linspace(0, 1, 11);
seed = 1;
parfor (i = 1:11)
    c = cs(i);
    [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed);
end
end
