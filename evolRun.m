function evolRunReactive();

starting_resident = [0, 0, 0, 0];
N = 100;
beta = 1;
b=1;
numberIterations = 10 ^ 7;

seeds = linspace(0, 10, 11);
c = 0.3;
parfor (i = 1:11)
    seed = seeds(i)
    [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed);
end
end