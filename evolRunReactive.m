function evolRunReactive();

starting_resident = [0, 0, 0, 0];
N = 100;
beta = 1;
b=1;
numberIterations = 10 ^ 3;

seeds = linspace(1, 10, 10);
c = 0.3;
parfor (i = 1:10)
    seed = seeds(i)
    [xDat]=evolReactiveSimulation(N, c, b, beta, numberIterations, starting_resident, seed);
end
end
