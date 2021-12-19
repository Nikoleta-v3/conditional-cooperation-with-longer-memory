function evolRunThreeBits();

starting_resident = [0, 0, 0, 0];
N = 100;
beta = 1;
b=1;
numberIterations = 10 ^ 3;
seed = 1;

cs = linspace(1, 10, 11);

parfor (i = 1:11)
    c = cs(i)
    [xDat]=evolSimulationThreeBits(N, c, b, beta, numberIterations, starting_resident, seed);
end
end
