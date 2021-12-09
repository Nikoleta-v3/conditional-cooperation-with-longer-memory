function evolRun();

starting_resident = [0, 0, 0];
N = 100;
beta = 1;
numberIterations = 10 ^ 3;

cs = linspace(0, 1, 11);

parfor (i = 1:11)
    c = cs(i);
    [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident);
end
end