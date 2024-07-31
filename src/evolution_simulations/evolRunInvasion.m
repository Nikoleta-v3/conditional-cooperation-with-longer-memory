function evolRunInvasion();

sdim = 4;

N = 100;
b = 1;
numberIterations = 10 ^ 3;
seed = 2;
beta = 1;
c = 0.3;
errorprobability = 0;
% counting = false;

reactive_two_resident = [0.999985	0.160136	0.553336	0.035629];

counting_two_resident = [0.998204	0.466715	0.466715	0.317605];
counting = false;

starting_resident = reactive_two_resident;


filename = append("../Understanding_counting/reactive_again");
[xDat]=Invasion(N, c, b, beta, numberIterations, starting_resident, seed, sdim, counting, errorprobability, filename);

end
