function evolRunInvasion();

sdim = 4;

N = 100;
b = 1;
numberIterations = 10 ^ 4;
seed = 2;
beta = 1;
c = 0.5;
errorprobability = 0;
% counting = false;

pplayer = [1 0.10 0.60 0.30];

pprimeplayer = [1 0.60 0.10 0.30];
counting = false;

starting_resident = pplayer;
p1 = starting_resident(1);
p2 = starting_resident(2);
p3 = starting_resident(3);
p4 = starting_resident(4);


filename = append("Invasion/sdim_", num2str(sdim), '_invasions_resident_p1_', num2str(p1), '_p2_', num2str(p2), '_p3_', num2str(p3), '_p4_', num2str(p4));
[xDat]=Invasion(N, c, b, beta, numberIterations, starting_resident, seed, sdim, counting, errorprobability, filename);

end
