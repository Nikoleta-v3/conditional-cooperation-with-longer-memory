function evolRunInvasion();

sdim = 4;

N = 100;
b = 1;
numberIterations = 10 ^ 3;
seed = 1;
beta = 1;
c = 0.5;
errorprobability = 0;

residents = [[1, .05, .6, .3]; [1, .6, .05, .3]; [1, .05, .6, .1]; [1, .6, .05, .1]];

for i=1:4
    starting_resident = residents(i, :);
    p1 = starting_resident(1);
    p2 = starting_resident(2);
    p3 = starting_resident(3);
    p4 = starting_resident(4);

%     % filename = append("Invasion/sdim_2_invasions_resident_p_", num2str(p), '_q_', num2str(q));
    filename = append("Invasion/sdim_4_invasions_resident_p1_", num2str(p1), '_p2_', num2str(p2), '_p3_', num2str(p3), '_p4_', num2str(p4));
    disp(filename)
    disp(starting_resident)
% 
%     [xDat]=Invasion(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);

end
end
