function evolRunInvading();

sdim = 4;

N = 100;
b = 1;
numberIterations = 10 ^ 3;
seed = 1;
beta = 1;
c = 0.5;
errorprobability = 0;

invaders = [[1, .05, .6, .3]; [1, .6, .05, .3]; [1, .05, .6, .1]; [1, .6, .05, .1]];

for i=progress(1:4)
    starting_invader = invaders(i, :);
    p1 = starting_invader(1);
    p2 = starting_invader(2);
    p3 = starting_invader(3);
    p4 = starting_invader(4);

    % filename = append("Invasion/sdim_2_invasions_resident_p_", num2str(p), '_q_', num2str(q));
    filename = append("Invading/sdim_4_invasions_by_p1_", num2str(p1), '_p2_', num2str(p2), '_p3_', num2str(p3), '_p4_', num2str(p4));

    [xDat]=Invading(N, c, b, beta, numberIterations, starting_invader, seed, sdim, errorprobability, filename);

end
end
