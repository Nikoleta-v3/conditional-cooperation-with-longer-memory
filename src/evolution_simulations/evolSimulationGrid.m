% Simulation Code for the Evolutionary Dynamics of Reactive Strategies
function [xDat]=evolSimulationGrid(N, c, b, beta, numberIterations, starting_resident, seed, sdim, errorprobability, filename);
rng(seed)
%% Preparations for the output
Data=['c=',num2str(c),'; b=',num2str(b),'; N=',num2str(N), '; beta=',num2str(beta), '; nIt=',num2str(numberIterations)];
AvCoop=0; AvPay=0; Res=starting_resident;

%% Initialization
xDat=zeros(numberIterations/100, sdim + 2);
xDat(1,:)=[Res, 0, 0];

u = repmat([b - c, -c, b, 0], 1, (sdim ^ 2) / 4);

grid = (0:0.001:1);

%% Running the evolutionary process
j = 2;
for t = 1:numberIterations
    Mut = zeros(1, sdim);
    for l = 1:sdim
        Mut(l) = grid(randi([1 1001]));
    end
    [phi, coopM, piM]=calcPhi(Mut, Res, N, u, beta, sdim, errorprobability);
    if rand(1) < phi
        Res=Mut; xDat(j,:)=[Res, t, coopM]; j=j+1;
    end
end

dlmwrite(filename + ".csv", xDat, 'precision', 9);
writematrix(Data, filename + ".txt");
end

function [phi, coopMM, piMM]=calcPhi(Mut, Res, N, u, beta, sdim, errorprobability);
%% Calculating the fixation probability

vMM=stationary(Mut, Mut, sdim, errorprobability);
vMR=stationary(Mut, Res, sdim, errorprobability);
vRM=stationary(Res, Mut, sdim, errorprobability);
vRR=stationary(Res, Res, sdim, errorprobability);

disp(vMM)
disp(vMR)
disp(vRM)
disp(vRR)
disp("=====")

piMM=vMM*u';
coopMM= sum(vMM(2:4:end)) + sum(vMM(1:4:end));

piMR=vMR*u';
piRM=vRM*u';
piRR=vRR*u';


laplus=zeros(1, N-1); laminus=laplus;
for k=1:N-1
    piM = (k-1) / (N-1) * piMM + (N-k) / (N-1) * piMR;
    piR = k / (N-1) * piRM + (N-k-1) / (N-1) * piRR;

    laplus(k) = 1 / (1 + exp(-beta * (piM - piR)));
    laminus(k) = 1 / (1 + exp(-beta * (piR - piM)));
end

phi = 1 / (1 + sum(cumprod(laminus./laplus)));

end
