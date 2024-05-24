% Simulation Code for the Evolutionary Dynamics of Reactive Strategies
function [xDat]=Invading(N, c, b, beta, numberIterations, invader, seed, sdim, errorprobability, filename);
rng(seed)
%% Preparations for the output
Data=['c=',num2str(c),'; b=',num2str(b),'; N=',num2str(N), '; beta=',num2str(beta), '; nIt=',num2str(numberIterations)];
AvCoop=0; AvPay=0; Mut=invader;

%% Initialization
xDat=zeros(numberIterations, sdim + 3);
xDat(1,:)=[Mut, 0, 0, 0];

u = repmat([b - c, -c, b, 0], 1, (sdim ^ 2) / 4);

%% Running the evolutionary process
j = 2;
for i = progress(1:numberIterations)
    Res=rand(1, sdim);
    [phi, coopR, piR]=calcPhi(Mut, Res, N, u, beta, sdim, errorprobability);
    if rand(1) < phi
        Invaded = 1;
    else
        Invaded = 0;
    end
    xDat(j,:)=[Res, i, coopR, Invaded]; j = j + 1;
end

dlmwrite(filename + ".csv", xDat, 'precision', 9);
writematrix(Data, filename + ".txt");
end

function [phi, coopRR, piRR]=calcPhi(Mut, Res, N, u, beta, sdim, errorprobability);
%% Calculating the fixation probability

vMM=stationary(Mut, Mut, sdim, errorprobability);
vMR=stationary(Mut, Res, sdim, errorprobability);
vRM=stationary(Res, Mut, sdim, errorprobability);
vRR=stationary(Res, Res, sdim, errorprobability);

piMM=vMM*u';
coopMM= sum(vMM(2:4:end)) + sum(vMM(1:4:end));
coopRR= sum(vRR(2:4:end)) + sum(vRR(1:4:end));

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
