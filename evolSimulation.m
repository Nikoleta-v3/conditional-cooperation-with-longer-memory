function [xDat]=evolSimulation(N, c, b, beta, numberIterations, starting_resident, seed);
rng(seed)
%% Preparations for the output
Data=['c=',num2str(c),'; b=',num2str(b),'; N=',num2str(N), '; beta=',num2str(beta), '; nIt=',num2str(numberIterations)];
AvCoop=0; AvPay=0; Res=starting_resident;
filename = "data/matlab_two_bits_beta_" + beta + "_seed_" + seed + "_c_" + c;

%% Initialization
sdim=4;
xDat=zeros(numberIterations/100, 6);
xDat(1,:)=[Res, 0, 0];


u = [b - c, -c, b, 0, b - c, -c, b, 0, b - c, -c, b, 0, b - c, -c, b, 0];

%% Running the evolutionary process
j = 2;
for t = progress(1:numberIterations)
    Mut=rand(1, sdim);
    [phi, coopM, piM]=calcPhi(Mut, Res, N, u, beta);
    if rand(1) < phi
        Res=Mut; xDat(j,:)=[Res, t, coopM]; j=j+1;
    end
end

dlmwrite(filename + ".csv", xDat, 'precision', 9);
writematrix(Data, filename + ".txt");

AvCoop = mean(xDat(:,end-2));
AvPay = mean(xDat(:,end-1));
end

function [phi, coopMM, piMM]=calcPhi(Mut, Res, N, u, beta);
%% Calculating the fixation probability

vMM=stationary(Mut, Mut);
vMR=stationary(Mut, Res);
vRM=stationary(Res, Mut);
vRR=stationary(Res, Res);

piMM=vMM*u';
coopMM=vMM(1) + vMM(5) + vMM(9) + vMM(13) + vMM(2) + vMM(6) + vMM(10) + vMM(14);


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
