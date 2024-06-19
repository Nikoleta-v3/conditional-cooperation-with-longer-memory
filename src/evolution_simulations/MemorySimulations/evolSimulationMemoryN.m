% Simulation Code for the Evolutionary Dynamics of Reactive Strategies
function [xDat]=evolSimulationMemoryN(N, c, b, beta, numberIterations, seed, memory, errorprobability, filename, counting);
rng(seed)

if memory == 1
    sdim = 4;
elseif memory == 2
    sdim = 16;
elseif memory == 3
    sdim = 64;
end

starting_resident = zeros(1, sdim);

%% Preparations for the output
Data=['c=',num2str(c),'; b=',num2str(b),'; N=',num2str(N), '; beta=',num2str(beta), '; nIt=',num2str(numberIterations)];
AvCoop=0; AvPay=0; Res=starting_resident;


%% Initialization
xDat=zeros(numberIterations/100, sdim + 2);
xDat(1,:)=[Res, 0, 0];

u = repmat([b - c, -c, b, 0], 1, sdim / 4);

%% Running the evolutionary process
j = 2;
for t = 1:numberIterations
    Mut=rand(1, sdim);
    if counting == 1
        if memory == 1
        Mut(3) = Mut(2);
        elseif memory == 2
            % Mut(1) 4 cooperations
            % Mut(2) 3 cooperations
            Mut(3) = Mut(2);
            % Mut(4) 2 cooperations
            Mut(5) = Mut(2);
            Mut(6) = Mut(4);
            Mut(7) = Mut(4);
            % Mut(8) 1 cooperation
            Mut(9) = Mut(2);
            Mut(10) = Mut(4);
            Mut(11) = Mut(4);
            Mut(12) = Mut(8);
            Mut(13) = Mut(4);
            Mut(14) = Mut(8);
            Mut(15) = Mut(8);
            % Mut(16) 0 cooperations
        end
    end
    disp(Mut)
    [phi, coopM, piM]=calcPhi(Mut, Res, N, u, beta, memory, errorprobability);
    if rand(1) < phi
        Res=Mut; xDat(j,:)=[Res, t, coopM]; j=j+1;
    end
end

dlmwrite(filename + ".csv", xDat, 'precision', 9);
writematrix(Data, filename + ".txt");
end

function [phi, coopMM, piMM]=calcPhi(Mut, Res, N, u, beta, memory, errorprobability);
%% Calculating the fixation probability

vMM=stationaryMemoryN(Mut, Mut, memory, errorprobability);
vMR=stationaryMemoryN(Mut, Res, memory, errorprobability);
vRM=stationaryMemoryN(Res, Mut, memory, errorprobability);
vRR=stationaryMemoryN(Res, Res, memory, errorprobability);

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
