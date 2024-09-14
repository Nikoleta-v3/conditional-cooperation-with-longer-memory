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
        if memory == 2
            Mut(2) = Mut(5);
            Mut(3) = Mut(9);
            Mut(4) = Mut(13);
            Mut(7) = Mut(13);
            Mut(10) = Mut(13);
            Mut(8) = Mut(14);
            Mut(12) = Mut(15);
        elseif memory == 3
            probs =  Mut(1:7);
            Mut(1) = probs(1);
            Mut(2) = probs(2);
            Mut(3) = probs(2);
            Mut(4) = probs(3);
            Mut(5) = probs(2);
            Mut(6) = probs(3);
            Mut(7) = probs(3);
            Mut(8) = probs(4);
            Mut(9) = probs(2);
            Mut(10) = probs(3);
            Mut(11) = probs(3);
            Mut(12) = probs(4);
            Mut(13) = probs(3);
            Mut(14) = probs(4);
            Mut(15) = probs(4);
            Mut(16) = probs(5);
            Mut(17) = probs(2);
            Mut(18) = probs(3);
            Mut(19) = probs(3);
            Mut(20) = probs(4);
            Mut(21) = probs(3);
            Mut(22) = probs(4);
            Mut(23) = probs(4);
            Mut(24) = probs(5);
            Mut(25) = probs(3);
            Mut(26) = probs(4);
            Mut(27) = probs(4);
            Mut(28) = probs(5);
            Mut(29) = probs(4);
            Mut(30) = probs(5);
            Mut(31) = probs(5);
            Mut(32) = probs(6);
            Mut(33) = probs(2);
            Mut(34) = probs(3);
            Mut(35) = probs(3);
            Mut(36) = probs(4);
            Mut(37) = probs(3);
            Mut(38) = probs(4);
            Mut(39) = probs(4);
            Mut(40) = probs(5);
            Mut(41) = probs(3);
            Mut(42) = probs(4);
            Mut(43) = probs(4);
            Mut(44) = probs(5);
            Mut(45) = probs(4);
            Mut(46) = probs(5);
            Mut(47) = probs(5);
            Mut(48) = probs(6);
            Mut(49) = probs(3);
            Mut(50) = probs(4);
            Mut(51) = probs(4);
            Mut(52) = probs(5);
            Mut(53) = probs(4);
            Mut(54) = probs(5);
            Mut(55) = probs(5);
            Mut(56) = probs(6);
            Mut(57) = probs(4);
            Mut(58) = probs(5);
            Mut(59) = probs(5);
            Mut(60) = probs(6);
            Mut(61) = probs(5);
            Mut(62) = probs(6);
            Mut(63) = probs(6);
            Mut(64) = probs(7);
        end
    end
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
