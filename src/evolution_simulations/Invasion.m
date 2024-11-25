% Simulation Code for the Evolutionary Dynamics of Reactive Strategies
function [xDat]=Invasion(N, c, b, beta, numberIterations, starting_resident, seed, sdim, counting, errorprobability, filename);
rng(seed)
%% Preparations for the output
Data=['c=',num2str(c),'; b=',num2str(b),'; N=',num2str(N), '; beta=',num2str(beta), '; nIt=',num2str(numberIterations)];
AvCoop=0; AvPay=0; Res=starting_resident;

%% Initialization
xDat=zeros(numberIterations, sdim + 2);
xDat(1,:)=[Res, 0, 0];

u = repmat([b - c, -c, b, 0], 1, (sdim ^ 2) / 4);

%% Running the evolutionary process
j = 2;
for i = progress(1:numberIterations)
    trials = 1;
    notInvaded = true;
    while notInvaded == true
        Mut=rand(1, sdim);
        if counting == true
            if sdim == 4
                Mut(3) = Mut(2);
            elseif sdim == 8
                Mut(3) = Mut(2);
                Mut(5) = Mut(2);
                Mut(6) = Mut(4);
             Mut(7) = Mut(4);
            end
        end
        
        trials = trials + 1;
        [phi, coopM, piM]=calcPhi(Mut, Res, N, u, beta, sdim, errorprobability);
        if rand(1) < phi
            xDat(j,:)=[Mut, trials, coopM]; j=j+1; notInvaded = false;
        end
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
