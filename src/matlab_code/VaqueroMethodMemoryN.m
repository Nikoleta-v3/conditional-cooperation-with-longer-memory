% This script is used to run the method for the strategies from Hilde 2017
% paper
function data=VaqueroMethodMemoryN(n, epsilon);

strategies = dec2bin(0:2 ^ (4 ^ n) - 1, 4 ^ n) - '0'; 

toCheck = [[1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];
 [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];
 [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1];
 [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1];
 [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1];
 [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1];
 [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1];
 [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1];
 [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1];
 [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];
 [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]];

indices = zeros(1, 11);

for i=progress(1:11)
    index = find(ismember(strategies, toCheck(i, :),'rows'));
    indices(i) = index;
end
    
for i=progress(1:11)
    data = zeros(2 ^ (4 ^ n), 7 + (4 ^ n) * 2);
    for j=(1:2 ^ (4 ^ n))
        p = (1 - epsilon) * toCheck(i, :) + epsilon * (1 - toCheck(i, :));
        q = (1 - epsilon) * strategies(j, :) + epsilon * (1 - strategies(j, :));
        v = stationaryMemoryN(p, q, n);
   
        xCC=sum(v(1:4:end)); xCD=sum(v(2:4:end));
        xDC=sum(v(3:4:end)); xDD=sum(v(4:4:end));
        
        data(j,:) = [indices(i), j, toCheck(i, :), strategies(j, :), epsilon, xCC, xCD, xDC, xDD];
    end
    filename = "cooperation_data/memone/error_" + epsilon + "_n_" + n + "_s" + i + ".csv";
    dlmwrite(filename, data, 'precision', 9);
end