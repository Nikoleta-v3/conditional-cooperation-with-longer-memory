function data=VaqueroMethodMemoryN(n, epsilon);

strategies = dec2bin(0:2 ^ (4 ^ n) - 1, 4 ^ n) - '0'; 

for i=progress(1:2 ^ (4 ^ n))
    data = zeros(2 ^ (4 ^ n), 7 + (4 ^ n) * 2);
    for j=(1:2 ^ (4 ^ n))
        p = (1 - epsilon) * strategies(i, :) + epsilon * (1 - strategies(i, :));
        q = (1 - epsilon) * strategies(j, :) + epsilon * (1 - strategies(j, :));
        v = stationaryMemoryN(p, q, n);
   
        xCC=sum(v(1:4:end)); xCD=sum(v(2:4:end));
        xDC=sum(v(3:4:end)); xDD=sum(v(4:4:end));
        
        data(j,:) = [i, j, strategies(i, :), strategies(j, :), epsilon, xCC, xCD, xDC, xDD];
    end
    filename = "cooperation_data/error_" + epsilon + "_n_" + n + "_s" + i + ".csv";
    dlmwrite(filename, data, 'precision', 9);
end