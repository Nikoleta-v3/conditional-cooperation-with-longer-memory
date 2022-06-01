function data=cooperationRates(sdim, epsilon);

strategies = dec2bin(0:2 ^ sdim - 1, sdim) - '0'; 

for i=progress(1:sdim ^ 2)
    data = zeros(sdim ^ 2, 7 + (sdim) * 2);
    for j=(1:sdim^2)
        p = (1 - epsilon) * strategies(i, :) + epsilon * (1 - strategies(i, :));
        q = (1 - epsilon) * strategies(j, :) + epsilon * (1 - strategies(j, :));
        v = stationary(p, q, sdim);
        xCC=sum(v(1:4:end)); xCD=sum(v(2:4:end)); 
        xDC=sum(v(3:4:end)); xDD=sum(v(4:4:end));
        data(j,:) = [i, j, strategies(i, :), strategies(j, :), epsilon, xCC, xCD, xDC, xDD];
    end
    filename = "cooperation_data/error_" + epsilon + "sdim_" + sdim + "_s" + i + ".csv";
    dlmwrite(filename, data, 'precision', 9);
end