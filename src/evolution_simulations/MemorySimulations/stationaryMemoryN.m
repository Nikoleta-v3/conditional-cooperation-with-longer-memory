function v=stationaryMemoryN(p1, p2, n, errorprobability);

if n == 1
    M = MemoryOneM(p1, p2);
elseif n == 2
    M = MemoryTwoM(p1, p2);
elseif n == 3
    M = MemoryThreeM(p1, p2);
end

ssp0 = null(eye(size(M))- M.');
v = ssp0./sum(ssp0);

v = v';
end