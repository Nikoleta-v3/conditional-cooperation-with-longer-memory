function v=stationary(p1, p2, sdim);

if sdim == 2
    M = oneBitM(p1, p2);
elseif sdim == 4
    M = twoBitsM(p1, p2);
elseif sdim == 8
    M = threeBitsM(p1, p2);
elseif sdim == 16
    M = fourBitsM(p1, p2);
end

ssp0 = null(eye(size(M))- M.');
v = ssp0./sum(ssp0);

v = v';
end