function v=stationary(p1, p2, sdim);

if sdim == 2
    M = reactiveOneM(p1, p2);
elseif sdim == 4
    M = reactiveTwoM(p1, p2);
elseif sdim == 8
    M = reactiveThreeM(p1, p2);
elseif sdim == 16
    M = reactiveFourM(p1, p2);
end

ssp0 = null(eye(size(M))- M.');
v = ssp0./sum(ssp0);

v = v';
end