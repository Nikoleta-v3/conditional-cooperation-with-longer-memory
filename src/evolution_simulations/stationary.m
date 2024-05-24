function v=stationary(p1, p2, sdim, errorprobability);

strategy_p1 = (1 - errorprobability) * p1 + errorprobability * (1 - p1); 
strategy_p2 = (1 - errorprobability) * p2 + errorprobability * (1 - p2);

if sdim == 2
    M = reactiveOneM(strategy_p1, strategy_p2);
elseif sdim == 4
    M = reactiveTwoM(strategy_p1, strategy_p2);
elseif sdim == 8
    M = reactiveThreeM(strategy_p1, strategy_p2);
elseif sdim == 16
    M = reactiveFourM(strategy_p1, strategy_p2);
end

ssp0 = null(eye(size(M))- M.');
v = ssp0./sum(ssp0);

v = v';
end