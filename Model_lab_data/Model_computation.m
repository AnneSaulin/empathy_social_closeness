%% This function generates the likelihood of each model/paramters
% This function contains the actual algorythms

function ll = Model_computation(params,s,a,r,model)

% "s"=state, corresponds to "con"
% "a"=action, corresponds to "cho"
% "r"=reward, corresponds to "out"
% "c"=counterfactual, corresponds to "cou"

%% Parameters
switch model
    case 1
        lr1     = params(1);
        
    case 2
        lr1_pos     = params(1);
        lr1_neg     = params(2);
        
    case 3
        w       = params(1);
        lr1     = params(2);        
    
end

%% Hidden variables
ncond   = length(unique(s));
Q          = a([1 1+25 1+25*2 1+25*3]);        % Initial option values (all Models) is selected as the first rating in each block
V          = a([1 1+25 1+25*2 1+25*3]);        % Initial Context value is selected as the first rating in each block
% V           = zeros(ncond,1);        % Initial Context values (Models 3) as a function of conditio ("s")
Ccount  = zeros(ncond,1);
ll        = 0;                                                                             % distance between prediction and actual


for n = 1:length(a)
    
    Ccount(s(n)) = Ccount(s(n)) +1;
    
    if model == 1
        ll   = ll - log(normpdf(a(n),Q(s(n)),0.4));%(a(n)-Q(s(n))).^2;
        deltaI      = r(n) - Q(s(n));                                        % the prediction error 
        Q(s(n))   = Q(s(n)) + lr1 * deltaI;                                  % the delta rule 
        
    elseif model==2
        ll   = ll - log(normpdf(a(n),Q(s(n)),0.4));%(a(n)-Q(s(n))).^2;
        deltaI          = r(n) - Q(s(n));
        if deltaI <0  % control condition
            Q(s(n))    = Q(s(n)) + lr1_neg * deltaI;                         % the delta rule for non-reinforced trials
        elseif deltaI >0 % decay condition
            Q(s(n))    = Q(s(n)) + lr1_pos * deltaI;                         % the delta rule for reinforced trials
        end
        
    elseif model == 3
        ll   = ll - log(normpdf(a(n),Q(s(n)),0.4));%(a(n)-Q(s(n))).^2;
        deltaI      = abs(r(n)-w) - Q(s(n));                                    % the prediction error including recalibration
        Q(s(n))   = Q(s(n)) + lr1 * deltaI;                                     % the delta rule for the factual choice        
   
end

