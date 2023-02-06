%% This function generates the likelihood of each model/paramters
% This function contains the actual algorythms

function [diffQ, Q, V, PEhist] = Model_timeseries_scan(params,s,a,r,model)

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
ntr = 25;
ncond   = length(unique(s));
Q          = zeros(ncond,ntr);        % Initial option values (all Models) as a function of conditio ("s")
V          = zeros(ncond,ntr);
trialc     = zeros(ncond);



for i = 1:length(a)
    % update trial counter
    trialc(s(i)) =  trialc(s(i)) +1;
    if mod(i,25) == 1
        Q(s(i),trialc(s(i))) = a(i);
    end
    
    if ~isnan(a(i)) && trialc(s(i)) <25
        if model == 1 % Fixed beta and lr
            %% Basics Q-learning
            % Prediction error
            deltaI                      = r(i) - Q(s(i),trialc(s(i)));                                        % the prediction error 
            Q(s(i),trialc(s(i))+1) = Q(s(i),trialc(s(i))) + lr1 * deltaI;                                     % the delta rule 
            
        elseif model==2
            deltaI                             =  r(i) - Q(s(i),trialc(s(i)));
            if deltaI <0
                Q(s(i),trialc(s(i))+1) = Q(s(i),trialc(s(i))) + lr1_neg * deltaI;
            elseif deltaI >0
                Q(s(i),trialc(s(i))+1) = Q(s(i),trialc(s(i))) + lr1_pos * deltaI;
            end
            
        elseif model == 3
            deltaI                      =  abs(r(i)-w) - Q(s(i),trialc(s(i)));                           % the prediction error incl recalibration
            Q(s(i),trialc(s(i))+1) = Q(s(i),trialc(s(i))) + lr1 * deltaI;                                % the delta rule 
            

        end
        
        
    elseif isnan(a(i))
        Q(s(i),trialc(s(i))+1)        = Q(s(i),trialc(s(i)))   ;
    end
    diffQ(s(i),trialc(s(i))) = Q(s(i),trialc(s(i))) - a(i);
    PEhist(i) = deltaI;
end

