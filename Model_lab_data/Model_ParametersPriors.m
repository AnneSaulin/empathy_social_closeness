%% This function attribute probability to parameters
% it is used to calculate the LPP
% the priors are taken from Daw et al. Neuron 2011
function [post]= Model_ParametersPriors(params,s,a,r,model)
% log prior of parameters
switch model
    case 1
        lr1     = params(1);                                                             % policy or factual learning rate
        
        plr1    = log(betapdf(lr1,1.1,1.1));
        p       = [plr1];
        
    case 2
        lr1_pos       = params(1);
        lr1_neg       = params(2);
        plr1_pos    = log(betapdf(lr1_pos,1.1,1.1));
        plr1_neg   = log(betapdf(lr1_neg,1.1,1.1));
        p       = [plr1_pos plr1_neg];
    case 3
        w       = params(1);
        lr1     = params(2);                                                             % policy or factual learning rate
        
        plrw    = log(betapdf(w,1.1,1.1));
        plr1    = log(betapdf(lr1,1.1,1.1));
        p       = [plrw plr1 ];
    
   
end

p               = -sum(p);
l               = Model_computation(params,s,a,r,model);

post            = p + l;



