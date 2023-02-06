% This function find the best fitting model/patamters
% This script prepare the data to be fed to computational model to obtain the LL and the LPP
% This script assumes all the data and codes are in the same folder
% This script assumes the subjects have performed only one session
% Using preservation rate !! https://www.jneurosci.org/content/38/48/10338
% The idea is: if people do keep closeness in mind when there is no reinforce
%, the PE might be smaller.

close all
clear
clc
rng('shuffle')
%%% set correct directory to where the files and functions are stored
cur_dir         = 'C:\~';
project_name    = '';
findnm          = strfind(pwd,project_name);
root_dir        = fullfile(cur_dir,project_name);
n_exp =1;

for k_exp = 1:n_exp
    
    subjects        = 46; % 
    n_sub           = 46;
    nTrial          = 25;
    nfpm            = [1 2 2];           % number of free parameters
    n_model       = length(nfpm);
    npmmax       = max(nfpm);
    
    % Pre-Allocate
    parameters      = NaN(n_exp,n_sub,n_model,npmmax);    parametersLPP   = NaN(n_exp,n_sub,n_model,npmmax);
    ll              = NaN(n_exp,n_sub,n_model);           LPP             = NaN(n_exp,n_sub,n_model);
    snll_learn      = NaN(n_exp,n_sub,n_model);
    snll_post       = NaN(n_exp,n_sub,n_model);
    loglik            = NaN(n_exp,n_sub,n_model);
    AIC              = NaN(n_exp,n_sub,n_model);
    BIC              = NaN(n_exp,n_sub,n_model);
    check_conv = NaN(n_exp,n_sub,1);
    
    % Optimization parameters
    % options = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000); % These increase the number of iterations to ensure the convergence
    options         = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000); % These increase the number of iterations to ensure the convergence
    lb = [0 0 0];        LBx = [0 0 0];
    ub = [1 1 1];       UBx = [1 1 1];
    ddb = ub - lb;
    
    %% preallocation
    % Condition: 1 = contral_B1; 2 = contral_B2; 3 = decay_B1;4 = decay_B2
    rating = NaN(n_exp,n_sub,4,nTrial); % subjcect, conditions, tirals.
    info    = NaN(n_exp,n_sub,4,nTrial); % subjcect, conditions, tirals.
end

for k_exp =1:n_exp
    clear proxempcon1scan46 proxempcon2scan46 reinforcerempcon1scan46 reinforcerempcon2scan46 proxemptreat1scan46 proxemptreat2scan46 reinforceremptreat1scan46 reinforceremptreat2scan46 ...
        
       
    data_dir        = fullfile(root_dir,'dat_fMRI_sample.mat');
    load(data_dir);
    
    [~, id] = sort(reinforcerempcon1scan46.trial_no);
reinforcerempcon1scan46 = reinforcerempcon1scan46(id,:);

[~, id] = sort(reinforcerempcon2scan46.trial_no);
reinforcerempcon2scan46 = reinforcerempcon2scan46(id,:);

[~, id] = sort(reinforceremptreat1scan46.trial_no);
reinforceremptreat1scan46 = reinforceremptreat1scan46(id,:);

[~, id] = sort(reinforceremptreat2scan46.trial_no);
reinforceremptreat2scan46 = reinforceremptreat2scan46(id,:);
   
    
    for k_sub       = 1:n_sub % subject loop
        
        rating(k_exp,k_sub,1,:) = proxempcon1scan46{1:end,k_sub+2}*-1+100;
        rating(k_exp,k_sub,2,:) = proxempcon2scan46{1:end,k_sub+2}*-1+100;
        rating(k_exp,k_sub,3,:) = proxemptreat1scan46{1:end,k_sub+2}*-1+100;
        rating(k_exp,k_sub,4,:) = proxemptreat2scan46{1:end,k_sub+2}*-1+100;
        
        % proxempcon1scan.Properties.VariableNames{k_sub+2} %prints name of
        % current vp
        
        % check outliers
        M_rating = mean(squeeze(rating(k_exp,k_sub,:)));
        S_rating = std(squeeze(rating(k_exp,k_sub,:)));
        for k_cond = 1:4
            for k_trial = 1:nTrial
                if  rating(k_exp,k_sub,k_cond,k_trial) > (M_rating-3*S_rating) && rating(k_exp,k_sub,k_cond,k_trial) < (M_rating+3*S_rating)
                else
                    rating(k_exp,k_sub,k_cond,k_trial) = NaN;
                end
                
                if  k_trial == 1 && isnan(rating(k_exp,k_sub,k_cond,1))
                    rating(k_exp,k_sub,k_cond,k_trial) = rating(k_exp,k_sub,k_cond,2);
                elseif  isnan(rating(k_exp,k_sub,k_cond,k_trial))
                    rating(k_exp,k_sub,k_cond,k_trial) = rating(k_exp,k_sub,k_cond,k_trial-1);
                end
            end
        end
        info(k_exp,k_sub,1,:)    = reinforcerempcon1scan46{1:end,k_sub+2};
        info(k_exp,k_sub,2,:)    = reinforcerempcon2scan46{1:end,k_sub+2};
        info(k_exp,k_sub,3,:)    = reinforceremptreat1scan46{1:end,k_sub+2};
        info(k_exp,k_sub,4,:)    = reinforceremptreat2scan46{1:end,k_sub+2};
        
        con              = [ones(nTrial,1);ones(nTrial,1)*2;ones(nTrial,1)*3;ones(nTrial,1)*4]'; %"con": 1 = contral_B1; 2 = contral_B2; 3 = decay_B1;4 = decay_B2;
        cho              = [squeeze(rating(k_exp,k_sub,1,:));squeeze(rating(k_exp,k_sub,2,:));squeeze(rating(k_exp,k_sub,3,:));squeeze(rating(k_exp,k_sub,4,:))]./100;   % calibrated rating (range between -1 and 1)
        out              = [squeeze(info(k_exp,k_sub,1,:));squeeze(info(k_exp,k_sub,2,:));squeeze(info(k_exp,k_sub,3,:));squeeze(info(k_exp,k_sub,4,:))];
        
        % estimate models
        for k_model = 1:n_model
            % This part requires the Matlab Optimization toolbox
            % prepare multiple starting points for estimation
            n_rep           = 5;
            parameters_rep  = NaN(n_rep,nfpm(k_model));     parametersLPP_rep  = NaN(n_rep,nfpm(k_model));
            ll_rep          = NaN(n_rep,1);          LPP_rep                  = NaN(n_rep,1);
            FminHess        = NaN(n_rep,nfpm(k_model),nfpm(k_model));
            
            for k_rep = 1:n_rep
                % set random staring points & params bounds
                x0 = lb + rand(1,3).*ddb;
                x0 = x0(1:nfpm(k_model));
                LB = LBx(1:nfpm(k_model));
                UB = UBx(1:nfpm(k_model));
                
                
                % run ML and MAP estimations
                [parameters_rep(k_rep,1:nfpm(k_model)),ll_rep(k_rep),~,~,~] = fmincon(@(x) Model_computation_scan(x,con,cho,out,k_model),x0,[],[],[],[],LB,UB,[],options);
                [parametersLPP_rep(k_rep,1:nfpm(k_model)),LPP_rep(k_rep),~,~,~,~,FminHess(k_rep,:,:)] = fmincon(@(x) Model_ParametersPriors_scan(x,con,cho,out,k_model),x0,[],[],[],[],LB,UB,[],options);
            end
            
            [~,pos]                                 = min(ll_rep);
            ll(k_exp,k_sub,k_model)  = min(ll_rep);
            parameters(k_exp,k_sub,k_model,1:nfpm(k_model))   = parameters_rep(pos(1),1:nfpm(k_model));
            
            [~,pos_LPP]                          = min(LPP_rep);
            
            parametersLPP(k_exp,k_sub,k_model,1:nfpm(k_model))   = parametersLPP_rep(pos_LPP(1),1:nfpm(k_model));
            LPP(k_exp,k_sub,k_model)  = LPP_rep(pos_LPP(1),:) - nfpm(k_model)*log(2*pi)/2 + log(abs(det(squeeze(FminHess(pos_LPP(1),:,:)))))/2;
            BIC(k_exp,k_sub,k_model)  = -2*-squeeze(ll(k_exp,k_sub,k_model)) + nfpm(k_model)*log(2*25*2); % l2 is already positive
            AIC(k_exp,k_sub,k_model)  = -2*-squeeze(ll(k_exp,k_sub,k_model)) + 2*nfpm(k_model);
            
            % Re-Calculate likelihood from actual choices
            [diffQ1, Q1, V1]    = Model_timeseries_scan_PE2(parameters(k_exp,k_sub,k_model,1:nfpm(k_model)), con(:), cho(:), out(:),k_model);
            err_learn(k_exp,k_sub,k_model)                   = nansum((Q1(:)-cho(:)).^2);
            Q1hist(k_exp,k_sub,k_model,:,:) = Q1;
            V1hist(k_exp,k_sub,k_model,:,:) = V1;
            %
            [diffQ2, Q2, V2, PEhist]    = Model_timeseries_scan(parametersLPP(k_exp,k_sub,k_model,1:nfpm(k_model)), con(:), cho(:), out(:),k_model);
            err_learn_LPP(k_exp,k_sub,k_model)                = -nansum(log(normpdf(diffQ2(:)))); % Using Q2-cho is wrong, because they are not aligned.
            sum_error_replication(k_exp,k_sub,k_model,:)      = sum(diffQ2.^2,2);
            Q2hist(k_exp,k_sub,k_model,:,:)                   = Q2;
            V2hist(k_exp,k_sub,k_model,:,:)                   = V2;
            PEhist_all(k_exp,k_sub,k_model,:)                 = PEhist;
        end
    end
    %         snll_learn(k_sub,k_model)   = -nansum(log(lik_learn));
    %         snll_post(k_sub,k_model)    = -nansum(log(lik_post));
    %         loglik(k_sub,k_model)       = snll_learn(k_sub,k_model) + snll_post(k_sub,k_model);
    %     end
    % end
    %
    
    %%  save('model1_2019_11_14_XP2','parameters','parametersLPP','ll','LPP','con','cho','BIC','AIC');
    %     save('model_scan_2021_0415','parameters','parametersLPP','ll','LPP','con','cho','BIC','AIC', ...
    %         'err_learn','err_learn_LPP','rating','Q1hist','Q2hist','sum_error_replication');
    save('model_scan','parameters','parametersLPP','ll','LPP','con','cho','BIC','AIC', ...
        'err_learn','err_learn_LPP','rating','sum_error_replication','Q2hist','PEhist_all'); 
end


%% get BIC values
for k = 1:4
    BICcompare(k) = sum(sum((squeeze(BIC(:,:,k)))));
    BICcompare_exp(1,k) = sum((squeeze(BIC(1,:,k))));
    %BICcompare_exp(2,k) = sum((squeeze(BIC(2,:,k))));
end

%% conduct Bayesian model comparison
for k_exp = 1:2
    [postBMCbic,outBMCbic]  = VBA_groupBMC(-squeeze(BIC(k_exp,:,:))');
    BMCbic_output.post(k_exp)    = postBMCbic;
    BMCbic_output.out(k_exp)        = outBMCbic;
    %
    % % % simplest models
    % [postBMClpp,outBMClpp]  = VBA_groupBMC(-squeeze(LPP(k_exp,:,[1 3 4 7]))');
    % % extension of model1 and model7
    % [postBMClpp,outBMClpp]  = VBA_groupBMC(-squeeze(LPP(k_exp,:,[1 2 7 8]))');
    % all models
%     [postBMClpp,outBMClpp]  = VBA_groupBMC(-squeeze(LPP(k_exp,:,1:3))');
%     
%     BMClpp_output.post(k_exp)    = postBMClpp;
%     BMClpp_output.out(k_exp)     = outBMClpp;
end
for k_con = 1:3
    estimatedDATA_emp{k_con} = squeeze(Q2hist(1,:,3,k_con,:))';
    %estimatedDATA_rec{k_con} = squeeze(Q2hist(2,:,1,k_con,:))';
end
save('model_scan','parameters','parametersLPP','ll','LPP','con','cho','BIC','AIC', ...
    'err_learn','err_learn_LPP','rating','PEhist_all', 'Q1hist','Q2hist','sum_error_replication','BMClpp_output',...
    'estimatedDATA_emp');

