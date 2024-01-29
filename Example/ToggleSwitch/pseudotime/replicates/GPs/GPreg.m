
function [meanVet,stanDev] = GPreg(DATA,temps,dt,oStep,T)

T = T-1;
x = temps';

for k=1:size(DATA,2)

    y = DATA(:,k);                  % concentrations   
    z = linspace(0, T, (T/dt)+1)';  % test points (to predict)
    n = length(y);                  % number of training points    
    
    covfunc = @covSEard;            % covariance function    
    hyp.cov = [0; 0];               % hyperparameters covariance function    
    meanfunc = [];                  % mean function is chosen to be zero    
    likfunc = @likGauss;            % likelihood function
    sn      = 2.5;
    hyp.lik = log(sn);
    
    % Optimization over the hyperparameters (by minimizing the negative log marginal likelihood w.r.t. the hyperparameters)
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);  
    
    % Predictive distribution
    [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, z);
    
    meanVet(k,:) = m;
    stanDev(k,:) = sqrt(s(1:(oStep/dt):end)); % to use in the likelihood for MCMC
    
end

for i = 1:size(meanVet,1)
    for j = 1:size(meanVet,2)
        if meanVet(i,j) < 0
            meanVet(i,j) = 0;
        end
    end
end