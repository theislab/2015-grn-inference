
function y = loglik(x,mu,sigma)

% Compute Gaussian log-likelihood

y = zeros(1,length(x));

for k = 1:length(x)

    xn = (x(k) - mu(k)) ./ sigma;
    y(k) = exp(-0.5 * xn .^2) ./ (sqrt(2*pi) .* sigma);

end

y = sum(log(y));

