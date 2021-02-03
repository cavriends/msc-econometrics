function [likelihood] = NegMarkovLogLikelihood(parameters, y)

    p11   = parameters(1,1);
    p22   = parameters(2,1);
    mu_1  = parameters(3,1);
    mu_2  = parameters(4,1);
    sigma_1 = parameters(5,1);
    sigma_2 = parameters(6,1);

    forecasted_state_probs = HamiltonFilter(y, p11, p22, [mu_1,mu_2], [sigma_1,sigma_2]);
    
    likelihood = -sum(log(forecasted_state_probs(:,1).'.*normpdf(y, mu_1, sigma_1)+forecasted_state_probs(:,2).'.*normpdf(y, mu_2, sigma_2)));
end