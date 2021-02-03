function [forecasted_state_probs] = HamiltonFilter(y, p11, p22, mu, sigma)
   
    forecasted_state_probs = zeros(length(y),2);
    forecasted_state_probs(1,1) = 1; %(1-p22)/(2-p11-p22);
    forecasted_state_probs(1,2) = 0; %(1-p11)/(2-p11-p22);
    
    updated_state_probs = zeros(length(y),2);
    complete_state_likelihood = zeros(length(y),2);
    P = [p11,1-p22;1-p11,p22];
    
    for i = 1:length(y)
        
        if i == 1
            
            for s = 1:2
                complete_state_likelihood(i,s) = normpdf(y(i), mu(s), sigma(s))*forecasted_state_probs(i,s);
            end
            
        else
            
            for s = 1:2
                forecasted_state_probs(i,s) = P(s,1)*updated_state_probs(i-1,1)+P(s,2)*updated_state_probs(i-1,2);
                complete_state_likelihood(i,s) = normpdf(y(i), mu(s), sigma(s))*forecasted_state_probs(i,s);     
            end
           
        end
        
        likelihood_contrib = sum(complete_state_likelihood(i,:));
        updated_state_probs(i,:) = complete_state_likelihood(i,:) ./ likelihood_contrib;
        
    end
    
end