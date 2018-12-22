%% test gradient descent in HMM model

% include adaptive learning rate

% seems that this works for MMs no hidden units, but cannot find the
% correct factor to mutate weights (as mutations should apply to the weights, not transition probabilities)

% Impose topological constraints

T = 100;
no_data_points = 4;

eps = 1e-8;
beta = 1;
nx = 2;
mx = 2;

lr = 0.0005;
steps = 1;

% simulated data
    environment = [zeros(no_data_points/2,T); 0.8*ones(no_data_points/2, T)];
    environment = abs(environment + 0.05*randn(no_data_points, T));
    environment = environment ./ repmat(sum(environment), [no_data_points 1]);

O = organism;
org = initialise_organism(O, eps, beta, nx, mx);


data = environment(:,1);


TT = 20000;
E = zeros(TT,1);


for t = 1 : TT
    
    t
    
    org.stat_dist = compute_stationary_distribution(org, eps, size(data,1));
    [org] = optimise_HMM_weights(org, data, steps, lr, eps, beta);
    
    
    norm(org.stat_dist - data)
    E(t) = norm(org.stat_dist - data);
    
    if mod(t, 10) == 0
        org.transition_matrix
    end

end
    