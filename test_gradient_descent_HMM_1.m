%% test gradient descent in HMM model

% include adaptive learning rate

% seems that this works for MMs no hidden units, but cannot find the
% correct approach to mutate weights directly (mutations should apply to the weights, not transition probabilities)


clear all

T = 100;
no_data_points = 4;

eps = 1e-8;
beta = 100;
nx = 2;
mx = 2;

graph_density = 1;
no_hidden_units = 0;

lr = 0.0005;
steps = 1;

% simulated data
    environment = [zeros(no_data_points/2,T); 0.8*ones(no_data_points/2, T)];
    environment = abs(environment + 0.05*randn(no_data_points, T));
    environment = environment ./ repmat(sum(environment), [no_data_points 1]);

O = organism;
org = initialise_organism(O, eps, beta, nx, mx, graph_density, no_hidden_units);

data = environment(:,1);


TT = 20000;
E = zeros(TT,1);

figure,


for t = 1 : TT
    
    t
    
    [org.stat_dist, org.stat_dist_full] = compute_stationary_distribution(org, eps, size(data,1));
    [org] = optimise_HMM_weights(org, data, steps, lr, eps, beta);
    
    E(t) = norm(org.stat_dist - data);
    
    if mod(t, 100) == 0
        org.transition_matrix;
    end
    
    if mod(t,100) == 0
        plot(E)
        drawnow
    end

end
    