%% Run basic first simulation

T = 5000;           % length of simulation
no_orgs = 100;       % no organisms
Kc = 5;              % number of k-cycle selections per timestep
no_data_points = 4;

% simulated data
% must be normalised
if 1
    environment = [zeros(no_data_points/2,T); 0.8*ones(no_data_points/2, T)];
    environment = abs(environment + 0.05*randn(no_data_points, T));
    environment = environment ./ repmat(sum(environment), [no_data_points 1]);
else
    environment = ones(no_data_points,T);
    environment = environment ./ repmat(sum(environment), [no_data_points 1]);

end

nx = 2;            % visible units x
mx = 2;
beta = 100;
eps = 1e-8;
eta = 0.5;

O = organism;
mf = [];

population = struct();

% initialise population
for i = 1 : no_orgs
    population(i).organism = initialise_organism(O, eps, beta, nx, mx);
    population(i).lineage = i;
end

figure(1), 
c = 1;

for t = 1 : T
        
    % interaction with environment - compute fitness
    data = environment(:,t);
    
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, 'kldiv', data, eps);
        population(j).num_hidden_nodes = population(j).organism.num_hidden_nodes;
    end
        
    [population] = evo_dynamics_step(population, 'Moran', Kc, eta, eps, beta);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    figure(1), plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps')
    drawnow,
        
    
    % display
    if mod(t,200) == 0
        
        c = c+1;
        
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        figure(c), subplot(1,3,1), imagesc(reshape(orr.stat_dist, [nx mx])), colorbar, subplot(1,3,2), imagesc(reshape(data, [nx mx])), colorbar
                   subplot(1,3,3), imagesc(abs(reshape(orr.stat_dist, [nx mx]) - reshape(data, [nx mx]))), colorbar
        
        t, orr.stat_dist, data
        
%         mf(end)
%         A = population(1).organism;
%         A.weights
%         B = population(2).organism;
%         B.weights
        
        mean([population(1:no_orgs).num_hidden_nodes])
        std([population(1:no_orgs).num_hidden_nodes])
        
    end
    
end