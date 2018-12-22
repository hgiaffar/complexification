%% Run basic first simulation

T = 50000;           % length of simulation
no_orgs = 100;       % no organisms
Kc = 5;              % number of k-cycle selections per timestep
no_data_points = 10;

% simulated data
% must be normalised
if 1
    DATA = [zeros(no_data_points/2,T); 0.8*ones(no_data_points/2, T)];
%     environment = abs(DATA + 0.05*randn(no_data_points, T));
    environment = environment ./ repmat(sum(environment), [no_data_points 1]);
else
    environment = ones(no_orgs,T);
end

nx = 5;            % visible units x
mx = 2;
beta = 1;
eps = 1e-8;
eta = 1;

O = organism;
mf = [];

population = struct();

% initialise population
for i = 1 : no_orgs
    population(i).organism = initialise_organism(O, eps, beta, nx, mx);
    population(i).lineage = i;
end

figure, 

for t = 1 : T
        
    % interaction with environment - compute fitness
    data = environment(:,t);
    
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, data, eps);
    end
        
    [population] = evo_dynamics_step(population, 'Moran', Kc, eta, eps, beta);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps')
    drawnow,
        
    if mod(t,100) == 0
        
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        %         figure, subplot(1,2,1), imagesc(orr), colorbar, subplot(1,2,2), imagesc(reshape(data, [2 2])), colorbar
        
        t
        orr.stat_dist
        data
        
%         mf(end)
%         A = population(1).organism;
%         A.weights
%         B = population(2).organism;
%         B.weights

    end
    
end