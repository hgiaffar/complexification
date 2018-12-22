%% Must fix issue - if connection probability is too low, end up with disconnected nodes
% one way to deal would be to add a self loop? or require a single connection to define a vertex (organism must be connected)

clear all

T = 2500;           % length of simulation
no_orgs = 100;       % no organisms
Kc = 1;              % number of k-cycle selections per timestep
no_data_points = 4;

% simulated data must be normalised
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
beta = 10;
eps = 1e-8;
eta = 0.005;

varop = 4;            % variational operator
P = [0 0 1];          % probability distribution over variational operators

O = organism;
mf = [];
no_hidden_units = 5;
graph_density = 1;

form = 'unstructured';
steps = 1,
lr = 0.00005;

population = struct();


% initialise population
for i = 1 : no_orgs
    population(i).organism = initialise_organism(O, eps, beta, nx, mx, graph_density, no_hidden_units);
    population(i).lineage = i;
end

hh =  findobj('type','figure');
nn = length(hh);
figure(nn+1), 
c = 1;
cc = 0;

for t = 1 : T
        
    % interaction with environment - compute fitness
    data = environment(:,t);
    
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, 'kldiv', data, eps);
        population(j).num_hidden_nodes = population(j).organism.num_hidden_nodes;
    end
        
%     [population] = evo_dynamics_step(population, 'Moran', Kc, eta, eps, beta, varop, 'structured', P);
    [population] = evo_dynamics_step(population, 'Moran', Kc, eta, eps, beta, varop, form, P, data, steps, lr);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    figure(nn+1), plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps')
    drawnow,
        
    
    if mod(t, 20) == 0
        cc = cc+1;
        mean_hidden_nodes(cc) = mean([population(1:no_orgs).num_hidden_nodes]);
        std_hidden_nodes(cc) = std([population(1:no_orgs).num_hidden_nodes]);
        ind2 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]));
        ind3 = find([population(1:no_orgs).fitness] == min([population(1:no_orgs).fitness]));
        fittest_hidden_fitness(cc) = population(ind2(1)).fitness;
        mean_fitness(cc) = mean([population(1:no_orgs).fitness]);
        weakest_hidden_fitness(cc) = population(ind3(1)).fitness;
        fittest_hidden_nodes(cc) = population(ind2(1)).num_hidden_nodes;
        weakest_hidden_nodes(cc) = population(ind3(1)).num_hidden_nodes;        
    end
    
    % display
    if mod(t,200) == 0    
        c = c+1;
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
%         figure(c), subplot(1,3,1), imagesc(reshape(orr.stat_dist, [nx mx])), colorbar, subplot(1,3,2), imagesc(reshape(data, [nx mx])), colorbar
%                    subplot(1,3,3), imagesc(abs(reshape(orr.stat_dist, [nx mx]) - reshape(data, [nx mx]))), colorbar
        
        t, orr.stat_dist, data
        mean([population(1:no_orgs).num_hidden_nodes])
        std([population(1:no_orgs).num_hidden_nodes])      
        
    end
    
end


figure,
subplot(1,2,1), errorbar(1:length(mean_hidden_nodes), mean_hidden_nodes, std_hidden_nodes), hold on, plot(fittest_hidden_nodes), hold on, plot(weakest_hidden_nodes) 
        legend('mean', 'fittest', 'weakest'), xlabel('times'), ylabel('number of hidden nodes')
        
subplot(1,2,2), plot(mean_fitness), hold on, plot(fittest_hidden_fitness), hold on, plot(weakest_hidden_fitness), legend('mean', 'fittest', 'weakest')
        xlabel('times'), ylabel('mean fitness')
