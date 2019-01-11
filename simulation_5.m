%% 

clear all

T = 20000;           % length of simulation
et = 300;           % number of steps during which enviroment is constant
no_orgs = 100;       % no organisms
Kc = 5;              % number of k-cycle selections per timestep
no_data_points = 3;

% simulated data must be normalised
[environment] = generate_environment(no_data_points, 'simple', T);

nx = 3;            % visible units 
mx = 1;
beta = 10;
eps = 1e-8;
eta = 0.005;

varop = 4;            % variational operator
P = [0.9 0 0.1];          % probability distribution over variational operators
copdelP = [0.5 0.5];  % distribution over copy vs delete.
numn = 1;               % size of substructure copied/deleted
form = 'structured';  % is the copyong structured or unstructured?

O = organism;
mf = [];
no_hidden_units = 1;
graph_density = 0.8;

steps = 1;
lr = 0.00005;

population = struct();
save_network = struct();            % save intermediate fittest weights for figure

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
ccc = 1;

for t = 1 : T
    
    % generate enviroments
    if t == 1
        data = environment(:,ccc)
    elseif mod(t,et) == 0
        ccc = ccc+1;
        data = environment(:,ccc)
    end
    
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, 'kldiv', data, eps);
        population(j).num_hidden_nodes = population(j).organism.num_hidden_nodes;
        population(j).genotype_complexity = population(j).organism.genotype_complexity;
    end
        
    [population] = evo_dynamics_step(population, 'Moran', Kc, eta, eps, beta, varop, form, P, copdelP, numn);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    figure(nn+1), plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps'), drawnow
        
    
    if mod(t, 50) == 0
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
       
        PC(cc) = sample_phenotypic_complexity( population, 50000, 0.1 );
        GC(cc) = mean([population(1:no_orgs).genotype_complexity]);
    end
    
    % display
    if mod(t,200) == 0 
        
        c = c+1;
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        t, orr.stat_dist, data
        mean([population(1:no_orgs).num_hidden_nodes])
        std([population(1:no_orgs).num_hidden_nodes])   
        
    end


    if t == et - 1

        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
    
        save_network.first_network = orr.weights;

    elseif t == T/2 - 1

        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        
        save_network.middle_network = orr.weights;

    elseif t == T - 1

        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;

        save_network.final_network = orr.weights;

    end

end


figure(nn+2),
subplot(1,2,1), errorbar(1:length(mean_hidden_nodes), mean_hidden_nodes, std_hidden_nodes), hold on, plot(fittest_hidden_nodes), hold on, plot(weakest_hidden_nodes) 
        legend('mean', 'fittest', 'weakest'), xlabel('times'), ylabel('number of hidden nodes')
        
subplot(1,2,2), plot(mean_fitness), hold on, plot(fittest_hidden_fitness), hold on, plot(weakest_hidden_fitness), legend('mean', 'fittest', 'weakest')
        xlabel('times'), ylabel('mean fitness')

figure(nn+3),

        yyaxis left, plot(PC), ylabel('phenotypic complexity')
        yyaxis right, plot(GC), ylabel('genotypic complexity')


G1 = digraph(save_network.first_network > 0.5);
G2 = digraph(save_network.middle_network > 0.5);
G3 = digraph(save_network.final_network > 0.5);


