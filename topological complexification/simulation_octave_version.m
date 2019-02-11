%% Octave friendly version

clear all

%% 0 - parameters

T = 30000;               % length of simulation
et = 5000;               % number of steps during which enviroment is constant
no_orgs = 100;          % number of organisms
Kc = 1;                 % number of k-cycle selections per timestep
no_data_points = 16;    % size of the data

varop = 4;            % variational operator
P = [0.7 0.2 0.1];          % probability distribution over variational operators
copdelP = [0.5 0.5];  % distribution over copy vs delete.
numn = 1;             % number of edges or nodes or size of substructure copied/deleted
no_hidden_units = 0;
graph_density = 0.5;
maxK = 5;               % maximum length path considered in the calculation of phenotype (A^k)

full_size = 200;            % PREALLOCATE MEMORY

% generative model of the environment

envt_type = 'complex';
no_samples = ceil(T/et);

tree_type = 'binary';           % type of tree structure for the 
sigma_envt = 0.2;               % variance of the gaussian centered at 0 for each edge in the graph
beta_envt = 1;                  % softmax parameter

%% 1 - generate environment

switch envt_type
    case 'simple'
        level_no = 1;
        branch_no = no_data_points;
    case 'intermediate'
        level_no = 2;
        branch_no = 4;
    case 'complex'
        level_no = 4;
        branch_no = 2;
end

%debug_on_error(1)
environment = generate_nonbinary_hierarchical_envt(no_data_points, branch_no, level_no, no_samples, sigma_envt, beta_envt, tree_type);


%% 2 - initialise population
O = organism;
mf = [];

population = struct();
save_network = struct();            % save intermediate fittest weights for figure

for i = 1 : no_orgs
    population(i).organism = initialise_organism(O, no_data_points, graph_density, no_hidden_units, full_size);
    population(i).lineage = i;
end

hh =  findobj('type','figure');
nn = length(hh);
figure(nn+1), 
c = 1;
cc = 0;
ccc = 1;


%% 3 - simulation

for t = 1 : T
    
    % sample environment
    if t == 1
        data = environment(ccc,:)'
    elseif mod(t-1,et) == 0
        ccc = ccc+1;
        data = environment(ccc,:)'
    end
    
    % update population
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, 'l1norm', data );
        population(j).num_hidden_nodes = population(j).organism.num_hidden_nodes;
        population(j).genotype_complexity = population(j).organism.genotype_complexity;
    end
        
    [population] = evo_dynamics_step(population, 'Moran', Kc, varop, numn, copdelP, maxK, P);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    figure(nn+1), plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps'), drawnow
    
%     mean([population(1:length(population)).num_hidden_nodes])
        

    %% PLOT FITNESS and GENOTYPIC/PHENOTYPIC COMPLEXITY
    
    % update some other parameters
%     if mod(t, 50) == 0
%         cc = cc+1;
%         mean_hidden_nodes(cc) = mean([population(1:no_orgs).num_hidden_nodes]);
%         std_hidden_nodes(cc) = std([population(1:no_orgs).num_hidden_nodes]);
%         ind2 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]));
%         ind3 = find([population(1:no_orgs).fitness] == min([population(1:no_orgs).fitness]));
%         fittest_hidden_fitness(cc) = population(ind2(1)).fitness;
%         mean_fitness(cc) = mean([population(1:no_orgs).fitness]);
%         weakest_hidden_fitness(cc) = population(ind3(1)).fitness;
%         fittest_hidden_nodes(cc) = population(ind2(1)).num_hidden_nodes;
%         weakest_hidden_nodes(cc) = population(ind3(1)).num_hidden_nodes; 
%        
%         PC(cc) = sample_phenotypic_complexity( population, 50000, 0.1 );
%         GC(cc) = mean([population(1:no_orgs).genotype_complexity]);
%     end
    
    % display
    if mod(t,200) == 0 
        c = c+1;
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        t, orr.phenotype, data
        mean([population(1:no_orgs).num_hidden_nodes])
        std([population(1:no_orgs).num_hidden_nodes])    
    end

    % save graph corresponding to the current fittest organism
    if t == et - 1
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        save_network.first_network = orr.adjacency;

    elseif t == T/2 - 1

        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        save_network.middle_network = orr.adjacency;

    elseif t == T - 1

        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        save_network.final_network = orr.adjacency;

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


