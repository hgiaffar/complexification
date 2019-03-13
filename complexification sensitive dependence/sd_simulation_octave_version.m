%% Octave friendly version

clear all

%% 0 - parameters

% simulation
T = 40000;              % length of simulation
et = 400;              % number of steps during which enviroment is constant
no_orgs = 100;          % number of organisms
Kc = 1;                 % number of k-cycle selections per timestep
no_data_points = 16;    % size of the data
no_hidden_units = 0;    % initial number of hidden nodes
graph_density = 0.5;    % initial graph density
maxK = 5;               % maximum length path considered in the calculation of phenotype (~ A^k)
slct_strength = 100;    % fitness distribution can be sharpened (have to change in evo_dynamics_step)
preall_size = 200;      % PREALLOCATE MEMORY

% mutation
varop = 4;              % variational operator
P = [0.6 0.2 0.2];      % probability distribution over variational operators
copdelP = [0.5 0.5];    % distribution over copy vs delete.
numn = 1;               % number of edges or nodes or size of substructure copied/deleted

% sampling phenotypic complexity
calc_pc = 0;            % compute phenotypic complexity as part of the Moran cycle step
num_samples_pc = 10;    % how many anchors in the phenotype simplex used to approximate the distribution          
num_mutations_pc = 10;  % number of phenotypes sampled through random mutation for a given genotype
sigma_pc = 0.1;         % variance of gaussian centered at each anchor point

% generative model of the environment
envt_type = 'complex';          % predefined variables for environment generation
no_samples = ceil(T/et);        % number of environments generated from the hierarchical generative model
tree_type = 'binary';           % type of tree structure  
sigma_envt = 0.2;               % variance of the gaussian centered at 0 for each edge in the graph
beta_envt = 2;                  % softmax parameter over generated environment

% report during simulation
reporting_interval = 100;       % for updating plots and calculating phenotypic complexity
check_adjacency = 100;          % check that the number of nodes isn't near the preallocated space for adjacency matrices
add_preacllocated_nodes = 50;   % if it is, by threshold measure below, add this many nodes
num_nodes_thresh = 0.75;        




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
mnhn = [];
maxnhn = [];
mgc= [];
maxgc = [];
population = struct();
save_network = struct();            % save intermediate fittest weights for figure

for i = 1 : no_orgs
    population(i).organism = initialise_organism(O, no_data_points, graph_density, no_hidden_units, preall_size);
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
        data = environment(ccc,:)';
    elseif mod(t-1,et) == 0
        ccc = ccc+1;
        data = environment(ccc,:)';
    end
    % update population
    for j = 1 : no_orgs
        [population(j).organism, population(j).fitness] = compute_fitness(population(j).organism, 'l1norm', data );
        population(j).num_hidden_nodes = population(j).organism.num_hidden_nodes;
        population(j).genotype_complexity = population(j).organism.genotype_complexity;
        population(j).number_nodes = population(j).organism.number_nodes;
    end
        
    [population] = evo_dynamics_step(population, 'Moran', Kc, varop, numn, copdelP, maxK, P, slct_strength, calc_pc);
    
    mf = [mf mean([population(1:no_orgs).fitness])];
    mnhn = [mnhn mean([population(1:no_orgs).num_hidden_nodes])];
    maxnhn = [maxnhn max([population(1:no_orgs).num_hidden_nodes])];
    mgc = [mgc mean([population(1:no_orgs).genotype_complexity])];
    maxgc = [maxgc max([population(1:no_orgs).genotype_complexity])];
    
    % display
    if mod(t,reporting_interval) == 0 
        c = c+1;
        ind1 = find([population(1:no_orgs).fitness] == max([population(1:no_orgs).fitness]))
        orr = population(ind1).organism;
        t, [data orr.phenotype] 
        
        figure(nn+1), 
            subplot(2,2,1), plot(1:length(mf), mf), ylabel('mean fitness'), xlabel('steps')
            subplot(2,2,2), plot(1:length(mnhn), mnhn), hold on, plot(maxnhn), ylabel('no. hidden nodes'), xlabel('steps'), legend('mean', 'max')
            subplot(2,2,3), plot(1:length(mgc), mgc), hold on, plot(maxgc), ylabel('genotypic complexity'), xlabel('steps'), legend('mean', 'max')
            
            drawnow
            
         % compute phenotypic complexity for the population every reporting interval
         if 0
         for i1 = 1 : no_orgs
            population(i1).phenotype_complexity(c) = ...
                sample_phenotypic_complexity( population(i1).organism, num_samples_pc, num_mutations_pc, sigma_pc, P, copdelP, numn, maxK );
         end, end
         
%         subplot(2,2,4), plot(1:c,mean([population(i1).phenotype_complexity(c)]))
    
    end
    
    % Are organisms near the limit of the preallocated memory for adjacency?  If so, increase the size fe agent.
    if mod(t,check_adjacency) == 0
        if max([population(1:no_orgs).number_nodes]) > num_nodes_thresh * preall_size    
            % add add_preacllocated_nodes
            for p1 = 1 : no_orgs
                adj = population(p1).organism.adjacency;
                population(p1).organism.adjacency = [[adj zeros(size(adj,1), add_preacllocated_nodes)]; ...
                    zeros(add_preacllocated_nodes, size(adj,2)+add_preacllocated_nodes)];
            end
            preall_size = preall_size + add_preallocated_nodes;
        end
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


% figure(nn+2),
% subplot(1,2,1), errorbar(1:length(mean_hidden_nodes), mean_hidden_nodes, std_hidden_nodes), hold on, plot(fittest_hidden_nodes), hold on, plot(weakest_hidden_nodes) 
%         legend('mean', 'fittest', 'weakest'), xlabel('times'), ylabel('number of hidden nodes')
%         
% subplot(1,2,2), plot(mean_fitness), hold on, plot(fittest_hidden_fitness), hold on, plot(weakest_hidden_fitness), legend('mean', 'fittest', 'weakest')
%         xlabel('times'), ylabel('mean fitness')
% 
% figure(nn+3),
% 
%         yyaxis left, plot(PC), ylabel('phenotypic complexity')
%         yyaxis right, plot(GC), ylabel('genotypic complexity')


G1 = digraph(save_network.first_network > 0.5);
G2 = digraph(save_network.middle_network > 0.5);
G3 = digraph(save_network.final_network > 0.5);


