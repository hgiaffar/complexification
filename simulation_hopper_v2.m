%% Octave friendly version

clear all

%% 0 - parameters

parameters_hopper_v2


%% 1 - generate environment

switch envt_type
    case 'simple'
        level_no = 1;
        branch_no = no_data_points;
        environment = abs(randn(no_data_points, no_samples));
        environment = (environment ./ repmat(sum(environment), [no_data_points 1]))';
    case 'complex'
        level_no = 4;
        branch_no = 2;
        environment = generate_nonbinary_hierarchical_envt(no_data_points, branch_no, level_no, no_samples, sigma_envt, beta_envt, tree_type);

end

%debug_on_error(1)

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
    population(i).phenotype_complexity = [];
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
         if 1
         for i1 = 1 : no_orgs
                sample_phenotypic_complexity( population(i1).organism, num_samples_pc, num_mutations_pc, sigma_pc, P, copdelP, numn, maxK, data );
                population(i1).phenotype_complexity = [population(i1).phenotype_complexity population(i1).organism.phenotype_complexity];
         end, end
         
%         subplot(2,2,4), plot(1:c,mean([population(i1).phenotype_complexity(c)]))
    
    end
    
    if mod(t, clear_leak_interval) == 0
        disp('clearing variables')
        save('current_sim_data', 'population', 'mf', 'mnhn', 'maxnhn', 'mgc', 'maxgc', 'data', 'environment', 'c', 'cc', 'ccc', 't', 'nn', '-v7.3')
        clear all
        parameters_hopper_v1
        load current_sim_data.mat
    end
    
    
    
    % Are organisms near the limit of the preallocated memory for adjacency?  If so, increase the size fe agent.
    if mod(t,check_adjacency) == 0
        if max([population(1:no_orgs).number_nodes]) > num_nodes_thresh * preall_size    
            % add add_preacllocated_nodes
            for p1 = 1 : no_orgs
                adj = population(p1).organism.adjacency;
                population(p1).organism.adjacency = [[adj zeros(size(adj,1), add_preallocated_nodes)]; ...
                    zeros(add_preallocated_nodes, size(adj,2)+add_preallocated_nodes)];
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


% G1 = digraph(save_network.first_network > 0.5);
% G2 = digraph(save_network.middle_network > 0.5);
% G3 = digraph(save_network.final_network > 0.5);


