classdef organism
    
    properties
        % weight matrix
        weights
        % transition matrix
        transition_matrix
        % static distribution
        stat_dist
        % number of hidden nodes
        num_hidden_nodes
        % adjacency matrix
        adjacency
        % full stat_dist - for optimisation via gd
        stat_dist_full
        % genotypic complexity
        genotype_complexity
    end
    
    methods
        
        % initialise weights and transition matrix
        
        function [obj] = initialise_organism(obj, eps, beta, nx, mx, rho, nhu)
            % obj - organism
            % eps - epsilon for weight to transition update
            % beta - temperature
            % nx, mx = size 1st, 2nd dim of environmental data
            % rho - hard density of connections in T/V - defines adjacency mat
            % nhu - number of hidden units upon initialisation
            
            states = (nx * mx) + nhu;
            T = rand(states);
            A = random_fsparse_mat(states, rho);
            T = T.*A;
            V = zeros(size(T));

            for j = 1 : states
                su = sum(exp(beta * T(:,j)));
                for i = 1 : states
                    if T(i,j) ~= 0
                        V(i,j) = ((1-eps) * exp(beta * T(i,j)) / su ) + eps; 
                    end
                end
                V(:,j) = V(:,j)/sum(V(:,j));
            end
            
            
            obj.weights = T;
            obj.transition_matrix = V;
            obj.adjacency = A;

            obj.genotype_complexity = length(find(obj.adjacency(:)));
            
        end
        
        
        % first variational operator - no change to topology, mutate weights
        function [obj] = var_operator_1(obj, eta, eps, beta)
            % eta - mutation scaling parameter
            
            wmax = 1;
            
            indx = randi(size(obj.weights,1), [2 1]);         % select indices to mutate
            obj.weights(indx(1), indx(2)) = obj.weights(indx(1), indx(2)) + (randn * eta);

            obj.weights(indx(1), indx(2)) = abs(min(obj.weights(indx(1), indx(2)), wmax));               % bound weights
            
            obj = update_transitions(obj, indx, eps, beta);

            obj.genotype_complexity = length(find(obj.adjacency(:)));

        end

        
        function [obj] = var_operator_2(obj, data, steps, eps, lr, beta, copdel)
            % add/delete edge
            % copdel - add or remove
            
            [ri, ci] = find(obj.adjacency);
            [rii, cii] = find(~obj.adjacency);
            
            switch copdel
                case 1
                    indy = randperm(length(ri), 1);
                    indx = randperm(length(rii), 1);
                    obj.weights(rii(indx), cii(indx)) = obj.weights(ri(indx), ci(indx));
                    obj.adjacency(rii(indx), cii(indx)) = 1;
                    obj = update_transitions(obj, 0, eps, beta);
                    
                case -1
                    indx = randperm(length(ri), 1);
                    obj.weights(ri(indx), ci(indx)) = 0;
                    obj.adjacency(ri(indx), ci(indx)) = 0;
                    obj = update_transitions(obj, 0, eps, beta);
            end                    

            obj.genotype_complexity = length(find(obj.adjacency(:)));

        end  
        
        
        
        % third variational operator - change in hidden unit topology by copying internal structure
        function [obj] = var_operator_3(obj, form, numn, eps, beta, copdel)
            % form: structured - whole subnetwork copied, inc. weights. unstructured - weights sampled from
            % network distribution
            % numn - size of substructure copied/deleted
            % copdel - copied/deleted
            
            num_nodes = size(obj.weights,1);
            wmax = 1;
            
            switch copdel
                
                case 1                  % add substructure
                    switch form
                        case 'structured'
                            % chose random index
                            ind_copy = randperm(num_nodes, numn);

                            while ~isempty(ind_copy)
                                % copy weights and adjacency as are, randomly generating the missing values
                                cent = obj.adjacency(ind_copy(1),ind_copy(1)); 
                                lef = [obj.adjacency(ind_copy(1),:) cent];
                                rig = obj.adjacency(:,ind_copy(1));
                                density = sum(obj.adjacency) / numel(obj.adjacency);
                                missing = rand(2,1) > (1-density);

                                obj.adjacency = [[obj.adjacency rig]; lef];
                                obj.adjacency(end, ind_copy(1)) = missing(1); 
                                obj.adjacency(ind_copy(1), end) = missing(2);

                                cent = obj.weights(ind_copy(1),ind_copy(1)); 
                                lef = [obj.weights(ind_copy(1),:) cent];
                                rig = obj.weights(:,ind_copy(1));

                                obj.weights = [[obj.weights rig]; lef];
                                obj.weights(end, ind_copy(1)) = rand(1); 
                                obj.weights(ind_copy(1), end) = rand(1);                                     

                                obj.weights = abs(min(obj.weights, wmax));
                                obj.weights = obj.weights .* obj.adjacency;

                                obj = update_transitions(obj, 0, eps, beta);

                                if sum(isnan(obj.weights(:))) > 0
                                    disp('fucked weights')
                                elseif sum(isnan(obj.transition_matrix(:))) > 0
                                    disp('fucked transitions')
                                end

                               ind_copy(1) = []; 
                            end

                        case 'unstructured'

                            ind_copy = randperm(num_nodes, numn);
                            density = sum(obj.adjacency(:)) / numel(obj.adjacency);
                            no_unts = size(obj.adjacency,1);

                            lef = [rand([(no_unts + 1) 1]) > (1-density)]';
                            rig = [rand(no_unts,1) > (1-density)];

                            obj.adjacency = [[obj.adjacency rig]; lef];

                            obj.weights = [[obj.weights rand(no_unts, 1)]; rand(1,no_unts+1)];
                            obj.weights = obj.weights .* obj.adjacency;

                            obj = update_transitions(obj, 0, eps, beta);
                           
                    end
                    
                    % change the number of hidden nodes
                    obj.num_hidden_nodes = obj.num_hidden_nodes + numn;
                    
                case -1             % delete substructure
                    
                    if numn <= obj.num_hidden_nodes
                        ind_copy = randperm(obj.num_hidden_nodes, numn) + length(obj.stat_dist);
                        obj.weights(ind_copy, :) = [];
                        obj.weights(:, ind_copy) = [];
                        obj.adjacency(ind_copy, :) = [];                        
                        obj.adjacency(:, ind_copy) = [];                        
                        obj = update_transitions(obj, 0, eps, beta);
                        obj.num_hidden_nodes = obj.num_hidden_nodes - numn;
                    else
                        numn = obj.num_hidden_nodes;
                        ind_copy = randperm(obj.num_hidden_nodes, numn) + length(obj.stat_dist);
                        obj.weights(ind_copy, :) = [];
                        obj.weights(:, ind_copy) = [];
                        obj.adjacency(ind_copy, :) = [];                                                
                        obj.adjacency(:, ind_copy) = [];                                                
                        obj = update_transitions(obj, 0, eps, beta);
                        obj.num_hidden_nodes = obj.num_hidden_nodes - numn;                  
                    end
                end

            obj.genotype_complexity = length(find(obj.adjacency(:)));

        end
        
        
        
        function [obj] = var_operator_distributed(obj, form, numn, eps, eta, beta, P, copdelP)
            % select above variational operators based on distribution P
            % form: structured vs unstructured
            % P - probability distribution over other operators in order
            % Pcopdel is the Pdist over delete/copy 
            
            indy = RouletteWheelSelection(copdelP);
            if indy == 1, copdel = -1; elseif indy == 2, copdel = 1; end 
            
            indx = RouletteWheelSelection(P);
                        
            switch indx
                case 1
                     obj = var_operator_1(obj, eta, eps, beta);
                     
                case 2
                     obj = var_operator_2(obj, data, steps, eps, lr, beta, copdel);
                     
                case 3
                    obj = var_operator_3(obj, form, numn, eps, beta, copdel);
            end       
        end
        
        
        
        % update transition matrix given points of mutation
        function [obj] = update_transitions(obj, indx, eps, beta)
            % indx - index of mutation - recompute only relevant column
            % indx == 0 - recompute whole transition matrix as mutation is topological
            
            if indx > 0
                 suma = sum(exp(beta * obj.weights(:,indx(2))));  
                 for i = 1 : size(obj.weights,1)
                    obj.transition_matrix(i,indx(2)) = ((1-eps) * exp(beta * obj.weights(i, indx(2))) / suma ) + eps; 
                 end
            elseif indx == 0
                V = zeros(size(obj.weights));
                for j = 1 : size(obj.weights,2)
                    su = sum(exp(beta * obj.weights(:,j)));
                    for i = 1 : size(obj.weights,1)
                        V(i,j) = ((1-eps) * exp(beta * obj.weights(i,j)) / su ) + eps; 
                    end 
                end
                obj.transition_matrix = V;
            end 
        end
        
        
        % compute stationary distribution of HMM
        function [stat_dist, stat_dist_full] = compute_stationary_distribution(obj, eps, num_vis_units)
           
            [Evec, Eval] = eigs(obj.transition_matrix);
            
            if 1-eps <= Eval(1) <= 1+eps
               Eval(1);
               stat_dist = real(Evec(:,1)); 
               stat_dist = stat_dist(1:num_vis_units)/sum(stat_dist(1:num_vis_units));
               stat_dist_full = real(Evec(:,1));
            else
                error('stationary distribution issue')
            end
        end
           
        
        % compute fitness
        function [obj, f] = compute_fitness(obj, funct, data, eps)
            
            num_vis_units = max(size(data));
        
            [obj.stat_dist, obj.stat_dist_full] = compute_stationary_distribution(obj, eps, num_vis_units);            % phenotype
            
            switch funct
                case 'kldiv'
                    f = exp(-real(sum((data+1e-8) .* log2((data+1e-8)./(obj.stat_dist+1e-8)))));

                case 'l1norm'
                    % L1 norm 
                    f = mean(abs(obj.stat_dist - data));
            end
            
            obj.num_hidden_nodes = size(obj.weights,1) - max(size(data));
            
            if f == Inf
                disp('fitness: inf')
            end
        end
            
    end
end
            