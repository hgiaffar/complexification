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
    end
    
    methods
        
        % initialise weights and transition matrix
        
        function [obj] = initialise_organism(obj, eps, beta, nx, mx)
            % obj - organism
            % eps - epsilon for weight to transition update
            % beta - temperature
            % nx, mx = size 1st, 2nd dim of environmental data
            
            states = nx * mx;
            T = rand(states);
            V = zeros(size(T));

            for j = 1 : states
                su = sum(exp(beta * T(:,j)));
                for i = 1 : states
                    V(i,j) = ((1-eps) * exp(beta * T(i,j)) / su ) + eps; 
                end 
            end
            
            obj.weights = T;
            obj.transition_matrix = V;
        end
        
        
        % first variational operator - no change to topology, mutate weights
        function [obj] = var_operator_1(obj, eta, eps, beta)
            % eta - mutation scaling parameter
            
            wmax = 1;
            
            indx = randi(size(obj.weights,1), [2 1]);         % select indices to mutate
            obj.weights(indx(1), indx(2)) = obj.weights(indx(1), indx(2)) + (randn * eta);

            obj.weights(indx(1), indx(2)) = abs(min(obj.weights(indx(1), indx(2)), wmax));               % bound weights
%                 obj.weights(:, indx(2)) = abs(min(obj.weights(indx(1), indx(2)), wmax));               % bound weights
            
            obj = update_transitions(obj, indx, eps, beta);
        end
        
        
        % second variational operator - change in hidden unit topology by copying internal structure
        function [obj] = var_operator_2(obj, form, numn, eps, beta)
            % form: structured - whole subnetwork copied, inc. weights. unstructured - weights sampled from
            % network distribution
            % numn - number of nodes sampled
            
            num_nodes = size(obj.weights,1);
            wmax = 1;
            
            switch form
                case 'structured'
                    ind_copy = randperm(num_nodes, 1);
                    ind2 = find(obj.weights(ind_copy,:));
                    ind_copy = [ind_copy ind2(randperm(length(ind2), numn-1))];
                    
                    % copy all weights                                     
                    ind3 = [];
                    for i = 1: numn, ind3 = [ind3; randperm(num_nodes, num_nodes)]; end
                    
                    E = [obj.weights; obj.weights(ind_copy,ind3)];
                    F = [obj.weights(ind3, ind_copy); obj.weights(ind_copy,ind3) * obj.weights(ind3, ind_copy) ]; 
                    
                    obj.weights = [E F];
                    
                    % normalise weights
%                     obj.weights = obj.weights ./ repmat(vecnorm(obj.weights), [size(obj.weights,1) 1]);
                    
                    obj.weights = abs(min(obj.weights, wmax));
                    
                    obj = update_transitions(obj, 0, eps, beta);
                    
                    if sum(isnan(obj.weights(:))) > 0
                        disp('fucked weights')
                    elseif sum(isnan(obj.transition_matrix(:))) > 0
                        disp('fucked transitions')
                    end
                    
                case 'unstructured'
            
                    ind_copy = randperm(num_nodes, numn);
                    
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
        function [stat_dist] = compute_stationary_distribution(obj, eps, num_vis_units)
           
            [Evec, Eval] = eigs(obj.transition_matrix);
            
            if 1-eps <= Eval(1) <= 1+eps
               Eval(1);
               stat_dist = real(Evec(:,1)); 
               stat_dist = stat_dist(1:num_vis_units)/sum(stat_dist(1:num_vis_units));
            else
                error('stationary distribution issue')
            end
        end
           
        
        % compute fitness as KLdiv btw organism and next environmental input 
        function [obj, f] = compute_fitness(obj, funct, data, eps)
            
            num_vis_units = max(size(data));
        
            [obj.stat_dist] = compute_stationary_distribution(obj, eps, num_vis_units);            % phenotype
            
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
            