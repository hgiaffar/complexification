classdef organism
    
    % topological variation only
    
    properties
        % number of visible nodes
        num_visible_nodes
        % number of hidden nodes
        num_hidden_nodes
        % adjacency matrix - 'genotype'
        adjacency
        % phenotype
        phenotype
        % genotypic complexity
        genotype_complexity
        % sample phenotypic complexity
        phenotype_complexity
        % 1-D hole spectrum
        num_cycles
        % number of nodes in network
        number_nodes
    end
    
    methods
        
        %% initialise organism        
        function [obj] = initialise_organism(obj, data_size, rho, nhu, full_size)
            % obj - organism
            % data_size - fixed size of data
            % rho - hard density of connections in T/V - defines adjacency mat
            % nhu - number of hidden units upon initialisation
            % full_size - preallocate memory
            
            obj.number_nodes = data_size + nhu;
            A = random_fsparse_mat(obj.number_nodes, rho);
            obj.adjacency = [[A zeros(obj.number_nodes, full_size-data_size)]; [zeros(full_size-data_size) zeros(full_size-data_size, data_size)]];
            obj.genotype_complexity = length(find(obj.adjacency(:)));
            obj.num_visible_nodes = data_size;
            obj.phenotype = zeros(data_size,1);
            
        end
        
        
        %% VAR I - add/delete edges
        function [obj] = var_operator_1(obj, num_edge, copdel)
            % num_edge - number of edges added/deleted
            % copdel - copy or delete
            
            % functional adjacency
            funct_adj = obj.adjacency(1:obj.number_nodes, 1:obj.number_nodes);
            
            switch copdel
                case 1

                    [indx, indy] = find(~funct_adj);
                    indi = randi(length(indx), num_edge);
                    obj.adjacency(indx(indi), indy(indi)) = 1;
                    
                case -1
                    
                    [indx, indy] = find(funct_adj);
                    indi = randi(length(indx), num_edge);
                    obj.adjacency(indx(indi), indy(indi)) = 0;
            end
            
            obj.genotype_complexity = length(find(obj.adjacency(:)));
            
        end

        %% VAR II - add/delete nodes
        function [obj] = var_operator_2(obj, num_nodes, copdel)
            % copdel - add or remove
            % num_nodes - number of nodes to be added/deleted

            onn = obj.number_nodes;
            funct_adj = obj.adjacency(1:onn, 1:onn);
            
            rho = sum(funct_adj(:)>0) / numel(funct_adj);
            edge_dist = [sum(funct_adj) sum(funct_adj,2)'];
            % how many new edges are created?
            num_new_edges = edge_dist(randi(length(edge_dist), 2*num_nodes));            % 
            
            switch copdel
                case 1
                    for i = 1 : num_nodes
                        obj.adjacency(onn+i, randi(onn, num_new_edges(i))) = 1;
                        obj.adjacency(randi(onn, num_new_edges(i+1)), onn+i) = 1;
                        obj.adjacency(onn+i, onn+i) = rand(1) > (1-rho);
                        onn = onn+1;
                    end
                    
                case -1             
                    % cannot delete kernel (visible) units - only hidden
                    if num_nodes <= obj.number_nodes - obj.num_visible_nodes 
                        indx = obj.num_visible_nodes + randi(obj.number_nodes - obj.num_visible_nodes, num_nodes);
                    else
                        num_nodes2 = obj.number_nodes - obj.num_visible_nodes;
                        if num_nodes2 == 0
                            [obj] = var_operator_2(obj, num_nodes2, 1);
                        else
                            indx = obj.num_visible_nodes + randi(obj.number_nodes - obj.num_visible_nodes, num_nodes2);
                            obj.adjacency(indx,indx) = [];
                            obj.adjacency = [[obj.adjacency; zeros(indx, size(obj.adjacency,1))] zeros(size(obj.adjacency,1) + num_nodes, num_nodes)];
                        end
                        num_nodes = num_nodes2;
                    end

                                  
            end                    

            obj.genotype_complexity = length(find(obj.adjacency(:)));
            
            obj.number_nodes = obj.number_nodes + copdel * num_nodes;

        end  
        
        
        
        %% Var III - change in hidden unit topology by copying internal structure
        function [obj] = var_operator_3(obj, num_nodes, copdel)

        % define stucture here as connected subgraph
        % can only delete hidden structure, but can copy any structure -
        % this asymmetry is subtle but might be important....
        
            onn = obj.number_nodes;
            funct_adj = obj.adjacency(1:onn, 1:onn);
            
            switch copdel
                
                case 1                  % add substructure
                    
                    % define substructure to be copied
                    ind1 = randi(obj.number_nodes,1);
                    if num_nodes > 1
                        J = funct_adj^num_nodes;
                        jj = J(ind1,:); jj(ind1) = 0; indj = find(jj); ind2 = randi(length(indj), num_nodes-1);
                        ind1 = [ind1 indj(ind2)];
                    end
                                        
                    overlap = obj.adjacency(ind1,ind1);
                    
                    % randomly assign connections between copied structure and nodes copied from
                    
                    
                    
                    X1 = funct_adj(ind1,:);
                    X1 = reshape(X1(randperm(numel(X1))), size(X1));
                    
                    X2 = funct_adj(:,ind1);
                    X2 = reshape(X2(randperm(numel(X2))), size(X2));
                    
                    
                    obj.adjacency(onn+1 : onn + num_nodes, 1 : onn) = X1;
                    obj.adjacency(1 : onn, onn+1 : onn + num_nodes) = X2;
                    obj.adjacency(onn+1 : onn + num_nodes, onn+1 : onn + num_nodes) = overlap;   
                    
                case -1             % delete substructure
                    
                    % cannot delete kernel (visible) units - only hidden
                    if num_nodes <= obj.number_nodes - obj.num_visible_nodes 
                        ind1 = obj.num_visible_nodes + randi(obj.number_nodes - obj.num_visible_nodes, 1);          % have to ensure connected subgraph
                        indz = [[1:obj.num_visible_nodes] ind1];
                        if num_nodes > 1
                            J = funct_adj^num_nodes;
                            jj = J(ind1,:); jj(indz) = 0; indj = find(jj); ind2 = randi(length(indj), num_nodes-1);
                            ind1 = [ind1 indj(ind2)];
                        end 
                        
                    elseif num_nodes == 0
                        
                        [obj] = var_operator_3(obj, num_nodes, copdel);
                        
                    else
                        num_nodes = obj.number_nodes - obj.num_visible_nodes;
                        ind1 = obj.num_visible_nodes + 1 : obj.number_nodes;
                    end
                   
                    
                    if num_nodes ~= 0
                        obj.adjacency(ind1, :) = [];
                        obj.adjacency(:, ind1) = [];

                        % repad 
                            obj.adjacency = [[obj.adjacency; zeros(num_nodes, size(obj.adjacency,1) + num_nodes - 1)] ...
                            zeros(size(obj.adjacency,1) + num_nodes, num_nodes)];
                    end
                    
            end

            obj.genotype_complexity = length(find(obj.adjacency(:)));
            obj.number_nodes = obj.number_nodes + copdel * num_nodes;
            
        end
        
        
        %%
        function [obj] = var_operator_distributed(obj, numn, maxK, P, copdelP)
            % select above variational operators based on distribution P
            % form: structured vs unstructured
            % P - probability distribution over other operators in order
            % Pcopdel is the Pdist over delete/copy 
            % maxK - max length of walks through visible kernel - used to compute phenotype
            
            
            % in updated version - have to feed probabilities through to
            % here - sample and pass forward to mutate
            
            indy = RouletteWheelSelection(copdelP);
            if indy == 1, copdel = -1; elseif indy == 2, copdel = 1; end 
            
            indx = RouletteWheelSelection(P);
                        
            switch indx
                case 1
                    obj = var_operator_1(obj, numn, copdel);
                case 2
                    obj = var_operator_2(obj, numn, copdel);
                case 3
                    obj = var_operator_3(obj, numn, copdel);
            end       
            
            % update phenotype following mutation
                        
            [obj] = update_phenotype(obj, maxK);
            
            % update phenotypic complexity                        
            
        end

           

        %% NEW - to remove paths which cross through the initial/terminal node
        function [obj] = update_phenotype(obj, maxK)
            % phenotype is size obj.num_visible_nodes - each element corresponds to the number of paths of length <maxK that
            % start and finish at node i (visible nodes), but don't pass through i as an intermediate node

            onn = obj.number_nodes;
            funct_adj = obj.adjacency(1:onn, 1:onn);
            
            Q = zeros(obj.num_visible_nodes,1);
            for k = 2 : maxK

                no_internal_mats = k - 2;
                if no_internal_mats == 0
                    A = diag(funct_adj^2);
                else
                    A = zeros(obj.num_visible_nodes,1);
                    for i = 1 : obj.num_visible_nodes
                        At = funct_adj;
                        At(i,:) = 0;
                        At(:,i) = 0;                    
                        G = funct_adj * (At)^(k-2) * funct_adj;             % no T
                        A(i) = G(i,i);
                    end
                end
                Q = Q + A(1:obj.num_visible_nodes);
            end
            
            obj.phenotype = Q / sum(Q)
            
        end
        
        
        %% compute fitness
        function [obj, f] = compute_fitness(obj, funct, data )
            
            switch funct
                case 'kldiv'
                    f = exp(-real(sum((data+1e-8) .* log2((data+1e-8)./(obj.phenotype+1e-8)))));

                case 'l1norm'
                    % L1 norm 
                    f = 1 - sum(abs(obj.phenotype - data));
            end
            
            obj.num_hidden_nodes = obj.number_nodes - obj.num_visible_nodes;
            
            if f == Inf
                disp('fitness: inf')
            end
        end
        
        
        
        %% Phenotypic complexity
        function [H] = sample_phenotypic_complexity( obj, Ntr, Nm, sigma, P, copdelP )
        % Ntr - number of anchor points (centroids) in phenotype space - Each centers a Gaussian
        % Nm - number of mutations
        % hierarchical sampling under copdel and P.

        % dimension of data simplex    
        D = obj.num_visible_nodes;

        % anchor points within the simplex corresponding to the phenotype
        Q = abs(randn(D, Ntr));
        anchors = Q ./ repmat(sum(Q), [D 1]);
        
        % generate many candidate mutations - what are the corresponding phenotypes
        prop_m = zeros(D, Nm);
        dist = zeros(Nm, Ntr);

        for nm = 1 : Nm
            % for each proposed mutation
            org_copy = obj;
            org_copy = var_operator_distributed(org_copy, num_nodes, maxK, P, copdelP);
            org_copy = update_phenotype(obj, maxK);
            prop_m(:, nm) = org_copy.phenotype;            
            dist(nm, :) = vecnorm(anchors - repmat(prop_m(:,nm), [1 Ntr]));
        end


        f = @(x, sigma, dim) (1/sqrt(sigma*((2*pi)^dim)) * exp(-(1/(2*sigma))*x.^2));

        X = f(dist, sigma, D);
        samples = sum(X,2);

        % Sturges? rule is popular due to its simplicity. It chooses the number of bins to be ceil(1 + log2(numel(X))).
        % h1 = histogram(samples, 'BinMethod', 'sturges');

        % octave implementation
        h1 = hist(samples, 20, 1);
        H = -sum(h1 .* log2(h1));

        end
        
        
            
    end
end
            