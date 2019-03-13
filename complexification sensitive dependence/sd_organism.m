classdef sd_organism
    
    % topological variation only
    
    properties
        % number of nodes in network
        number_nodes
        % adjacency matrix - 'genotype'
        adjacency
        % phenotype
        phenotype
        % genotypic complexity
        genotype_complexity
        % sample phenotypic complexity
        phenotype_complexity
        % inter phenotype distance / inter sample distance
        phenotype_complexity_dist
        % 1-D hole spectrum - number of paths of length k
        num_cycles
        % globality - measure of the distribution of node sensitivities
        
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
            obj.adjacency = [[A zeros(obj.number_nodes, full_size-obj.number_nodes)]; ...
                [zeros(full_size-obj.number_nodes) zeros(full_size-obj.number_nodes, obj.number_nodes)]];
            obj.genotype_complexity = length(find(obj.adjacency(:)));
            obj.num_visible_nodes = data_size;
            obj.phenotype = zeros(data_size,1);
            
        end
        
        
        %% VAR I - add/delete edges
        function [obj] = var_operator_1(obj, num_edge, copdel)
            
%             disp(sprintf('VO_I_%d', copdel))        
            
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
            
%             disp(sprintf('VO_II_%d', copdel))

            % copdel - add or remove
            % num_nodes - number of nodes to be added/deleted

            onn = obj.number_nodes;
            funct_adj = obj.adjacency(1:onn, 1:onn);
            
            rho = sum(funct_adj(:)>0) / numel(funct_adj);
            edge_dist = [sum(funct_adj) sum(funct_adj,2)'];
            % how many new edges are created?
            num_new_edges = edge_dist(randi(length(edge_dist), [2*num_nodes 1]));            
            
            switch copdel
                case 1
                    for i = 1 : num_nodes
                        obj.adjacency(onn+i, randi(onn, [num_new_edges(i) 1])) = 1;
                        obj.adjacency(randi(onn, [num_new_edges(i+1) 1]), onn+i) = 1;
                        obj.adjacency(onn+i, onn+i) = rand(1) > (1-rho);
                        onn = onn+1;
                    end
                                        
                case -1                      
                    % if no hidden nodes, add rather than delete
                    if  obj.num_hidden_nodes == 0
                        [obj] = var_operator_2(obj, num_nodes, 1);  
                        copdel = 1;
                    else
                        % cannot delete kernel (visible) units - only hidden
                        if num_nodes <= obj.num_hidden_nodes 
                            indx = obj.num_visible_nodes + randi(obj.num_hidden_nodes , [num_nodes 1]);
                            
                        % if num_hidden_nodes < num_nodes    
                        elseif num_nodes > obj.num_hidden_nodes 
                            num_nodes = obj.num_hidden_nodes;
                            indx = obj.num_visible_nodes + randi(obj.num_hidden_nodes, [num_nodes 1]);
                        end
                        
                        obj.adjacency(indx,:) = [];
                        obj.adjacency(:,indx) = [];
                        obj.adjacency = [[obj.adjacency; zeros(length(indx), size(obj.adjacency,2))] zeros(size(obj.adjacency,1) + length(indx), length(indx))];
                            
                    end
                                                   
            end                    

            obj.genotype_complexity = length(find(obj.adjacency(:)));
            
            obj.number_nodes = obj.number_nodes + copdel * num_nodes;
            obj.num_hidden_nodes = obj.number_nodes - obj.num_visible_nodes;

        end  
        
        
        
        %% Var III - change in hidden unit topology by copying internal structure
        function [obj] = var_operator_3(obj, num_nodes, copdel)
            
%         disp(sprintf('VO_III_%d', copdel))

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
                        copdel = 1;
                        
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
        function [obj] = var_operator_distributed(obj, numn, maxK, P, copdelP, calc_pc)
            % select above variational operators based on distribution P
            % form: structured vs unstructured
            % P - probability distribution over other operators in order
            % Pcopdel is the Pdist over delete/copy 
            % maxK - max length of walks through visible kernel - used to compute phenotype
            % calc_pc - controls calculation of phenotypic complexity - 1 compute, 0 don't
                        
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
            
%             if calc_pc
%                 [obj] = sample_phenotypic_complexity( obj, 100, 10, 0.1, P, copdelP, numn , maxK, data);
%             end    
            
        end

           

        %% NEW - to remove paths which cross through the initial/terminal node
        function [obj] = update_phenotype(obj, data)
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
            
            obj.phenotype = Q / sum(Q);
            
            D = numel(data);
            
            for k = 2 : D+1
                
                

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
        function [obj] = sample_phenotypic_complexity( obj, Ntr, Nm, sigma, P, copdelP, num_nodes, maxK, data )
        % Ntr - number of anchor points (centroids) in phenotype space - Each centers a Gaussian
        % Nm - number of mutations
        % hierarchical sampling under copdel and P.

        % dimension of data simplex    
        D = obj.num_visible_nodes;

        % sample points within the simplex corresponding to the phenotype
        Q = abs(randn(D, Ntr));
        anchors = Q ./ repmat(sum(Q), [D 1]);
                
        % generate many candidate mutations - what are the corresponding phenotypes
        dist1 = zeros(Nm, Ntr);
        dist2 = zeros(Nm, Ntr);
        orgphen = zeros(Nm,D);

        for nm = 1 : Nm
            % for each proposed mutation
            org_copy = var_operator_distributed(obj, num_nodes, maxK, P, copdelP, 0);
            org_copy = update_phenotype(org_copy, maxK);
            orgphen(nm,:) = org_copy.phenotype;
            dist1(nm, :) = vecnorm(anchors - repmat(org_copy.phenotype, [1 Ntr]), 1);
            dist2(nm, :) = vecnorm(anchors - repmat(org_copy.phenotype, [1 Ntr]), 2);
        end        
        
        f = @(x, sigma, dim) (1/sqrt(sigma*((2*pi)^dim)) .* exp(-(1/(2*sigma))*x.^2));
               
        X1 = f(dist1, sigma, D);
        samples = sum(X1,2);

        if rand(1) > 0.99
        figure, scatter3(orgphen(:,1),orgphen(:,2),orgphen(:,3), 10, 'r'), hold on,
                scatter3(anchors(1,:),anchors(2,:),anchors(3,:), 5, 'b'), hold on
                scatter3(obj.phenotype(1),obj.phenotype(2),obj.phenotype(3), 100, 'g', 'Filled'), hold on,
                scatter3(data(1),data(2),data(3), 100, 'c', 'Filled')
        figure, subplot(1,2,1), imagesc(X1), subplot(1,2,2), plot(samples)        
        end
        
        % Sturges? rule is popular due to its simplicity. It chooses the number of bins to be ceil(1 + log2(numel(X))).        
        sturges_bins = ceil(1 + log2(numel(samples)));

        % octave implementation
        h1 = hist(samples, sturges_bins, 1) + 1e-8;
        h1 = h1/sum(h1);
        Zn = -sum((ones(size(h1))/length(h1)) .* log2((ones(size(h1))/length(h1))));
        obj.phenotype_complexity = -sum(h1 .* log2(h1)) / Zn;
        
        %% should normalise inter phenotype distance by the average distance between sample points
        obj.phenotype_complexity_dist = mean(pdist(orgphen))/mean(pdist(anchors'));
        
        end        
        
        
        %% APPE and globality
        function [obj] = sample_adjacent_phenotype( obj, n_samples, P, copdelP, num_nodes, maxK, data )
        % compute adjacent possible phenotype entropy (Sp) and globality ( G - std of pw distance btw phenotype distribution )
        % n_samples - number of mutation samples from each node

        % dimension of data simplex    
        D = numel(data);
        
        % distribution over pairwise distances between nodes
        mean_pwdist = zeros(obj.number_nodes, 1);
        phenotype_list = [];

        for no = 1 : obj.number_nodes
            orgphen = zeros(n_samples, D);
            for nm = 1 : n_samples
                % for each proposed mutation
                org_copy = var_operator_distributed(obj, num_nodes, maxK, P, copdelP, 0);
                org_copy = update_phenotype(org_copy, maxK);
                orgphen(nm,:) = org_copy.phenotype;
            end
            
            mean_pwdist(no) = mean(pdist(orgphen'));
            phenotype_list = [phenotype_list; orgphen];
            
        end  
        
        % high globality score means small variance in btw-sampled phenotypes for each node
        org.globality = 1 / std(mean_pwdist) + 1e-8;
        
        % APPE is computed from eigenvalues of covariance matrix - in latent
        [~, ~, latent] = pca(phenotype_list');
        
        % log of the product of the positive eigenvalues
        org.phenotype_complexity = log(prod(latent(latent>0)));
        
        
%         if rand(1) > 0.99
%         figure, scatter3(orgphen(:,1),orgphen(:,2),orgphen(:,3), 10, 'r'), hold on,
%                 scatter3(anchors(1,:),anchors(2,:),anchors(3,:), 5, 'b'), hold on
%                 scatter3(obj.phenotype(1),obj.phenotype(2),obj.phenotype(3), 100, 'g', 'Filled'), hold on,
%                 scatter3(data(1),data(2),data(3), 100, 'c', 'Filled')
%         figure, subplot(1,2,1), imagesc(X1), subplot(1,2,2), plot(samples)        
%         end

        end        
        
        
        
            
    end
end
            