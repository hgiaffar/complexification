function [M] = generate_nonbinary_hierarchical_envt(size_envt, b, l, no_samples, sigma, beta, type)
% size_envt - number of leaf nodes - size of environment
% b - branch number
% l - number of levels in the hierarchy
% sigma - sd of gaussian P distribution over edge 
% type - type of graph topology (NOT output)
% beta - controls softmax

if isempty(beta), beta = 1; end

%check
if b^l ~= size_envt
    disp('error')
    return
end

edge_mean = 0;              % mean value for each edge is zero

M = zeros(no_samples, size_envt);

if type == 'binary'

    for ntr = 1 : no_samples

        T = zeros(l, size_envt);
        for i = 1 : l
            g = normrnd(0, sigma, [b^i 1]);            % draw a series of edge weights from a normal dist.  can be changed s.t. each edge is contolled by its own sigma
            T(i, :) = repelem(g, size_envt/(2^i));
        end
        
        % softmax over output
        M(ntr,:) = exp(beta*sum(T)) / sum(exp(beta*sum(T)));   
    end
    
% allow struture to be randomly assigned, consistent with depth, average branch number and number of outputs
elseif type == 'non-binary'
    
    % construct matrix: rows level in hierarchy, columns nodes at level, 
    % entries out degree of the node.  branching number should equal average out degree of nodes
    
    for ntr = 1 : no_samples
        
        
        
    end
    
end

end
