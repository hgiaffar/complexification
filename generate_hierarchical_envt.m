function [M] = generate_hierarchical_envt(size_envt, b, l, no_samples)
% size_envt - number of leaf nodes - size of environment
% b - branch number
% l - number of levels in the hierarchy

%check
if b^l ~= size_envt
    disp('error')
    return
end

% generate graph
G = cell(l,1);

for i = 1 : l
    G{i} = rand(b^l, 2);            % sample from uniform (p(o==1 | i==0) and p(o==1 | i==1))
end

% sample:

M = zeros(no_samples, size_envt);

for ntr = 1 : no_samples
    for j = 1 : l + 1
        if j == 1
            N = rand(1)>0.5;
        else
            tN = zeros(1,b^(j-1));
            for k = 1 : length(N)            
                if N(k) == 0
                    tN((k-1)*b+1 : k*b) = rand(1,b) > G{j-1}((k-1)*b+1 : k*b);          
                else
                    tN((k-1)*b+1 : k*b) = rand(1,b) > 1 - G{j-1}((k-1)*b+1 : k*b);
                end
            end
            N = tN;
        end 
    end
    
    N(find(N==0)) = 1e-3;
    N = N/sum(N);
    
    M(ntr,:) = N;

end

end
