function [M] = random_fsparse_mat(n,density)

% in this version - remove self loops

M = rand(n) > (1-(density + density/n));
M = tril(M, -1) + triu(M,1);

% must ensure that there is at least one connection per column or corresponding row
while sum(sum(M)==0) ~= 0 || sum(sum(M,2)==0) ~= 0
 
    [ri, ci] = find(M);
    [rii, cii] = find(~M);
    % remove self loops from consideration
    rii = rii(~(rii==cii));
    cii = cii(~(cii==rii));
    indx = randi(length(ri), 1);
    indx2 = randi(length(rii), 1);

    M(ri(indx), ci(indx)) = 0;
    M(rii(indx2), cii(indx2)) = 1;

end



end