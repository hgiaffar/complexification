function [M] = random_fsparse_mat(n,density)

M = zeros(n);

row = find(rand((n^2), 1) > (1-density));

M(row) = 1;

% must ensure that there is at least one connection per column or corresponding row

while sum(sum(M)==0) ~= 0 || sum(sum(M,2)==0) ~= 0
 
    [ri, ci] = find(M);
    [rii, cii] = find(~M);
    indx = randi(length(ri), 1);
    indx2 = randi(length(rii), 1);

    M(ri(indx), ci(indx)) = 0;
    M(rii(indx2), cii(indx2)) = 1;

end



end