function [M] = random_fsparse_mat(n,density)

M = zeros(n);

row = find(rand((n^2), 1) > (1-density));

M(row) = 1;


end