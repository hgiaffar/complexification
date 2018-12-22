% create dtmc

N = 28;
M = 28;

states = N*M;

T = rand(states);

eps = 1e-7;
beta = 100;

V = zeros(size(T));

tic
for j = 1 : states
    su = sum(exp(beta * T(:,j)));
        
    for i = 1 : states
        V(i,j) = ((1-eps) * exp(beta * T(i,j)) / su ) + eps; 
    end 
    
end
toc

figure, subplot(1,2,1), imagesc(T), colorbar, subplot(1,2,2), imagesc(V), colorbar

[Evec, Eval] = eigs(V);
K = real(Evec(:,1));        % stationary distributino of dtmc is equal to the eigenvector corresponding to the largest eigenvalue

phenotype = K / sum(K);

% calculate fitness as DKL btw phenotype and next data point.

[KL] = KLdiv(q,p, 2)


