function [H] = sample_phenotypic_complexity( population, Ntr, sigma )

no_orgs = length(population);
siz_data = length(population(1).organism.stat_dist);
Z = zeros(no_orgs, siz_data);
D = zeros(Ntr, no_orgs);

% sampling points
Q = randn(siz_data, Ntr);
Q = Q ./ repmat(sum(Q), [siz_data 1]);

for i = 1 : no_orgs
    % centroids defined by organism fitness
    Z(i,:) = population(i).organism.stat_dist;
    % distance between centroids and sampling points
    D(:,i) = vecnorm(Q - repmat(Z(i,:)', [1 Ntr]));
end

f = @(x, sigma, dim) (1/sqrt(sigma*((2*pi)^dim)) * exp(-(1/(2*sigma))*x.^2));

X = f(D, sigma, siz_data-1);

samples = sum(D,2);

% Sturges? rule is popular due to its simplicity. It chooses the number of bins to be ceil(1 + log2(numel(X))).
h1 = histogram(samples, 'BinMethod', 'sturges');

vals = (h1.Values+1e-8)/sum(h1.Values+1e-8);

H = -sum(vals .* log2(vals));
close

end