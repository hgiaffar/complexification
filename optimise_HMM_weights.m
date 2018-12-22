function [obj] = optimise_HMM_weights(obj, data, steps, lr, eps, beta)

% backpropogate weights through the eigenvector equation and take steps to
% minimise error in the model

% steps - number of steps taken in direction of desired soln
% lr - length of step

M = obj.transition_matrix;

visible_nodes = max(size(data));
hidden_nodes = size(M,1) - visible_nodes;

% error vector
e = (data - obj.stat_dist(1:visible_nodes))/sum(obj.stat_dist(1:visible_nodes));

e = [e; randn(hidden_nodes,1)*0.001];

% initialise guess error in M
Mtil = randn(size(M));

if max(size(obj.stat_dist_full)) ~= max(size(Mtil))
    obj.stat_dist_full = [obj.stat_dist_full; 0];
end

D = ((eye(size(M,1)) - M) * e) \ (Mtil * obj.stat_dist_full);
% D = ((eye(size(M,1)) - M) * e) ./ (Mtil * obj.stat_dist_full);

Mt = Mtil .* repmat(D, [1 size(M,2)]);

X = max(M + (lr * Mt), 1e-6);
% X = max(M + (lr * Mt), 0);

for i = 1 : size(M,2)
    X(:,i) = X(:,i) / sum(X(:,i));
end

% log(X)
% 
% obj.weights = obj.weights + log(X);
% obj = update_transitions(obj, 0, eps, beta);

obj.transition_matrix = X;


end