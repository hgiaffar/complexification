% simulation
T = 400000;              % length of simulation
et = 5000;              % number of steps during which enviroment is constant
no_orgs = 100;          % number of organisms
Kc = 1;                 % number of k-cycle selections per timestep
no_data_points = 16;    % size of the data
no_hidden_units = 0;    % initial number of hidden nodes
graph_density = 0.5;    % initial graph density
maxK = 5;               % maximum length path considered in the calculation of phenotype (~ A^k)
slct_strength = 100;    % fitness distribution can be sharpened (have to change in evo_dynamics_step)
preall_size = 200;      % PREALLOCATE MEMORY

% mutation
varop = 4;              % variational operator
P = [0.6 0.2 0.2];      % probability distribution over variational operators
copdelP = [0.5 0.5];    % distribution over copy vs delete.
numn = 1;               % number of edges or nodes or size of substructure copied/deleted

% sampling phenotypic complexity
calc_pc = 0;            % compute phenotypic complexity as part of the Moran cycle step
num_samples_pc = 10;    % how many anchors in the phenotype simplex used to approximate the distribution          
num_mutations_pc = 10;  % number of phenotypes sampled through random mutation for a given genotype
sigma_pc = 0.1;         % variance of gaussian centered at each anchor point

% generative model of the environment
envt_type = 'complex';          % predefined variables for environment generation
no_samples = ceil(T/et);        % number of environments generated from the hierarchical generative model
tree_type = 'binary';           % type of tree structure  
sigma_envt = 0.2;               % variance of the gaussian centered at 0 for each edge in the graph
beta_envt = 2;                  % softmax parameter

% report during simulation
reporting_interval = 1000;      % for updating plots and calculating phenotypic complexity
check_adjacency = 100;          % check that the number of nodes isn't near the preallocated space for adjacency matrices
add_preallocated_nodes = 50;    % if it is, by threshold measure below, add this many nodes
num_nodes_thresh = 0.75;        
clear_leak_interval = 5000;    % clear parameters and reload 