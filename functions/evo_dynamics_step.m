function [population] = evo_dynamics_step(population, dynamics, varargin)

no_orgs = max(size(population));

pop_fitness = [population(1:no_orgs).fitness];          % vector of organism fitness

rng('shuffle')                                          % shuffle random no gen

switch dynamics
    
    case 'Moran'
        % selection cycle: one individual is selected for reproduction with probability proportional  
        % to fitness and one is selected for death uniformly. reproduction means copying with mutation.
   
        K = varargin{1};                    % number of k-cycles (birth/death selection) per timestep
        eta = varargin{2};
        eps = varargin{3};
        beta = varargin{4};
        
        ind_death = randi(no_orgs, [K 1]);
        
        for ki = 1 : K
             
            % death and reproduction
            % fitness proportionate selection
            ind_repro = RouletteWheelSelection(pop_fitness);
            
%             ind_repro = find(pop_fitness==max(pop_fitness));
                
            repro_org = population(ind_repro).organism;
                        
%             population(ind_death(ki)).organism = var_operator_1(repro_org, eta, eps, beta);
            population(ind_death(ki)).organism = var_operator_2(repro_org, 'structured', 1, eps, beta);

        end       
                      
end
        
end