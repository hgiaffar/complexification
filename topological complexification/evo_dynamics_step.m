function [population] = evo_dynamics_step(population, dynamics, varargin)

no_orgs = max(size(population));

pop_fitness = [population(1:no_orgs).fitness];          % vector of organism fitness

if mean(pop_fitness) == 0
    pop_fitness = pop_fitness + rand(size(pop_fitness));
end

% rng('shuffle')                                          % shuffle random no gen

switch dynamics
    
    case 'Moran'
        % selection cycle: one individual is selected for reproduction with probability proportional  
        % to fitness and one is selected for death uniformly. reproduction means copying with mutation.
   
        K = varargin{1};                    % number of k-cycles (birth/death selection) per timestep
        varop = varargin{2};
        numn = varargin{3};
        copdel = varargin{4};
                
        ind_death = randi(no_orgs, [K 1]);
        
        for ki = 1 : K
                         
            % death and reproduction, fitness proportionate selection
            ind_repro = RouletteWheelSelection(abs(pop_fitness));
            
%            ind_repro = find(pop_fitness==max(pop_fitness));
            repro_org = population(ind_repro).organism;
                        
            switch varop
                case 1
                    population(ind_death(ki)).organism = var_operator_1(repro_org, numn, copdel);
                                       
                case 2
                    population(ind_death(ki)).organism = var_operator_2(repro_org, numn, copdel);

                case 3                
                    population(ind_death(ki)).organism = var_operator_3(repro_org, numn, copdel);                    
                    
                case 4
                    maxK = varargin{5};
                    P = varargin{6};
                    population(ind_death(ki)).organism = var_operator_distributed(repro_org, numn, maxK, P, copdel);
            end
        end       
                      
end
        
end