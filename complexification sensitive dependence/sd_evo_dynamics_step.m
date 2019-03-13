function [population] = evo_dynamics_step(population, dynamics, varargin)

no_orgs = max(size(population));

pop_fitness = [population(1:no_orgs).fitness];          % vector of organism fitness

if mean(pop_fitness) == 0
    pop_fitness = pop_fitness + rand(size(pop_fitness));
else
    pop_fitness = pop_fitness - min(pop_fitness);
    pop_fitness = pop_fitness / max(pop_fitness);
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
            if 1
                ind_repro = RouletteWheelSelection(abs(pop_fitness/sum(pop_fitness)));            
            else        % distort the above probabilities
                selection = varargin{7};
                select_pop_fitness = (exp(selection * pop_fitness) / sum(exp(selection * pop_fitness)));
                ind_repro = RouletteWheelSelection(abs(select_pop_fitness));
                
            end
            % strongest possible selection
            ind_repro = find(pop_fitness==max(pop_fitness));

%                 find([population(1:length(population)).fitness] == max([population(1:length(population)).fitness]))
%                 ind_repro
%                 find([population(1:length(population)).fitness] == max([population(1:length(population)).fitness])) == ind_repro
                                                
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
                    calc_pc = varargin{8};
                    population(ind_death(ki)).organism = var_operator_distributed(repro_org, numn, maxK, P, copdel, calc_pc);
%                     population(ind_death(ki)).lineage = [population(ind_death(ki)).lineage population(ind_repro).lineage];
            end
        end       
                      
end
        
end