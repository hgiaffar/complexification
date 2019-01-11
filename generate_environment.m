function [environment] = generate_environment(siz_data, type, num_data)

% generate centoids of the environment.  These can be randomised within the
% main program

switch type
    case 'complex'
        % no structure - random walk in the simplex
        
        Q = abs(randn(siz_data, num_data));
        environment = Q ./ repmat(sum(Q), [siz_data 1]);
        
        
    case 'simple'
       % include a constraint - i.e. constrain two coordinates to be equal
        
       environment = zeros(siz_data, num_data);
       environment(1,:) = rand(1, num_data);
       for i = 1 : num_data
           environment(2:3,i) = (1-environment(1,i))/2;
       end
       
       
end

end