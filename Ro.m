classdef Ro < plant_type
    % sub-class of plant_type which contains all info about the
    % Rhus ovata species
    
    properties (SetAccess = private)
    end
    
    methods
        function obj = Ro()
        % Arguments to Rhus ovata constructor are        
        %
        % 'Ro'      name as string
        % 'FS'      life history type as string
        %           resprouter = 1
        
        %  0.1914    mean seedling height growth parameter
        %  0.0698    std seedling height growth
        %  0.2378    mean seedling crown growth parameter 
        %  0.0978    std seedling crown growth 
        %            ******taken from css data
                
        %  -1.493    seedling avg initial height parameter
        %  0.3522    seedling avg std initial height parameter
        %  -2.540    seedling avg initial crown parameter
        %  0.4938    seedling avg std crown parameter
        %            ******taken from css data       
        
        %  0.0545   resprout mean height growth parameter
        %  0.0282   resprout standard height growth parameter
        %  0.0884   resprout mean crown growth parameter
        %  0.0367   resprout standard crown growth parameter
        
        %  [1 0 1]  species color (magenta)
        
        %  2000     release rate
        %  0.9      predation rate
        %  0.15     germination rate
            
        %  0.021    seedling competition factor
        %  0.000    resprout competition factor
        %  [0.002 0.001 1.0 1.0] drought tolerance
        
            % Rhus ovata constructor
            obj = obj@plant_type('Ro', 'FS', 1, .1914, .0698, .2378, .0978, -1.493, .3522, -2.540, .4938, 0.0545, 0.0282, 0.0884, 0.0367, [1 0 1], 2000, 0.9, 0.15, 0.000, 0.0, [0.002 0.001 1.0 0.99]);  
        end
    end
end
