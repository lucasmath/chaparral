  classdef Ml < plant_type
    % sub-class of plant_type which contains all info about the
    % Malosma laurina species
    
    properties (SetAccess = private)
    end
    
    methods
        function obj = Ml()
        % Arguments to Malosma laurina constructor are    
        %
        %  'Ml'      name as string
        %  'FS'      life history type as string
        %   1        resprouter = 1
        
        %  0.1914    mean seedling height growth parameter
        %  0.0698    std seedling height growth
        %  0.2378    mean seedling crown growth parameter 
        %  0.0978    std seedling crown growth 
        %            ******taken from Css data
        
        %  -1.493    seedling avg initial height parameter
        %  0.3522    seedling std initial height parameter
        %  -2.540    seedling avg initial crown parameter
        %  0.4938    seedling std initial crown parameter
        %            *******taken from Css data
        
        %  0.0663   resprout mean height growth parameter
        %  0.0140   resprout standard height growth parameter
        %  0.0774   resprout mean crown growth parameter
        %  0.0347   resprout standard crown growth parameter
        
        %  [0 1 0]  species color (green)
        
        %  2000     release rate
        %  0.9      predation rate
        %  0.18     germination rate
        
        %  0.021    seedling competition factor
        %  0.000    resprout competition factor
        %  [0.002 0.001 1.0 1.0] drought tolerance
   
            % Malosma laurina constructor
            obj = obj@plant_type('Ml', 'FS', 1, .1914, .0698, .2378, .0978, -1.493, .3522, -2.540, .4938, 0.0663, 0.0140, 0.0774, 0.0347, [0 1 0], 2000, 0.9, 0.18, 0.000, 0.0, [0.002 0.001 1.0 0.99]);
            
        end
    end
end
