classdef Cs < plant_type
    % sub-class of plant_type which contains all info about 
    % Ceanothus spinosus species
    
    properties (SetAccess = private)
    end
    
    methods
        function obj = Cs()
        % Arguments to constructor are
        %    
        %  'Cs'       name as string
        %  'FS'       life history type as string
        %   1         reprouter = 1
        
        %  0.1914    mean seedling height growth parameter
        %  0.0698    std seedling height growth
        %  0.2378    mean seedling crown growth parameter 
        %  0.0978    std seedling crown growth 
        
        %  -1.493    seedling avg initial height parameter
        %  0.3522    seedling avg stdv initial height parameter
        %  -2.540    seedling avg initial crown parameter
        %  0.4938    seedling avg initial crown parameter
        
        %  0.0721   mean resprout height growth
        %  0.0192   std resprout height
        %  0.0907   mean resprout crown growth
        %  0.0231   std resprout crown
           
        %  [0 0 1] species color (blue)
       
        %  2000      release rate
        %  0.9       predation rate
        %  0.0396    germination rate
        
        %  0.002     seedling competition factor
        %  3.5      resprout competition factor
        %  [1.0 0.60 1.0 0.95] drought tolerance

           % Ceanothus spinosus constructor
           obj = obj@plant_type('Cs','FS', 1, .1914, .0698, .2378, .0978, -1.493, 0.3522, -2.540, 0.4938, 0.07210, 0.0192, 0.0607, 0.0350, [0 0 1], 2000, 0.9, 0.0396, 0.002, 0.0375, [1.0 0.55 1.0 0.95]);
        end
       
    end
end