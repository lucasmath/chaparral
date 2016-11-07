classdef Cm < plant_type
    % sub-class of plant_type which contains all info about 
    % Ceanothus megacarpus species
    
    properties (SetAccess = private)
    end
    
    methods
        function obj = Cm()
        % Arguments to Ceanothus megacarpus constructor are
        %  
        %  'Cm'     name as string
        %  'NS'     life history type as string%
        %   0       reprouter = 0
        
        %  0.2301   mean seedling height growth parameter
        %  0.0602   std seedling height growth
        %  0.3237   mean seedling crown growth parameter  
        %  0.0931   std seedling crown growth 
       
        %  -1.439   seedling mean initial height parameter
        %  0.3040   seedling std initial height parameter
        %  -2.621   seedling mean initial crown parameter
        %  0.4701   seedling std initial crown parameter
        
        %  Borrowed from Cs to be used for adult growth
        %  0.07210   average resprout height growth
        %  0.02400   std resprout height
        %  0.0907    average resprout crown growth
        %  0.03500   std resprout crown
        
        %  [1 0 0] species color (red)
        
        %  2000      release rate
        %  0.9       predation rate
        %  0.0195    germination rate
        
        %  0.002    seedling competition factor
        %  0.0        resprout competition factor
        %  [1.0 0.8 0.0 0.0] drought tolerance  
        
            % Ceanothus megacarpus constructor
            obj = obj@plant_type('Cm','NS',0,0.2301,0.0602,0.3237,0.0931,-1.439,0.3040,-2.621,0.4701, 0.07210, 0.02400, 0.0907, 0.03500,[1 0 0],2000,0.9,0.0195,0.002,0.0,[1.0 0.8 0.0 0.0]);
            
        end
        
        end
            
end