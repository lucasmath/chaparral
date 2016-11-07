classdef plant_type < handle
    %This Class calls a specific species and gives its various parameters
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        name % string
        type % string NS, FS, OS
        resprouter  % 1 if resprouter, 0 if nonsprouter
            
        % growth parameters for seedlings
        avg_height_growth
        std_height_growth
        avg_crown_growth
        std_crown_growth
        
        % initial parameters for seedlings
        s_avg_height
        s_std_height
        s_avg_crown
        s_std_crown
       
        % growth parameters for resprouts
        r_avg_height_growth
        r_std_height_growth
        r_avg_crown_growth
        r_std_crown_growth
        
        color
        
        release_rate 
        predation_rate
        germination_rate
        
        effective_release_rate
        longterm_release_rate
        
        seedling_competition_factor
        resprout_competition_factor
        drought_tolerance    
    end
    
    methods
        function obj = plant_type(n,t,r,ahg,shg,acg,scg,sah,ssh,sac,ssc,rahg,rshg,racg,rscg,c,rr,pr,gr,scf,rcf,dt)
            obj.name = n;    %species
            obj.type = t;    %life history type
            obj.resprouter = r;
            obj.avg_height_growth = ahg;
            obj.std_height_growth = shg;
            obj.avg_crown_growth = acg;
            obj.std_crown_growth = scg;
            obj.s_avg_height = sah;
            obj.s_std_height = ssh;
            obj.s_avg_crown = sac;
            obj.s_std_crown = ssc;
            obj.r_avg_height_growth = rahg;
            obj.r_std_height_growth = rshg;
            obj.r_avg_crown_growth = racg;
            obj.r_std_crown_growth = rscg;
            obj.release_rate = rr;
            obj.predation_rate = pr;
            obj.germination_rate = gr;
            obj.effective_release_rate = gr*(1.0-pr)*rr;
            obj.longterm_release_rate = gr*(1.0-pr)*rr/pr;
            obj.seedling_competition_factor=scf;
            obj.resprout_competition_factor=rcf;
            obj.drought_tolerance=dt;
            R=fliplr(linspace(0,c(1),200))';
            G=fliplr(linspace(0,c(2),200))';
            B=fliplr(linspace(0,c(3),200))';
            colormap=horzcat(R,G,B);
            obj.color=colormap;
        end
            
    end
        
end


