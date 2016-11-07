classdef seedbank < handle
   % A list of all seeds with locations and species
    
    properties (SetAccess = private)
        xdim;
        ydim;
        plant_species
        num_species
        plant_name
        plant_type
        location                % x and y coordinates
        predation_rate = 0.9;
        num_seeds;
    end    
    methods
        function obj = seedbank(x,y,species_list)
            obj.xdim = x;
            obj.ydim = y;
            obj.num_species = length(species_list);
            obj.plant_species = species_list;
            % cell array with matrices of seed locations
            obj.location = cell(1,obj.num_species);         
            for i=1:obj.num_species
                obj.plant_name{i} = species_list{i}.name;
                obj.plant_type{i} = species_list{i}.type;
                obj.location{i} = zeros(50,2);
            end
            obj.num_seeds = zeros(1,obj.num_species);
        end
        
        function add_seed(obj,sp,xy)
          % find index of plant species (p)
          index=find(ismember(obj.plant_name, sp));
          % j = row of the last seed in location array
          j = obj.num_seeds(index);
          
          obj.location{index}(j+1,:)=xy; %new seed added
          obj.num_seeds(index) = j+1; 
        end
        
        function add_seed_list(obj,sp,xy)
          % find index of plant species (p)
          index=find(ismember(obj.plant_name, sp));
          % j = row of the last seed in location array
          j = obj.num_seeds(index);
          s = size(xy,1);
          
          obj.location{index}(j+1:j+s,:)=xy; %new seed added
          obj.num_seeds(index) = j+s; 
        end
            
        
        function add_seeds(obj,p,r) % inputs a plant
          sp=p.species; 
          % find index of plant species (p)
          index=find(ismember(obj.plant_name, sp));
          % get location of plant p
          px = p.location(1);
          py = p.location(2);
          % j = row of the last seed in location array
          j = obj.num_seeds(index);
          
          % older exponential model
          % get mean seed dist for plant p   
          % mean_dist = (-log(.65))/p.crown_max; 
          % generate random numbers for radius
          % zr=-log(rand(r,1))/mean_dist; 
          
          % 34% under the canopy
          zr = zeros(r,1);
          u = rand(r,1);
          indices = find( u < 0.34 );
          n = length(indices);
          zr(indices) = p.crown_max*rand(n,1);
          indices = find( u > 0.34 );
          n = length(indices);
          zr(indices) = -log(rand(n,1))/0.423;
          
          ztheta = 2*pi*rand(r,1); %generates random numbers for angle 
          xy=zeros(r,2); %preallocating seed coordinate space
         
          xy(:,1)= px + zr.*cos(ztheta); %xcoordinate
          xy(:,2)= py + zr.*sin(ztheta); %ycoordinate
          
          obj.location{index}(j+1:j+r,:)=xy; %new seeds added
          obj.num_seeds(index) = j+r; 

        end
        
        % Add seeds from plant array
        function add_seeds_from_parray(obj, parray)
            %add an array of plants
            for i=1:length(parray)
                obj.add_seeds(parray(i))
            end
        end
        
        % Add seeds from plant grid
        function add_seeds_from_grid(obj, g)
            for i=1:g.m
                for j=1:g.n
                    if(~isempty(g.plants{i, j}))
                        for k=1:g.pole_counts(i,j)
                            add_seeds(obj, g.plants{i, j}(k));
                        end
                    end
                end
            end
        end
        
        % Get rid of 90% of previously accumulated seeds
        function predate_seeds(obj)
            rho = obj.predation_rate;
            for i=1:obj.num_species
                z = ceil(rand(1,obj.num_seeds(i)) - rho);
                index = find(z);
                obj.location{i} = obj.location{i}(index,:);
                obj.num_seeds(i) = length(index);
            end
        end
       
       %THIS FUNCTION IS FOR FIRE-STIMULATED SEEDS.  IT ADDS THE SEEDS IN
       %THE SEEDBANK--WHICH ARE ALREADY ASSUMED TO HAVE BEEN PREDATED AND
       %GERMINATED--AS PLANTS ON THE GRID.  IT IS TO BE USED WHEN FIRE
       %OCCURS
       function fire_germination(obj, g)
           %for i=1:(obj.first_os-1)
           for i=1:obj.num_species
                    % numseeds = obj.num_seeds(i)
                    % obj.location{i}
                    for j=1:obj.num_seeds(i)
                            p = plant(obj.plant_species{i},[obj.location{i}(j, 1), obj.location{i}(j, 2)]);
                            add_plant(g, p);
                    end
                obj.num_seeds(i) = 0;
           end
       end
       
       %THIS FUNCTION COUNTS THE NUMBER OF SEEDS FOR EACH SPECIES
       function count_seeds(obj)
           for i=1:length(obj.plant_species)
               obj.num_seeds(i) = length(obj.location{i});
           end
       end
       
    end
                
end

