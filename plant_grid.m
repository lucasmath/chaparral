classdef plant_grid < handle
    %This class indexes plants and seeds to particular grids.  
    
    properties
        n % number of poles west to east
        m % number of poles north to south
        poles %total number of poles
        dimensions %size of grid length by width
        plot_dimensions %length by width of each plot
        pole_locations %where the poles are placed on the grid.
        pole_counts % number of plants at each pole
        plants %assigns plants to poles.  
        plant_locations %gives ordered list of where plants are indexed. 
        species_list
        plant_name
        plant_counts % number of each species
        seed_resp_counts %number of each species with seedling and resprouts
        tot_area %total area of each species
        tot_crown
        tot_height
        overlap
        dead_plants
    end
    
    methods
        %THIS FUNCTION CREATES A PLANT GRID WITH THE ABOVE PROPERTIES
        function obj = plant_grid(n, m, xdim, ydim, sl)
            obj.n = n; %number of poles west to east
            obj.m = m; %number of poles north to south
            obj.poles = m*n; %total number of poles
            obj.dimensions = [xdim, ydim]; %length by width of grid
            obj.plot_dimensions = [xdim/n, ydim/m]; %length by width of each plot
            obj.pole_locations = cell(m, n); %allocates memory for poles
            obj.pole_counts = zeros(m,n);
            obj.plants = cell(m, n); %allocate memory for plants
            obj.plant_locations = cell(100, 2);
            obj.species_list = sl; %makespecieslist;
            obj.plant_counts = zeros(1, length(obj.species_list));
            obj.seed_resp_counts = zeros(1, 2*length(obj.species_list));
            obj.tot_area = zeros(1, length(obj.species_list));
            obj.tot_crown = zeros(1, length(obj.species_list));
            obj.tot_height = zeros(1, length(obj.species_list));
            obj.dead_plants=cell(m,n); %saves a list of plants marked for death
            for i=1:m
                for j=1:n
                    obj.pole_locations{i, j} = [j*xdim/n - 0.5*xdim/n, i*ydim/m - 0.5*ydim/m];
                    obj.plants{i,j} = plant.empty(50,0);
                end
            end
            obj.plant_name=cell(1,length(sl));
            for i=1:length(sl)
                obj.plant_name{i} = sl{i}.name;
            end
            
        end
        
        %THIS FUNCTION RETURNS THE INTEGER POLE-LOCATION OF THE PLANT
        %X-COORD
        function ix = column(obj, p) %p is the plant
            ix = round((p.location(1) + 0.5*obj.plot_dimensions(1))/obj.plot_dimensions(1));
        end
        
        %THIS FUNCTION RETURNS THE INTEGER POLE-LOCATION OF THE PLANT
        %Y-COORD
        function iy = row(obj, p) %p is the plant
            iy = round((p.location(2) + 0.5*obj.plot_dimensions(2))/obj.plot_dimensions(2));
        end
        
        
        %THIS FUNCTION ADDS A SINGLE PLANT TO A POLE LOCATION
        function add_plant(obj, p)
            r = obj.row(p);
            c = obj.column(p);
            if r>0 && r<=obj.m && c>0 && c<=obj.n
                obj.pole_counts(r,c) = obj.pole_counts(r,c)+1;
                %num_plant = length(obj.plants{r, c});
                obj.plants{r, c}(obj.pole_counts(r,c)) = p;
            end
        end
        
        %THIS FUNCTION ADDS AN ARRAY OF PLANTS TO CORRESPONDING POLE
        %LOCATIONS
        function add_plants(obj, parray)
            % add an array of plants
            for k=1:length(parray)
               obj.add_plant(parray(k));
            end
        end      
        
        %THIS FUNCTION CALCULATES A LIST OF THE PLANTS AND THEIR LOCATIONS
        function plant_list(obj)
            u=1;
            for i=1:obj.m
                for j=1:obj.n
                    for k=1:obj.pole_counts(i,j)
                        obj.plant_locations{u, 1} = [i, j];
                        obj.plant_locations{u, 2} = [obj.plants{i, j}(k).location(1), obj.plants{i, j}(k).location(2)];
                        u=u+1;
                    end
                end
            end
            obj.plant_locations = obj.plant_locations(1:u-1, :);
        end
        
        %THIS FUNCTION GENERATES AN IMAGE OF THE PLANTS AND THE GRID
        function show_plants(obj)
            hold all
            grid on
            % go through every pole and every plant per pole
             for i=1:obj.m
                for j=1:obj.n
                   parray = obj.plants{i, j};           
                   for k=1:obj.pole_counts(i,j);
                      if (parray(k).status)
                        hci=round((parray(k).height)*30.0);
                        if hci>200
                            hci=200;
                        elseif hci<1
                            hci=1;
                        end
                         color = parray(k).color(hci,:);
                        Edge=[0 0 0];
                        if parray(k).adult==1
                            Edge=[1 1 1];
                        end
                        ellipsefill(parray(k).crown_east, parray(k).crown_north, parray(k).crown_west, parray(k).crown_south, parray(k).location(1), parray(k).location(2), color, Edge);
                      end
                   end
                end
             end

            axis([-.5 obj.dimensions(1)+.5 -.5 obj.dimensions(2)+.5]);
            hold off
        end
    
        function  resprout_plants(obj)        
            % Go through every pole
            for i=1:obj.m
                for j=1:obj.n
                   obj.dead_plants{i,j}=[];
                   % Go through every plant at pole i,j
                   if(~isempty(obj.plants{i, j}))
                        for k=1:obj.pole_counts(i,j)
                            index=find(ismember(obj.plant_name, obj.plants{i,j}(k).species));
                            obj.plants{i,j}(k).resprout_plant(obj.species_list{index}); 
                            if ~obj.plants{i,j}(k).status
                               obj.dead_plants{i,j}(end+1)=k;
                            end
                        end
                        obj.plants{i,j}(obj.dead_plants{i,j}) = [];
                        obj.pole_counts(i,j)=length(obj.plants{i,j});
                   end
                end
            end
        end
        
        function  grow_plants(obj, sb, rain) %, spk, rpk)        
            % Gathering statistics on area
            obj.tot_area = zeros(1,length(obj.species_list));
           
            % Go through every pole
            for i=1:obj.m
                for j=1:obj.n
                   % reset the list of dead plants to be empty
                   obj.dead_plants{i,j}=[]; 
                   % Go through every plant at pole i,j
                   if(~isempty(obj.plants{i, j})) 
                        for k=1:obj.pole_counts(i,j)
                            p = obj.plants{i,j}(k);
                            if ( p.adult || strcmp(p.species,'Cm') || strcmp(p.species,'Cs') )
                            % Go through all neighboring poles
                            for u=max(1, i-1):min(obj.m, i+1) 
                                for v=max(1, j-1):min(obj.n, j+1)
                                    % Go through plant at each neighbor pole           
                                    for w=1:obj.pole_counts(u,v)
                                       % Calculate dist once and use it twice 
                                       [dist,alpha] = p.distance(obj.plants{u,v}(w));
                                       if(dist > 0.001)
                                        % How do I make this more efficient?
                                        % ~(adult && seedling) = seedling || adult
                                        if (~p.adult || obj.plants{u,v}(w).adult)
                                            % Calculate collision with plant w
                                            obj.plants{i,j}(k).collide(obj.plants{u,v}(w),dist,alpha);
                                        end
                                        % Calculate competition with plant w
                                        if p.competition_factor
                                            obj.plants{i,j}(k).compete(obj.plants{u,v}(w),dist);
                                        end
                                       end
                                    end
                                end
                            end
                            end
                            
                            %  Check each plant for large negative
                            %  collisions that cause death
                            death_check=[obj.plants{i,j}(k).crown_east obj.plants{i,j}(k).crown_north obj.plants{i,j}(k).crown_west obj.plants{i,j}(k).crown_south]+obj.plants{i,j}(k).collision;
                            
                            flag1 = min(death_check)<(obj.plants{i,j}(k).crown_max/(-8.0));
                            
                            drought_tol = obj.plants{i,j}(k).get_drought_tolerance(rain); 
                           
                            flag2 = rand > drought_tol/(1.0 + obj.plants{i,j}(k).competition);
                            
                            if (flag1 || flag2) 
                                obj.dead_plants{i,j}(end+1)=k;
                            end
                        end
                   end
                end
            end
            
            %Set plant counts to 0
            for i=1:length(obj.plant_counts)
                obj.plant_counts(i) = 0;
                obj.tot_area(i)=0.0;
            end
            
            for i=1:obj.m
                for j=1:obj.n
                   % remove premarked dead plants
                   obj.plants{i,j}(obj.dead_plants{i,j}) = [];
                   obj.pole_counts(i,j)=obj.pole_counts(i,j)-length(obj.dead_plants{i,j});
                   % Go through every plant at pole i,j
                   if(~isempty(obj.plants{i, j}))
                        for k=1:obj.pole_counts(i,j)
                            
                            %Count plants
                            index=find(ismember(obj.plant_name, obj.plants{i,j}(k).species));
                            
                            if (obj.plants{i,j}(k).age == 6) 
                                if (~obj.plants{i,j}(k).adult)
                                    obj.plants{i,j}(k).maturation(obj.species_list{index});
                                end
                            end
                            
                            %Grow plant 
                            obj.plants{i,j}(k).grow(rain);
                            
                            %Add seeds to the seedbank from plant p, scale by age
                            r = round(obj.plants{i,j}(k).release_rate*(1.0 - exp(-0.5*obj.plants{i,j}(k).mature_age)));
                            if r 
                                sb.add_seeds(obj.plants{i,j}(k),r); 
                            end
                           
                            obj.plant_counts(index)=obj.plant_counts(index)+1;

                            %Add area
                            obj.tot_area(index)=obj.tot_area(index)+obj.plants{i,j}(k).area+obj.plants{i,j}(k).overlap_area;
                            obj.plants{i,j}(k).overlap_area=0.0;
                        end
                   end
                end
            end
        end
        
  
        function count_plants(obj) %This function counts seeds/resprouts separately
            num_species = length(obj.species_list);
            for i=1:num_species
                 obj.seed_resp_counts(i) = 0;
                 obj.seed_resp_counts(i+num_species) = 0;
                 obj.plant_counts(i) = 0;
                 obj.tot_area(i) = 0.0;
            end
            for i=1:obj.m
                for j=1:obj.n
                    for k=1:obj.pole_counts(i,j)
                        p = obj.plants{i,j}(k);
                        s = find(ismember(obj.plant_name, p.species));
                        if(~p.adult)
                            obj.seed_resp_counts(s) = obj.seed_resp_counts(s) + 1;
                        else
                            obj.seed_resp_counts(num_species+s) = obj.seed_resp_counts(num_species+s) + 1;
                        end
                        obj.plant_counts(s) = obj.plant_counts(s) + 1;
                        obj.tot_area(s)=obj.tot_area(s)+p.area;
                    end
                end
            end
        end
      
        %THIS FUNCTION CREATES N RANDOM PLANTS OF TYPE CS AND ASSINGS A DISTRIBUTION OF HEIGHTS AND CROWNS 
        function random_plants(obj,sb,numplants,overlap,r_density)  
            % Go through every pole
            prob=zeros(1,length(r_density));
            for i=1:length(r_density)
                prob(i)=sum(r_density(1:i));
            end
            for i=1:obj.m
                for j=1:obj.n
                    for z=1:numplants
                        k=rand;
                         for p=1:length(prob);
                            if (k<prob(p))
                                s=p;
                                break
                            end
                         end
                         species = obj.species_list{s};
                         c=-1.0; h=-1.0;
                         while c<0.0 || h<0.0
                            h = species.r_avg_height_growth+species.r_std_height_growth*randn; 
                            % fits were done for the crown diam so we mult by 1/2
                            c = 0.5*(species.r_avg_crown_growth+species.r_std_crown_growth*randn);
                         end
                         tot_collision = 10.0;
                         iter = 1;
                         while(tot_collision>overlap)
                            l=[obj.pole_locations{i,j}(1)+(rand-.5)*obj.plot_dimensions(1), obj.pole_locations{i,j}(2)+(rand-.5)*obj.plot_dimensions(2)];
                            temp_plant = plant(obj.species_list{s}, l, h, c, 23.5, 7); %temporarily gf is 23.6
                            temp_plant.adult = 1;
                            collision = zeros(1,4);
                            % Go through all neighboring poles
                            for u=max(1, i-1):min(obj.m, i+1) 
                                for v=max(1, j-1):min(obj.n, j+1)
                                % Go through plant at each neighbor pole
                                    for w=1:obj.pole_counts(u,v)
                                        [dist,alpha] = temp_plant.distance(obj.plants{u,v}(w));
                                        collision=collision-temp_plant.collide(obj.plants{u,v}(w),dist,alpha);          
                                        %temp_plant.collide(obj.plants{u,v}(w),dist,alpha)
                                    end
                                end
                            end
                            
                            if(iter < 101) 
                                tot_collision = sum(collision);
                                iter = iter+1;
                            else
                                tot_collision = 0.0; 
                            end             
                         end
                         
                         if iter < 101
                            r = round(temp_plant.release_rate);%/obj.species_list{s}.predation_rate);
                            if r
                                sb.add_seeds(temp_plant,r);
                            end 
                            p = obj.pole_counts(i,j)+1;
                            obj.pole_counts(i,j) = p;
                            obj.plants{i,j}(p) = temp_plant; 
                                    
                            %Count plants
                            index=find(ismember(obj.plant_name, obj.plants{i,j}(p).species));
                            obj.plant_counts(index)=obj.plant_counts(index)+1;

                            %Add area
                            obj.tot_area(index)=obj.tot_area(index)+obj.plants{i,j}(p).area;
                           
                         end
                    end
                end
            end
        end
                  
    end
    
end

