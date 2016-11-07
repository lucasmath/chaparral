function  [totals,inits,tot_percentarea_yr] = spatial_sim_nofig(n, m, xdim, ydim, species_list, r_density, numplants, collision, fire_schedule, rain)
% Seed the random number generator using the clock
rng('shuffle');

% make a plant grid with n by m poles, xdim by ydim dimensions
f=plant_grid(n, m, xdim, ydim, species_list);
% create a seedbank of size xdim by ydim
sb=seedbank(xdim, ydim, species_list);
% total time of simulation
T = length(fire_schedule);
% total number of fires
num_fires = sum(fire_schedule);
% number of species
numspecies = length(species_list);

%Statistics collected for each year
num_species_yr=zeros(T+1+num_fires,numspecies);
num_seed_resp_yr=zeros(T+num_fires,2*numspecies);
% total area of each species for each year
species_tot_area_yr=zeros(T+1,numspecies);
% percent ground cover for each year
tot_percentarea_yr=zeros(T+1,1);
% average crown radius by species for each year
average_crown_radius_yr=zeros(T+1,numspecies);
% average height by species for each year
average_height_yr=zeros(T+1,numspecies);
% relative density of each species for each year
species_percentarea_yr=zeros(T+1,numspecies);
% total crown area
tot_crown_yr=zeros(T+1,numspecies);

f.random_plants(sb, numplants, collision, r_density);
f.count_plants;
num_species_yr(1,:)=f.plant_counts;
inits = num_species_yr(1,:);
init_percents = inits/sum(inits)
num_seed_resp_yr(1,:)=f.seed_resp_counts;
species_tot_area_yr(1,:)=f.tot_area;
species_percentarea_yr(1,:)=f.tot_area/(xdim*ydim);
average_crown_radius_yr(1,:)=f.tot_crown./f.plant_counts;
average_height_yr(1,:)=f.tot_height./f.plant_counts;
tot_percentarea_yr(1)=sum(f.tot_area)/(xdim*ydim);

totals=num_species_yr(1,:);

for y=1:T 
    if (fire_schedule(y)==1)
        f.resprout_plants;
        sb.fire_germination(f);
        f.count_plants();
        num_seed_resp_yr(y+sum(fire_schedule(1:y)),:)=f.seed_resp_counts;
        num_species_yr(y+sum(fire_schedule(1:y)),:)=f.plant_counts; 
        average_crown_radius_yr(y+sum(fire_schedule(1:y)),:)=f.tot_crown;
    else
        sb.predate_seeds;
    end
      f.grow_plants(sb,rain(y));
      
      f.count_plants;     
      num_species_yr(y+sum(fire_schedule(1:y))+1,:)=f.plant_counts;
      totals(y+1,:) = num_species_yr(y+sum(fire_schedule(1:y))+1,:);
      num_seed_resp_yr(y+sum(fire_schedule(1:y))+1,:)=f.seed_resp_counts;
      species_tot_area_yr(y+1,:)=f.tot_area;
      species_percentarea_yr(y+1,:)=f.tot_area/(xdim*ydim);
      tot_percentarea_yr(y+1)=sum(f.tot_area)/(xdim*ydim);
      average_crown_radius_yr(y+sum(fire_schedule(1:y))+1,:)=f.tot_crown./f.plant_counts;
      average_height_yr(y+1,:)=f.tot_height./f.plant_counts;
end


savefile='simulation.mat';
save(savefile,'fire_schedule', 'rain', 'species_list', 'r_density', 'num_species_yr','species_tot_area_yr', 'tot_percentarea_yr', 'average_crown_radius_yr','average_height_yr', 'species_percentarea_yr');

end



