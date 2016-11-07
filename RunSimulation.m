clear all;
clc;
rng('shuffle');

sl=makespecieslist;
numspecies=4;
year = 60;
fire_schedule=FireSchedule(8,2.5,year);

r_density=[.2601 .374 .187 .1789];

% Actual rain from 1985-2014
% rain=[19.94 4.64 8.84 6.72 5.83 8.26 14.85 23.34 8.21 22.86 10.23 13.57 31.02 9.26 10.17 15.50 4.24 10.32 8.59 26.76 10.75 3.02 9.74 8.13 12.42 17.86 7.60 6.92 4.5];  % 12 is the rain amount
% year = length(rain)

% Random distribution of rain values using log-normal distribution
% Add rain distribution information here
pd = makedist('Lognormal','mu',2.35,'sigma',0.5);
rain = random(pd,[1,year]);


% Simulation with the following parameters
% 20 by 10 poles for 80 by 40 m domain, species list 
% initial relative density, 5 plants per pole, 
% 0.1 m overlap allowed initially, fire schedule and rain schedule

% Run the simulation and create a movie of the simulation
[tots,its,cov] = spatial_sim(20,10,80,40,sl,r_density,5,.1,fire_schedule,rain);

% Run the simulation without generating figures
%[tots,its,cov] = spatial_sim_nofig(20,10,80,40,sl,r_density,5,.1,fire_schedule,rain);

