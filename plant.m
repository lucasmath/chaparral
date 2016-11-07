classdef plant < handle
   % Individual plants with attributes and functions for growth/seed
   % dispersal, etc
    properties %(SetAccess = private)
        species  %plant_type class
        age % age/time since previous wildfire
        mature_age %
        resprouter % capable of resprouting
        adult % is a resprout or reproductively mature
        location %vector with xcoordinate and ycoordinate
        height %in m
        % crown radii, allows for elliptical plants
        crown_north
        crown_south
        crown_east
        crown_west
        crown_max
        crown_avg
        area 
        overlap_area
        
        % growth parameters
        height_growth
        crown_growth
        
        % competition factors
        collision
        competition
        competition_factor
        competition_distance
        drought_tolerance
        
        mature_year % number of years until reproductively mature 
        release_rate % annual release rate of seeds per plant
              
        status  %1 for alive, 0 for dead
        color % color for output
        % resprout_color
    end 
            
    methods     
                      
    function obj = plant(s,l,h,c,gf,a)
        %s = species, l=location, h=height_growth, c = crown_growth
        %gf = growth_factor, a = age
        obj.species=s.name;  %species name, e.g. Cs,Ml,Ro,Cm
        obj.resprouter=s.resprouter;   
        obj.location=l; 
        obj.mature_year = 6;
     
        % If we are specifying a new mature plant
        if(nargin>2)
          obj.age=a;
          obj.mature_age=a;
          if(a > 6)
              obj.adult = 1;
          end
          obj.height = h*gf;
          cr = c*gf;
          obj.crown_north = cr;
          obj.crown_south = cr;
          obj.crown_east = cr;
          obj.crown_west = cr;
          obj.crown_max = cr;
          obj.crown_avg = cr;
          obj.area = pi*cr*cr;
          obj.overlap_area = 0.0;
          % growth factors for adults
          obj.height_growth = h;
          obj.crown_growth = c;
           
        if obj.crown_growth > 1.0
            disp('adult crown_growth is large')
        end
          % number of seedlings per adult just prior to first fire
          obj.release_rate = s.longterm_release_rate;
        else
          obj.age = 0;
          obj.mature_age = 0;
          obj.adult = 0; % initially all plants are seedlings
          % fitted height and crown diameter using log y = a+b*x
          % so y0 = exp(a)
          obj.height = exp(s.s_avg_height + s.s_std_height*randn);
          % fits were done for the crown diam so we mult by 1/2
          sc = 0.5*exp(s.s_avg_crown + s.s_std_crown*randn);        
          obj.crown_north = sc;
          obj.crown_south = sc;
          obj.crown_east = sc;
          obj.crown_west = sc;
          obj.crown_max = sc;
          obj.crown_avg = sc;
          obj.area = pi*sc*sc; % area = pi*r^2
          obj.overlap_area = 0.0;
          % seedlings initially cannot produce seeds
          obj.release_rate = 0;
          % seedling growth parameters
          obj.height_growth = exp(s.avg_height_growth+s.std_height_growth*randn);
          obj.crown_growth = exp(s.avg_crown_growth+s.std_crown_growth*randn);
        end
        
        obj.collision=zeros(1,4);
        obj.competition_factor = s.seedling_competition_factor;
        obj.competition_distance = 1.5;
        obj.competition = 0.0;
        obj.drought_tolerance = s.drought_tolerance(1:2);
          
        obj.status = 1;  %plant starts as alive
        obj.color=s.color;
end

function [dist,alpha] = distance(obj,p)               
       vec1 = p.location(1)-obj.location(1);
       vec2 = p.location(2)-obj.location(2);   %vector between two plant center
       dist = realsqrt(vec1*vec1+vec2*vec2); % distance between two plants
       alpha = mod(atan2(vec2,vec1),2.0*pi); %angle of vector in radians
end
 
function compete(obj,p,d)
    % effect of plant p at a distance d away
    if(d < obj.competition_distance)
        obj.competition = obj.competition + obj.competition_factor*p.height*p.area/(obj.height*obj.area*d*d);
    end
end      
      
% Compute the overlap between two plants and how much they push back on
% each other
function impede = collide(obj,p,dist,alpha)
    % north, south, east, west impedance
    %impede(1) = 0.0;   %east
    %impede(2) = 0.0;   %north
    %impede(3) = 0.0;   %west
    %impede(4) = 0.0;   %south
    
    impede = zeros(1,4);
    
    crown_sum = obj.crown_max + p.crown_max;
    
    % crown overlap
    if (dist < crown_sum)    
        v=0.7; %impede percentage per plant
        quadrant = ceil(alpha*2.0/pi);
      
        calpha = cos(alpha);
        salpha = sin(alpha);
        
        switch quadrant
            case 1 % East & North
                b=obj.crown_north*salpha;
                a=obj.crown_east*calpha;
                d=p.crown_south*salpha;
                c=p.crown_west*calpha;
                
                overlap = dist - ( realsqrt( a*a + b*b ) + realsqrt( c*c + d*d ) );
                
                if(overlap<0)
                    scale1=p.crown_south/(p.crown_south+obj.crown_north);
                    scale2=p.crown_west/(p.crown_west+obj.crown_east);
                    impede(1)= calpha*v*scale1*overlap;
                    impede(2)= salpha*v*scale2*overlap;
                    r=max(p.crown_south,p.crown_west);
                    R=max(obj.crown_north,obj.crown_east);
                    a = (-dist+r+R)*(dist+r-R)*(dist-r+R)*(dist+r+R);
                    if(a>0)
                        obj.overlap_area = obj.overlap_area + 0.25*realsqrt(a);
                    end
                end
            case 2 % North & West
                b=obj.crown_north*salpha;
                a=obj.crown_west*calpha;
                d=p.crown_south*salpha;
                c=p.crown_east*calpha;
               
                overlap = dist - ( realsqrt( a*a + b*b ) + realsqrt( c*c + d*d ) );
                if(overlap<0)
                    scale1=p.crown_south/(p.crown_south+obj.crown_north);
                    scale2=p.crown_east/(p.crown_east+obj.crown_west);
                    impede(2)= salpha*v*scale2*overlap;
                    impede(3)= -calpha*v*scale1*overlap;
                    r=max(p.crown_south,p.crown_east);
                    R=max(obj.crown_north,obj.crown_west);
                    a = (-dist+r+R)*(dist+r-R)*(dist-r+R)*(dist+r+R);
                    if(a>0)
                        obj.overlap_area = obj.overlap_area + 0.25*realsqrt(a);
                    end
                end
            case 3 % West & South
                b=obj.crown_south*salpha;
                a=obj.crown_west*calpha;
                d=p.crown_north*salpha;
                c=p.crown_east*calpha;
                
                overlap = dist - ( realsqrt( a*a + b*b ) + realsqrt( c*c + d*d ) );
                
                if(overlap<0)
                    scale1=p.crown_north/(p.crown_north+obj.crown_south);
                    scale2=p.crown_east/(p.crown_east+obj.crown_west);
                    impede(3)= -calpha*v*scale1*overlap;
                    impede(4)= -salpha*v*scale2*overlap;
                    r=max(obj.crown_south,obj.crown_west);
                    R=max(p.crown_north,p.crown_east);
                    a = (-dist+r+R)*(dist+r-R)*(dist-r+R)*(dist+r+R);
                    if(a>0)
                        obj.overlap_area = obj.overlap_area + 0.25*realsqrt(a);
                    end
                end
            case 4 % East & South
                b=obj.crown_south*salpha;
                a=obj.crown_east*calpha;
                d=p.crown_north*salpha;
                c=p.crown_west*calpha;
                
                overlap = dist - ( realsqrt( a*a + b*b ) + realsqrt( c*c + d*d ) );
                
                if(overlap<0)
                    scale1=p.crown_north/(p.crown_north+obj.crown_south);
                    scale2=p.crown_west/(p.crown_west+obj.crown_east);
                    impede(1)= calpha*v*scale1*overlap ;
                    impede(4)= -salpha*v*scale2*overlap;
                    r=max(obj.crown_south,obj.crown_east);
                    R=max(p.crown_north,p.crown_west);
                    a = (-dist+r+R)*(dist+r-R)*(dist-r+R)*(dist+r+R);
                    if(a>0)
                        obj.overlap_area = obj.overlap_area + 0.25*realsqrt(a);
                    end
                end
        end
        obj.collision = min(obj.collision, impede);       
    end
end
 
function maturation(obj,species)
   % change to an adult plant
   obj.adult = 1;
   % update to growth parameters of adult plants
   if (obj.resprouter)
       obj.height_growth = obj.height_growth/16.795;
       obj.crown_growth = obj.crown_growth/15.24;
   else
       obj.height_growth = obj.height_growth/17.458;
       obj.crown_growth = obj.crown_growth/15.24;
   end
   
   % update germanation rate to adult rate
   obj.release_rate = species.effective_release_rate;
   % update competition rates
   obj.drought_tolerance = ones(1,2);
   obj.competition_factor = 0.0;
end


function grow(obj,rain)    
    obj.age = obj.age + 1; 
               
    %if plant is a seeding
    if (~obj.adult)
        h_growth = (obj.height_growth - 1.0)*obj.height;
        c_growth = (obj.crown_growth - 1.0)*obj.crown_max;
    % else plant is a adult
    else
        ra=rain/obj.age;
        h_growth = obj.height_growth*ra;
        c_growth = obj.crown_growth*ra;
        obj.mature_age = obj.mature_age+1;
    end
    
    % Grow crown sizes
    obj.crown_east = obj.crown_east + c_growth + obj.collision(1);
    obj.crown_north = obj.crown_north + c_growth + obj.collision(2);
    obj.crown_west = obj.crown_west + c_growth + obj.collision(3);
    obj.crown_south = obj.crown_south + c_growth + obj.collision(4);
    % Grow heights
    obj.height=obj.height+h_growth;
    % Update crown max and avg
    obj.crown_max = max([obj.crown_east,obj.crown_north, obj.crown_west, obj.crown_south]);
    obj.crown_avg=0.25*(obj.crown_east+obj.crown_north+obj.crown_west+obj.crown_south);
    % Update area
    obj.area = 0.25*pi*(obj.crown_north*obj.crown_east + obj.crown_north*obj.crown_west + obj.crown_south*obj.crown_east+obj.crown_south*obj.crown_west);

    obj.collision = zeros(1,4);
    obj.competition = 0.0;
 end    
                                     
 function resprout_plant(obj,species)
                
             if (obj.resprouter && obj.adult)
               obj.crown_north = 0.2;
               obj.crown_south = 0.2;
               obj.crown_east = 0.2;
               obj.crown_west = 0.2;
               obj.crown_max = 0.2;
               obj.height = 0.2;
               obj.age = 0;
               obj.competition = 0.0;
               obj.competition_factor = species.resprout_competition_factor;
               obj.competition_distance = 5.0;
               obj.drought_tolerance = species.drought_tolerance(3:4); 
               obj.release_rate = species.effective_release_rate;
             else
                obj.status = 0; %dead
                obj.height = 0; 
                obj.crown_east = 0;
                obj.crown_north = 0;
                obj.crown_west = 0;
                obj.crown_south = 0;
                obj.crown_max = 0;
                obj.age = 0;
             end
 end

 function dt = get_drought_tolerance(obj,rain)
    if (rain > 8.0)
        dt = obj.drought_tolerance(1);
    else
        dt = obj.drought_tolerance(2);
    end
 end
         
    end   
end
