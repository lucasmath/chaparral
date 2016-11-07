function ellipsefill(e, n, w, s , x0, y0, color, Edge)
%A function that computes 400 x and y points to make and fill four
%elliptial quadrants with radii e,n,w and s centered at (x0,y0).

 ta = linspace(0,0.5*pi);
 tb = linspace(0.5*pi,pi);
 tc = linspace(pi,1.5*pi);
 td = linspace(1.5*pi, 2*pi);
 total = length(ta) + length(tb) + length(tc) + length(td);
 x = zeros(1, total);
 y = zeros(1, total);
 
 for i=1:total;
     if(i <= 0.25*total)
         x(i) = x0 + e*cos(ta(i));
         y(i) = y0 + n*sin(ta(i));
     elseif(i <= 0.5*total);
         x(i) = x0 + w*cos(tb(i-length(ta)));
         y(i) = y0 + n*sin(tb(i-length(ta)));
     elseif(i <= 0.75*total);
         x(i) = x0 + w*cos(tc(i-length(tb)-length(ta)));
         y(i) = y0 + s*sin(tc(i-length(tb)-length(ta))); 
     elseif(i <= total)
         x(i) = x0 + e*cos(td(i-length(tc) - length(tb) - length(ta)));
         y(i) = y0 + s*sin(td(i-length(tc) - length(tb) - length(ta)));
     end
 end
 
 patch(x,y, color,'EdgeColor',Edge,'LineWidth',2)

end