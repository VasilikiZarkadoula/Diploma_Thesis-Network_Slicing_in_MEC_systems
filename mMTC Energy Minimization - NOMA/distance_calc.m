function [dist,Xcoord,Ycoord] = distance_calc(d_min_users, d_min_bs, radius , users)

rng shuffle;

d_min_users = d_min_users^2;
d_min_bs = d_min_bs^2;

Xcoord = zeros(users, 1);
Ycoord = zeros(users, 1);
R = zeros(users, 1);

xc = 0;           
yc = 0;  

uValid=0;
iteration = 1; 

while uValid < users && iteration < 1e8
    theta = 2*pi*rand;
    r = sqrt(rand);
    XcoordNew = (radius*r)*cos(theta) + xc;
    YcoordNew = (radius*r)*sin(theta) + yc;
    
    if all(((Xcoord(1:uValid) - XcoordNew).^2 + (Ycoord(1:uValid) - YcoordNew).^2) > d_min_users)
        if (XcoordNew^2+YcoordNew^2)>d_min_bs
            uValid    = uValid + 1; 
            Xcoord(uValid) = XcoordNew;
            Ycoord(uValid) = YcoordNew;
            R(uValid) = sqrt((XcoordNew^2 + YcoordNew^2));
        end
    end
  iteration = iteration + 1;
end

if uValid < users
  error('D istance_calc cannot create wanted number of points in %d iterations.', iteration)
end
    
dist = R;
end