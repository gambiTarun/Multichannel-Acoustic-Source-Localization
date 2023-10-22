function [bound,sen_pos_rad] = semicircgrid_rad(Q,R,radial_units,divs)

sen_pos_rad = zeros(Q,2);
for i=1:Q

%     Center of semi circle considered at 4r/(3pi)
        
    dis = i*(pi+2)*radial_units/Q;
    
    sen_pos_rad(i,:) = [radial_units dis/(pi*radial_units)*180];   
    if(dis/(pi*radial_units)*180 > 180)
        edis = (dis/(pi*radial_units)*180 - 180)*pi*radial_units/180;
        sen_pos_rad(i,:) = [radial_units-edis 180];
    end
    
end

rshift = [4*radial_units/(3*pi) 270];
radial_shift = [4*R/(3*pi) 270*pi/180];

for i=1:Q
    
    r = sen_pos_rad(i,1);
    t = sen_pos_rad(i,2);
    
    sen_pos_rad(i,1) = sqrt((r*sind(t)+rshift(1)*sind(rshift(2)))^2 + (r*cosd(t)+rshift(1)*cosd(rshift(2)))^2);
    sen_pos_rad(i,2) = mod(atan2((r*sind(t)+rshift(1)*sind(rshift(2))),(r*cosd(t)+rshift(1)*cosd(rshift(2)))),2*pi)*180/pi;
    
end

angle_res = 1/divs;

bound = zeros(2,1/angle_res);
theta_orginal = linspace(0,2*pi,1/angle_res);
for i=1:length(theta_orginal)
    
    t = theta_orginal(i);
    bound(2,i) = mod(atan2((R*sin(t)+radial_shift(1)*sin(radial_shift(2))),(R*cos(t)+radial_shift(1)*cos(radial_shift(2)))),2*pi);
    if(t<=pi)
        bound(1,i) = sqrt((R*sin(t)+radial_shift(1)*sin(radial_shift(2)))^2 + (R*cos(t)+radial_shift(1)*cos(radial_shift(2)))^2);
    else
        bound(1,i) = radial_shift(1)/cos(pi/2-(bound(2,i)-pi)); 
    end
end

end

