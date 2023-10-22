function [bound,sen_pos_rad] = ellipsogrid_rad(Q,R,a,b,radial_units,bound_divs)

a = a*radial_units;
b = b*radial_units;
sen_pos_rad = zeros(Q,2);
for i=1:Q
    sen_pos_rad(i,:) = [sqrt((a*b)/(b*cosd(i*360/Q)^2+a*sind(i*360/Q)^2)) i*360/Q];
end

angle_res = 1/bound_divs;

bound = zeros(2,1/angle_res);
bound(2,:) = linspace(0,2*pi,1/angle_res);
bound(1,:) = R*sqrt((a*b)./(b*cos(bound(2,:)).^2+a*sin(bound(2,:)).^2))/radial_units;

end

