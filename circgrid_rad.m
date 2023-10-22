function [bound,sen_pos_rad] = circgrid_rad(Q,R,radial_units,divs)

sen_pos_rad = zeros(Q,2);
for i=1:Q
    sen_pos_rad(i,:) = [radial_units i*360/Q];
end

angle_res = 1/divs;

bound = zeros(2,1/angle_res);
bound(2,:) = linspace(0,2*pi,1/angle_res);
bound(1,:) = R*(radial_units)/radial_units;

end

