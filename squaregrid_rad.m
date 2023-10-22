function [bound,sen_pos_rad] = squaregrid_rad(Q,R,radial_units,divs)

sen_pos_rad = zeros(Q,2);
for i=1:Q
    sen_pos_rad(i,:) = [(radial_units)/max(abs(cos(i*2/Q*pi)),abs(sin(i*2/Q*pi))) i*360/Q];
end

angle_res = 1/divs;

bound = zeros(2,1/angle_res);
bound(2,:) = linspace(0,2*pi,1/angle_res);
size(bound,2)
for i=1:size(bound,2)
    bound(1,i) = R/max(abs(cos(bound(2,i))),abs(sin(bound(2,i))))*(radial_units)/radial_units;
end

end

