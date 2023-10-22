function [bound,sen_pos_rad] = rectgrid_rad(Q,R,lr,br,radial_units,divs)

length = lr/R*radial_units;
breadth = br/R*radial_units;
l = lr;
b = br;

sen_pos_rad = zeros(Q,2);
for i=1:Q
    dis = i*2*(length+breadth)/Q;
    if(dis <= breadth/2)
        sen_pos_rad(i,:) = [sqrt((length/2)^2 + dis^2) atan((dis)/(length/2))*180/pi];
    elseif(dis <= (length+breadth)/2)
        sen_pos_rad(i,:) = [sqrt((breadth/2)^2 + ((length+breadth)/2-dis)^2) 90-atan(((length+breadth)/2-dis)/(breadth/2))*180/pi];
    elseif(dis <= length+breadth/2)
        sen_pos_rad(i,:) = [sqrt((breadth/2)^2 + (dis-(length+breadth)/2)^2) 90+atan((dis-(length+breadth)/2)/(breadth/2))*180/pi];
    elseif(dis <= length+breadth)
        sen_pos_rad(i,:) = [sqrt((length/2)^2 + ((length+breadth)-dis)^2) 180-atan(((length+breadth)-dis)/(length/2))*180/pi];
    elseif(dis <= length+breadth*3/2)
        sen_pos_rad(i,:) = [sqrt((length/2)^2 + (dis-(length+breadth))^2) 180+atan((dis-(length+breadth))/(length/2))*180/pi];
    elseif(dis <= (length+breadth)*3/2)
        sen_pos_rad(i,:) = [sqrt((breadth/2)^2 + ((length+breadth)*3/2-dis)^2) 270-atan(((length+breadth)*3/2-dis)/(breadth/2))*180/pi];
    elseif(dis <= length*2+breadth*3/2)
        sen_pos_rad(i,:) = [sqrt((breadth/2)^2 + (dis-(length+breadth)*3/2)^2) 270+atan((dis-(length+breadth)*3/2)/(breadth/2))*180/pi];
    elseif(dis <= 2*(length+breadth))
        sen_pos_rad(i,:) = [sqrt((length/2)^2 + (2*(length+breadth)-dis)^2) 360-atan((2*(length+breadth)-dis)/(length/2))*180/pi];
    end  

end

angle_res = 1/divs;

bound = zeros(2,1/angle_res);
bound(2,:) = linspace(0,2*pi,1/angle_res);
for i=1:size(bound,2)
    if(bound(2,i)<=atan(b/l) || bound(2,i)>2*pi-atan(b/l))
        bound(1,i) = l*0.5/cos(bound(2,i));
    elseif(atan(b/l)<bound(2,i) && bound(2,i)<=pi-atan(b/l))
        bound(1,i) = b*0.5/cos(pi/2-bound(2,i));
    elseif(pi-atan(b/l)<bound(2,i) && bound(2,i)<=pi+atan(b/l))
        bound(1,i) = l*0.5/cos(pi-bound(2,i));
    elseif(pi+atan(b/l)<bound(2,i) && bound(2,i)<=2*pi-atan(b/l))
        bound(1,i) = b*0.5/cos(3*pi/2-bound(2,i));
    end
end

end

