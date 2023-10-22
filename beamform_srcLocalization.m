clear
close all
clc

% =========================================================================
% APPARATUS
% =========================================================================

beamformer = 1;
Q = 40;
R = 0.15;
radial_units = 6;
canvas_units = 9;
speed = 25;
T = 1000;
wavelength_units = 4;

freq = speed/(wavelength_units*R/radial_units);
% freq = 15e3;

Nx = 129;           % number of grid points in the x (row) direction
Ny = 129;           % number of grid points in the y (column) direction
dx = 2*(R/radial_units*canvas_units)/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]

% =========================================================================
% SOURCES
% =========================================================================

% radius and angle of the source locations in [r(units),theta(deg)] form
% src_pos_rad = [1 180; 2 290; 3 220; 4 90; 5 18; 6 24; 7 45; 8 180];
src_pos_rad = [1 180; 3 340; 4 90;];
    
airway_src = rand_lung_cart(0,0,R,radial_units,Nx,Ny,dx,dy);
% airway_src(1) = airway_src(1)-15;

num = 0;
seed = 21;
rng(seed);
locX = randi([round(Nx/2-Nx/2*(7)/9),round(Nx/2-Nx/2*(1)/9)],num,1);
rng(seed+1);
locY = randi([round(Ny/2+Nx/2*(-3)/9), round(Ny/2+Nx/2*(3)/9)],num,1);
fault_srcl = [locX locY];

% =========================================================================
% SENSORS
% =========================================================================

% [radial_bound,sen_pos_rad] = rectgrid_rad(Q,R,2*R*20/30,2*R*30/30,radial_units,360);
[radial_bound,sen_pos_rad] = circgrid_rad(Q,R,radial_units,360);
% [radial_bound,sen_pos_rad] = ellipsogrid_rad(Q,R*0.6,2,9,radial_units*0.6,360);
% [radial_bound,sen_pos_rad] = semicircgrid_rad(Q,R,radial_units,360);

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = speed ;   % [m/s]
medium.alpha_coeff = 0;
medium.alpha_power = 1.5;
% medium.sound_speed = 25 * ones(Nx, Ny);   % [m/s]
% medium.sound_speed(1:Nx/2, :) = 50;       % [m/s]
% medium.density = 1000 * ones(Nx, Ny);       % [kg/m^3]
% medium.density(:, Ny/4:Ny) = 1200;          % [kg/m^3]

% create the time array
kgrid.makeTime(medium.sound_speed);
% t_end = 9;
% kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);
record_start = round(2*R*sqrt(2)/(speed*kgrid.dt));
kgrid.Nt = T + record_start;

% define a the source points by setting 1 value at corresponding
% coordinates of the p_mask
source.p_mask = zeros(Nx, Ny);
type_mat = zeros(Nx,Ny);
source.p = [];

for loc=1:size(src_pos_rad,1)
    
    % introducing a source point in the p_mask grid
    r = src_pos_rad(loc,1)*R/radial_units;
    t = src_pos_rad(loc,2)*pi/180;
    [x,y] = pol2cart(t,r);
    xind = round((x/dx)+Nx/2);
    yind = round((y/dy)+Ny/2);
    source.p_mask(xind,yind) = 1;
    type_mat(xind,yind) = 1;

end

for loc=1:size(airway_src,1)

    x = airway_src(loc,1);
    y = airway_src(loc,2);
    source.p_mask(x,y) = 1;
    type_mat(x,y) = 2;
    
end

for loc=1:size(fault_srcl,1)
    
    % introducing a source point in the p_mask grid
    x = fault_srcl(loc,1);
    y = fault_srcl(loc,2);
    source.p_mask(x,y) = 1;
    type_mat(x,y) = 3;

end

[srcX,srcY,val] = find(type_mat);

for i = 1:length(srcX)
    
    if val(i)==1
        
        % define a time varying sinusoidal source
        source_freq = freq;   % [Hz]
        source_mag = 20;
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

    elseif val(i)==2
        
        % define a time varying sinusoidal source
        source_freq = freq;   % [Hz]
        source_mag = 50;
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

    elseif val(i)==3
        
        % define a time varying sinusoidal source
        source_freq = freq;   % [Hz]
        source_mag = 2;
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

    end
    
    if isempty(source.p)
        source.p(1,:) = filterTimeSeries(kgrid, medium, s);
    else
        source.p = [source.p; filterTimeSeries(kgrid, medium, s)];
    end
    
end

% define a centered circular sensor
sen_pos_cart = zeros(2,Q);
for loc=1:Q  
    % introducing a source point in the p_mask grid
    r = sen_pos_rad(loc,1)*R/radial_units;
    t = sen_pos_rad(loc,2)*pi/180;
    [x,y] = pol2cart(t,r);   
    sen_pos_cart(:,loc) = [x y];

end    
% sensor.mask = makeCartCircleplt.scatter(yfr_test[:len_*96,0],yfr_test[:len_*96,1],c='k')(R, Q);
sensor.mask = sen_pos_cart;

% define the acoustic parameters to record
sensor.record = {'p', 'p_final'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the final wave-field
figure;
imagesc(kgrid.x_vec, kgrid.y_vec, flip((sensor_data.p_final + source.p_mask + cart2grid(kgrid, sensor.mask))',1), [-1, 1]);
colormap(getColorMap);
ylabel('y-position [m]');
xlabel('x-position [m]');
grid on;
axis image;

% plot the simulated sensor data
figure;
subplot(2,1,1);
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
plot(kgrid.t_array, source.p, 'k-');
xlabel('Time [s]');
ylabel('Signal Amplitude');
axis tight;
title('Input Pressure Signal');

subplot(2,1,2);
imagesc(sensor_data.p, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;

% the sensor location is stored in sensor.mask in their x and y
% coordinates. These are converted to the polar form.
[theta,radius] = cart2pol(sensor.mask(1,:),sensor.mask(2,:));

% actual location of the source
[sx_loc,sy_loc] = find(type_mat==1);
[hx_loc,hy_loc] = find(type_mat==2);
[fx_loc,fy_loc] = find(type_mat==3);
sx_loc = (sx_loc-Nx/2)*dx;
sy_loc = (sy_loc-Ny/2)*dy;
hx_loc = (hx_loc-Nx/2)*dx;
hy_loc = (hy_loc-Ny/2)*dy;
fx_loc = (fx_loc-Nx/2)*dx;
fy_loc = (fy_loc-Ny/2)*dy;

% wavenumber k and the maximum mode M
k = 2*pi*source_freq/medium.sound_speed;
m = ceil(11/9*k*R);

% the sensor readings
z = awgn(sensor_data.p(:,record_start:end),40);
figure;
imagesc(z, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;

Ns = 760;
fs = 1/kgrid.dt;
windows = size(z,2)-Ns+1;
zf = zeros(Q,windows);
for senum=1:Q
    for shift=1:windows
        yf = 2/Ns*(fft(z(senum,shift:Ns+shift-1)));
        zf(senum,shift) = max(yf(1:Ns/2));
    end
end    
% xf = linspace(0,fs/2,Ns/2);
% plot(xf,abs(yf(1:Ns/2)));

z_ = zf;

% The matrix for fourier transformation
ind1 = repmat(-m:m,Q,1);    % varying M
ind21 = repmat(theta',1,2*m+1);     % varying theta
ind22 = repmat(radius',1,2*m+1);    % varying radii

% if(alg==1)
%     
%     disp('change alg');
%     
%     xi = exp(1i*ind1.*ind21);
% 
%     Hinv = diag(besselh(-m:m,1,k*R).^-1);
%     
%     % the spatial fourier coefficients
%     alpha = Hinv * 4/(Q*1i)*ctranspose(xi)*z;
% 
% elseif(alg==2)
%     
% end
    
gamma = besselh(ind1,1,k*ind22).*exp(1i*ind1.*ind21);

% computing alpha
alpha = 4/(1i)*pinv(gamma)*z_;

% the covariance of the fourier coefficients
Ra = 1/size(alpha,2)*alpha*ctranspose(alpha);

% the polar coordinates vectors, y and phi
NNx = 129;
NNy = 129;
ddx = 2*(R/radial_units*canvas_units)/NNx;
ddy = ddx;

X_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,NNx);
Y_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,NNy);

[X_grid,Y_grid] = meshgrid(X_var,Y_var);

[P,Y] = cart2pol(X_grid,Y_grid);
P = mod(P,2*pi);

% initialising the MV spectrum matrix 
Z = zeros(NNx,NNy);

RaI = pinv(Ra);

for i=1:NNx
    for j=1:NNy
        [~,min_index] = min(abs(radial_bound(2,:)-P(i,j)));
        if(Y(i,j)<radial_bound(1,min_index))
            c = besselj(-m:m,k*Y(i,j)).*exp(1i*(-m:m)*P(i,j));
            if(beamformer==1)   % MV beamformer
                Z(i,j) = (c*RaI*ctranspose(c))^-1;
            elseif(beamformer==2)   % Bartlett beamformer
                Z(i,j) = (c*Ra*ctranspose(c));
            end
        end
    end
end

[A,B,C] = pol2cart(P,Y,real(Z));
[x_bound,y_bound] = pol2cart(radial_bound(2,:),radial_bound(1,:));

% predicted location of the source
pks = find(imregionalmax(C));
[val_sorted, ind] = sort(C(pks), 'descend');
valk = val_sorted(1:min(length(pks),sum(sum(source.p_mask))));
pks_sorted = pks(ind);
pksk = pks_sorted(1:min(length(pks),sum(sum(source.p_mask))));
% xmax = xmax(1);
% ymax = ymax(1);

figure;
mesh(A,B,C);
ylabel('x-position [m]');
xlabel('y-position [m]');
zlabel('Z');

fig = figure;
mesh(A,B,C);
ylabel('y-position [m]');
xlabel('x-position [m]');
hold on;
scatter3(sx_loc,sy_loc,ones(length(sx_loc),1)*valk(1),150,'X','b');
% scatter3(A(pks),B(pks),C(pks),150,'X','y');
scatter3(A(pksk),B(pksk),valk,150,'X','k');
scatter3(x_bound,y_bound,zeros(length(x_bound),1),100,'o','k');
legend('MV Spec', 'Source', 'Top peaks', 'Sensors');
view(2);
