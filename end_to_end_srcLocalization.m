clear
close all
clc

clearvars;

% =========================================================================
% APPARATUS
% =========================================================================
samples = 20;
source_mag = [5 25 50];

alg = 1;
Q = 40;
R = 0.15;
radial_units = 8;
T = 1000;
speed = 25;
wavelength_units = 2;

freq = speed/(wavelength_units*R/radial_units);
% freq = 1000;

Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 2*(R/radial_units*9)/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]

src_locations = zeros(2,samples);
for rand_locs = 1:samples
    
    % radius and angle of the source locations in [r(units),theta(deg)] form
    % src_pos_rad = [1 180; 2 290; 3 220; 4 90; 5 18; 6 24; 7 45; 8 180];
    src_pos_rad = [rand()*8 rand()*360];

    % [radial_bound,sen_pos_rad] = rectgrid_rad(Q,R,2*R*30/30,2*R*15/30,radial_units,360);
    [radial_bound,sen_pos_rad] = circgrid_rad(Q,R,radial_units,360);
    % [radial_bound,sen_pos_rad] = semicircgrid_rad(Q,R,radial_units,360);
    
    for sm=1:length(source_mag)
        % =========================================================================
        % SIMULATION
        % =========================================================================

        % create the computational grid
        kgrid = kWaveGrid(Nx, dx, Ny, dy);

        % define the properties of the propagation medium
        medium.sound_speed = speed;  % [m/s]
        % medium.alpha_coeff = 0;
        % medium.alpha_power = 0;
        medium.alpha_mode = 'no_dispersion';

        % create the time array
        kgrid.makeTime(medium.sound_speed);
        % t_end = 9;
        % kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);
        record_start = round(2*R*sqrt(2)/(speed*kgrid.dt));
        kgrid.Nt = T + record_start;

        % define a the source points by setting 1 value at corresponding
        % coordinates of the p_mask
        source.p_mask = zeros(Nx, Ny);
        source.p = [];

        for loc=1:size(src_pos_rad,1)

            % introducing a source point in the p_mask grid
            r = src_pos_rad(loc,1)*R/radial_units;
            t = src_pos_rad(loc,2)*pi/180;
            [x,y] = pol2cart(t,r);
            xind = round((x/dx)+Nx/2);
            yind = round((y/dy)+Ny/2);
            source.p_mask(xind,yind) = 1;

            % define a time varying sinusoidal source
            source_freq = freq;
            s = source_mag(sm) * sin(2 * pi * source_freq * kgrid.t_array);

            if isempty(source.p)
                source.p(1,:) = filterTimeSeries(kgrid, medium, s);
                continue;
            end
            source.p = [filterTimeSeries(kgrid, medium, s); source.p];
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
        % sensor.mask = makeCartCircle(R, Q);
        sensor.mask = sen_pos_cart;

        % define the acoustic parameters to record
        sensor.record = {'p', 'p_final'};

        % run the simulation
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

        % the sensor location is stored in sensor.mask in their x and y
        % coordinates. These are converted to the polar form.
        [theta,radius] = cart2pol(sensor.mask(1,:),sensor.mask(2,:));

        % actual location of the source
        [x_loc,y_loc] = find(source.p_mask);
        x_loc = (x_loc-Nx/2)*dx;
        y_loc = (y_loc-Ny/2)*dy;

        % wavenumber k and the maximum mode M
        k = 2*pi*source_freq/medium.sound_speed;
        m = ceil(11/9*k*R);

        % the sensor readings
        z = sensor_data.p(:,record_start+1:end);
        z = awgn(z,10);

        save("ZdataGen/"+rand_locs+"_mag"+source_mag(sm)+".mat",'z');

        src_locations(:,rand_locs) = src_pos_rad;
    end
end

save("ZdataGen/loc.mat",'src_locations');
return;

% The matrix for fourier transformation
ind1 = repmat(-m:m,Q,1);    % varying M
ind21 = repmat(theta',1,2*m+1);     % varying theta
ind22 = repmat(radius',1,2*m+1);    % varying radii

if(alg==1)

    xi = exp(1i*ind1.*ind21);

    Hinv = diag(besselh(-m:m,1,k*R).^-1);
    
    % the spatial fourier coefficients
    alpha = Hinv * 4/(Q*1i)*ctranspose(xi)*z;

elseif(alg==2)
    
    gamma = besselh(ind1,1,k*ind22).*exp(1i*ind1.*ind21);

    % computing alpha
    alpha = 4/(1i)*pinv(gamma)*z;
    
end
    
% the covariance of the fourier coefficients
Ra = 1/T*alpha*ctranspose(alpha);

% the polar coordinates vectors, y and phi
Y = linspace(0, R, 200);
P = linspace(0, 2*pi, 180);

[P_grid,Y_grid] = meshgrid(P,Y);

% initialising the MV spectrum matrix 
Z = zeros(length(Y),length(P));

RaI = pinv(Ra);

for i=1:length(Y)
    for j=1:length(P)
        c = besselj(-m:m,k*Y(i)).*exp(1i*(-m:m)*P(j));
        Z(i,j) = 10*log(1/real(c*RaI*ctranspose(c)));
    end
end

[A,B,C] = pol2cart(P_grid,Y_grid,Z);

% predicted location of the source
pks = find(imregionalmax(C));
[val_sorted, ind] = sort(C(pks), 'descend');
valk = val_sorted(1:sum(sum(source.p_mask)));
pks_sorted = pks(ind);
pksk = pks_sorted(1:sum(sum(source.p_mask)));
% xmax = xmax(1);
% ymax = ymax(1);

figure(6);
mesh(A,B,C);
ylabel('x-position [m]');
xlabel('y-position [m]');
zlabel('Z');

fig = figure;
mesh(A,B,C);
ylabel('y-position [m]');
xlabel('x-position [m]');
zlabel('Z(db)');
hold on;
axis tight;
scatter3(x_loc,y_loc,ones(length(x_loc),1)*valk(1),150,'X','r');
scatter3(A(pks),B(pks),C(pks),150,'X','y');
scatter3(A(pksk),B(pksk),valk,150,'X','k');
legend('MV Spec', 'Source','All peaks', 'Top peaks');
view(2);
