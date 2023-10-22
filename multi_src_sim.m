clear
close all
clc

% Monopole Point Source In A Homogeneous Propagation Medium Example
%
% This example provides a simple demonstration of using k-Wave for the
% simulation and detection of a time varying pressure source within a
% two-dimensional homogeneous propagation medium.  It builds on the
% Homogeneous Propagation Medium and Recording The Particle Velocity
% examples.
%
% author: Bradley Treeby
% date: 2nd December 2009
% last update: 4th May 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% =========================================================================
% APPARATUS
% =========================================================================
alg = 2;
Q = 40;
R = 0.15;
radial_units = 8;
T = 1000;
speed = 25;
wavelength_units = 4;

freq = speed/(wavelength_units*R/radial_units);
% freq = 1000;

Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 2*(R/radial_units*9)/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]

% radius and angle of the source locations in [r(units),theta(deg)] form
% src_pos_rad = [1 180; 2 290; 3 220; 4 90; 5 18; 6 24; 7 45; 8 180];
src_pos_rad = [round(Nx/2+Nx/2*(3)/9) round(Nx/2+Nx/2*(2)/9);
                round(Nx/2+Nx/2*(7)/9) round(Nx/2+Nx/2*(2)/9);
                round(Nx/2+Nx/2*(5)/9) round(Nx/2+Nx/2*(-4)/9);
                round(Nx/2+Nx/2*(3)/9) round(Nx/2+Nx/2*(-2)/9)];

sen_pos_rad = zeros(Q,2);
for i=1:Q
    sen_pos_rad(i,:) = [radial_units i*360/Q];
%     sen_pos_rad(i,:) = [radial_units/max(abs(cos(i*2/Q*pi)),abs(sin(i*2/Q*pi))) i*360/Q];
end

% src_pos_rad = [src_pos_rad; healty_srcl];

% =========================================================================
% SIMULATION
% =========================================================================

for skipper = 1:5
    
    % create the computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    % define the properties of the propagation medium
    medium.sound_speed = 25;  % [m/s]
    % medium.alpha_coeff = 0;
    % medium.alpha_power = 0;
    % medium.alpha_mode = 'no_dispersion';

    % create the time array
    kgrid.makeTime(medium.sound_speed);
    % t_end = 9;
    % kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);
    kgrid.Nt = T;

    % define a the source points by setting 1 value at corresponding
    % coordinates of the p_mask
    source.p_mask = zeros(Nx, Ny);
    source.p = [];

    for loc=1:size(src_pos_rad,1)

        if(loc==skipper)
            continue
        end

        % introducing a source point in the p_mask grid
        x = src_pos_rad(loc,1);
        y = src_pos_rad(loc,2);
        source.p_mask(x,y) = 1;

        % define a time varying sinusoidal source
        source_freq = 333;   % [Hz]
        source_mag = 5;
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

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
    [t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));

    plot(kgrid.t_array, source.p, 'k-');
    xlabel('Time [s]');
    ylabel('Signal Amplitude');
    axis tight;
    title('Input Pressure Signal');

    figure;
    imagesc(sensor_data.p, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    colorbar;

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
    z = sensor_data.p;
    z = awgn(z,10);

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
            Z(i,j) = (c*RaI*ctranspose(c))^-1;
        end
    end

    [A,B,C] = pol2cart(P_grid,Y_grid,abs(Z));

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
    hold on;
    scatter3(x_loc,y_loc,ones(length(x_loc),1)*valk(1),150,'X','r');
    scatter3(A(pks),B(pks),C(pks),150,'X','g');
    scatter3(A(pksk),B(pksk),valk,150,'X','k');
    view(2);
    saveas(fig,char("multi_src_sim/MVspec"+skipper+".fig"));
    saveas(fig,char("multi_src_sim/MVspec"+skipper+".jpeg"));
end