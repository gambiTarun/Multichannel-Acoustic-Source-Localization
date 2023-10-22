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

noise_levs = 12:2:16;
source_radius = 0:0.5:8;
iterations = 100;

alg = 1;
Q = 80;
R = 0.15;
radial_units = 8;
T = 1000;
speed = 25;
wavelength_units = 4;

freq = speed/(wavelength_units*R/radial_units);
% freq = 1000;

% radius and angle of the source locations in [r(units),theta(deg)] form
% Source_locations = [1 180; 2 290; 3 220; 4 90; 5 18; 6 24; 7 45; 8 180];
% Source_locations = [6 180];


% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = (2*R+0.02)/128;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = speed;  % [m/s]
medium.alpha_coeff = 0;  % [dB/(MHz^y cm)]
medium.alpha_power = 0;

% create the time array
kgrid.makeTime(medium.sound_speed);
% t_end = 9;
% kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);
kgrid.Nt = T;

for radi = 1:length(source_radius)
      
    % define a the source points by setting 1 value at corresponding
    % coordinates of the p_mask
    source.p_mask = zeros(Nx, Ny);
    
    Source_locations = [source_radius(radi) 180];

    for loc=1:size(Source_locations,1)

        % introducing a source point in the p_mask grid
        r = Source_locations(loc,1)*R/radial_units;
        t = Source_locations(loc,2)*pi/180;
        [x,y] = pol2cart(t,r);
        xind = round((x/dx)+Nx/2);
        yind = round((y/dy)+Ny/2);
        source.p_mask(xind,yind) = 1;

        % define a time varying sinusoidal source
        source_freq = freq;   % [Hz]
        source_mag = 5;         % [Pa]
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
    %     s = awgn(s,10);
    %     s = source_mag*awgn(zeros(1,length(kgrid.t_array)),10);

        source.p(loc,:) = filterTimeSeries(kgrid, medium, s);
    end

    % define a centered circular sensor
    sensor_radius = R;   % [m]
    num_sensor_points = 80;
    sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

    % define the acoustic parameters to record
    sensor.record = {'p', 'p_final'};

    % run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

    % the sensor location is stored in sensor.mask in their x and y
    % coordinates. These are converted to the polar form.
    [theta,~] = cart2pol(sensor.mask(1,:),sensor.mask(2,:));

    % actual location of the source 
    [x_loc,y_loc] = find(source.p_mask);
    x_loc = (x_loc-Nx/2)*dx;
    y_loc = (y_loc-Ny/2)*dy;


    % wavenumber k and the maximum mode M
    k = 2*pi*source_freq/medium.sound_speed;
    m = ceil(11/9*k*R);

    ind1 = repmat(-m:m,num_sensor_points,1);
    ind2 = repmat(theta',1,2*m+1);

%     dis = zeros(1,length(noise_levs));
    
    disp("src_radius: "+source_radius(radi)+" units");
    
    for n=1:length(noise_levs)

        distance = zeros(1,iterations);

        disp("snr: "+noise_levs(n)+" db");

        for iter=1:iterations

            % the sensor readings
            z = sensor_data.p;
            
%             save("Zdata/Z_rad"+source_radius(radi)+"u_freq"+freq+"Hz.mat",'z');
            
            z = awgn(z,noise_levs(n));

            if(alg==1)

                xi = exp(1i*ind1.*ind2);

                % Taking the inverse of H(sensor radial component) if the condition number
                % is small enough. Took a threshold randomly
                if(cond(diag(besselh(-m:m,1,k*R))) < 100)
                    Hinv = diag(besselh(-m:m,1,k*R).^-1);
                else
                    Hinv = 1;
                end

                % the spatial fourier coefficients
                alpha = Hinv * 4/(num_sensor_points*1i)*ctranspose(xi)*z;

            elseif(alg==2)

                gamma = besselh(ind1,1,k*R).*exp(1i*ind1.*ind2);

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

            [A,B,C] = pol2cart(P_grid,Y_grid,Z);

            % predicted location of the source
            [xm,ym] = find(abs(C) == max(max(abs(C))));

            xmax = A(xm(1),ym(1));
            ymax = B(xm(1),ym(1));

            distance(iter) = sqrt((xmax-x_loc)^2 + (ymax-y_loc)^2);

            disp("iter: "+iter);
            
            save("accVnoiseVradius/MV_iter"+iter+"_rad"+source_radius(radi)+"u_snr"+noise_levs(n)+"db.mat",'C');

        end

        save("accVnoiseVradius/dist_rad"+source_radius(radi)+"u_snr"+noise_levs(n)+"db.mat",'distance');

    end
    
%     plot_mat = [noise_levs; dis];
    
%     save("accVnoiseVradius/plot_rad"+source_radius(radi)+".mat",'plot_mat');
    
%     figure(1);
%     fig = imagesc(kgrid.x_vec, kgrid.y_vec, flip((sensor_data.p_final + source.p_mask + cart2grid(kgrid, sensor.mask))',1), [-1, 1]);
%     colormap(getColorMap);
%     ylabel('y-position [m]');
%     xlabel('x-position [m]');
%     grid on;
%     title(source_radius(radi)+"units");
%     axis image;
%     saveas(fig,"accVnoiseVradius/fig_"+source_radius(radi)+"_units.fig");
% 
%     figure(2);
%     fig = plot(noise_levs,dis);
%     xlabel('noise(db)');
%     ylabel('dis');
%     title(source_radius(radi)+"units");
%     saveas(fig,"accVnoiseVradius/plot_"+source_radius(radi)+"_units.fig");
end
    