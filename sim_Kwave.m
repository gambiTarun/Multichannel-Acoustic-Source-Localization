clear
close all
clc

% =========================================================================
% APPARATUS
% =========================================================================

num_pos = 10;
noise_levs = 40:-5:0;
iterations = 100;
src_num = 2;
db_size = 60;

beamformer = 1;
Q = 40;
R = 0.15;
radial_units = 6;
canvas_units = 9;
speed = 25;
wavelength_units = 4;
peaks_recorded = 10;

freq = speed/(wavelength_units*R/radial_units);
% freq = 1000;

Nx = 129;           % number of grid points in the x (row) direction
Ny = 129;           % number of grid points in the y (column) direction
dx = 2*(R/radial_units*canvas_units)/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]

T = peaks_recorded*wavelength_units*R/(radial_units*speed);

X = -radial_units/sqrt(2):sqrt(2)*radial_units/(sqrt(db_size+4)-1):radial_units/sqrt(2);
Y = -radial_units/sqrt(2):sqrt(2)*radial_units/(sqrt(db_size+4)-1):radial_units/sqrt(2);

[X_grid, Y_grid] = meshgrid(X,Y);

[t,r] = cart2pol(reshape(X_grid,[],1), reshape(Y_grid,[],1));

sim_loc_db_rad = [r t*180/pi];

% deleting corners
for l=size(sim_loc_db_rad):-1:1
    if(abs(sim_loc_db_rad(l,1)-radial_units)<=1e-7)
        sim_loc_db_rad([l],:) = [];
    end
end

for pos=1:num_pos
    
    % =========================================================================
    % SOURCES
    % =========================================================================

    % radius and angle of the source locations in [r(units),theta(deg)] form
%     sim_loc_rad = [3 220; 2 90];
    sim_loc_rad = sim_loc_db_rad(randi(50,2,1),:);

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
    medium.sound_speed = 25 ;   % [m/s]
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
    snaps = round(T/kgrid.dt);
    kgrid.Nt = snaps + record_start;

    % define a the source points by setting 1 value at corresponding
    % coordinates of the p_mask
    source.p_mask = zeros(Nx, Ny);
    source.p = [];

    for loc=1:size(sim_loc_rad,1)

        % introducing a source point in the p_mask grid
        r = sim_loc_rad(loc,1)*R/radial_units;
        t = sim_loc_rad(loc,2)*pi/180;
        [x,y] = pol2cart(t,r);
        xind = round((x/dx)+Nx/2);
        yind = round((y/dy)+Ny/2);
        source.p_mask(xind,yind) = 1;

        source_freq = freq;   % [Hz]
        source_mag = 20;
        s = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

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

    % the sensor location is stored in sensor.mask in their x and y
    % coordinates. These are converted to the polar form.
    [theta,radius] = cart2pol(sensor.mask(1,:),sensor.mask(2,:));

    % actual location of the source
    [sx_loc,sy_loc] = find(source.p_mask);

    sx_loc = (sx_loc-Nx/2)*dx;
    sy_loc = (sy_loc-Ny/2)*dy;
  
    % wavenumber k and the maximum mode M
    k = 2*pi*source_freq/medium.sound_speed;
    m = ceil(11/9*k*R);

    % the polar coordinates vectors, y and phi
    X_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,Nx);
    Y_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,Ny);

    [X_grid,Y_grid] = meshgrid(X_var,Y_var);

    [P,Y] = cart2pol(X_grid,Y_grid);
    P = mod(P,2*pi);
    
%     src_pos_rad = [];
%     for i=1:Nx
%         for j=1:Ny
%             [~,min_index] = min(abs(radial_bound(2,:)-P(i,j)));
%             if(Y(i,j)<radial_bound(1,min_index))
%                 src_pos_rad = [src_pos_rad;Y(i,j)*radial_units/R P(i,j)*180/pi];
%             end
%         end
%     end

    for noise=noise_levs
        for iter = 1:iterations
            % the sensor readings
            z = awgn(sensor_data.p(:,record_start:end),noise);

            Ns = snaps/2;
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
           
            save("sig_dataGen/sim_signal_"+pos+"_"+noise+"_"+iter+"_.mat",...
            'z_','k','m','Q','theta','radius','radial_units','radial_bound',...
            'canvas_units','R','sen_pos_rad','sim_loc_rad','Y',...
            'P','sim_loc_db_rad');
        
           
        end
    end    
end
