clear
close all
clc

% =========================================================================
% APPARATUS
% =========================================================================

snr = 20;
beamformer = 1;
sen = 40;
R = 0.15;
radial_units = 6;
canvas_units = 9;
speed = 25;
windows = 1000;

freq = 250;
% freq = 1000;

% wavenumber k and the maximum mode M
k = 2*pi*freq/speed;
m = ceil(11/9*k*R);

% =========================================================================
% GRIDSPACE
% =========================================================================

Nx = 129;           % number of grid points in the x (row) direction
Ny = 129;           % number of grid points in the y (column) direction
dx = 2*(R/radial_units*canvas_units)/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]

X_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,Nx);
Y_var = linspace(-R*canvas_units/radial_units,R*canvas_units/radial_units,Ny);

[X_grid,Y_grid] = meshgrid(X_var,Y_var);

[P,Y] = cart2pol(X_grid,Y_grid);
P = mod(P,2*pi);

% =========================================================================
% SENSORS 
% =========================================================================

% [radial_bound,sen_pos_rad] = rectgrid_rad(sen,R,2*R*30/30,2*R*15/30,radial_units,360);
[radial_bound,sen_loc_rad] = circgrid_rad(sen,R,radial_units,360/40*sen);
% [radial_bound,sen_pos_rad] = ellipsogrid_rad(sen,R*0.6,15,3,radial_units*0.6,360);
% [radial_bound,sen_pos_rad] = semicircgrid_rad(sen,R,radial_units,360);

% =========================================================================
% SOURCE POSITIONS
% =========================================================================
% src_loc_rad = [3 220;];
src_loc_rad = [4 160;5 280; 3 0; 4 310; 5 200; 2 45; 2.5 90];
grid_loc_rad = [];

for i=1:Nx
    for j=1:Ny
        [~,min_index] = min(abs(radial_bound(2,:)-P(i,j)));
        if(Y(i,j)<radial_bound(1,min_index))
            grid_loc_rad = [grid_loc_rad;Y(i,j)*radial_units/R P(i,j)*180/pi];
        end
    end
end

% [x,y] = pol2cart(src_pos_rad(:,2)*pi/180, src_pos_rad(:,1)*R/radial_units);
% scatter(x,y);

% =========================================================================
% RECORDINGS 
% =========================================================================

grids = size(grid_loc_rad,1);
src = size(src_loc_rad,1);
[x_grid,y_grid] = pol2cart(grid_loc_rad(:,2)*pi/180, grid_loc_rad(:,1)*R/radial_units);
[x_sen,y_sen] = pol2cart(sen_loc_rad(:,2)*pi/180, sen_loc_rad(:,1)*R/radial_units);
[x_src,y_src] = pol2cart(src_loc_rad(:,2)*pi/180, src_loc_rad(:,1)*R/radial_units);

% A = zeros(sen,grids);
% for k_ind=1:grids
%     A(:,k_ind) = 1i/4*besselh(0,1,k*vecnorm(([x_sen y_sen] - [x_grid(k_ind) y_grid(k_ind)]),2,2));
% end
som_data = load("SOM-Acoustics-2(Normal)\forward_scatter_lungs.mat");
A = som_data.z_u;

% imagesc(abs(A));
% colorbar;
% return;

S = zeros(grids,windows);
for p=1:src
    [~,src_ind] = min(vecnorm(([x_src(p) y_src(p)] - [x_grid y_grid]),2,2));
    var = 1;
    S(src_ind,:) = sqrt(var/2)*(randn(1,windows)+1i*randn(1,windows));
end

% for c = 10:5:40
y = awgn(A*S,snr);

% the covariance of the fourier coefficients
Ra = y*y';
% figure;
% imagesc(real(Ra1));
% colorbar;
% return;

% initialising the MV spectrum matrix 
Z = zeros(Nx,Nx);

RaI = pinv(Ra);
% subplot(2,1,1);
% imagesc(real(Ra));
% colorbar;
% subplot(2,1,2);
% imagesc(real(RaI));
% colorbar;
% end

for i=1:Nx
    for j=1:Ny
        [~,min_index] = min(abs(radial_bound(2,:)-P(i,j)));
        if(Y(i,j)<radial_bound(1,min_index))
            [x_loc,y_loc] = pol2cart(P(i,j),Y(i,j));
            [~,min_i] = min(vecnorm([x_loc,y_loc]-[x_grid, y_grid],2,2));
            c = A(:,min_i);
            Z(i,j) = (c'*RaI*c)^-1;
        end
    end
end

[A,B,C] = pol2cart(P,Y,real(Z));
[x_bound,y_bound] = pol2cart(radial_bound(2,:),radial_bound(1,:));

% predicted location of the source
pks = find(imregionalmax(C));
[val_sorted, ind] = sort(C(pks), 'descend');
valk = val_sorted(1:min(length(pks),src));
pks_sorted = pks(ind);
pksk = pks_sorted(1:min(length(pks),src));

figure;
mesh(A,B,C);
ylabel('x-position [m]');
xlabel('y-position [m]');
zlabel('Z');

[~,xorder] = sort(x_src);
[~,porder] = sort(A(pksk));
px_src = A(pksk);
py_src = B(pksk);

fig = figure;
mesh(A,B,C);
ylabel('y-position [m]');
xlabel('x-position [m]');
hold on;
scatter3(px_src(porder),py_src(porder),valk(porder),150,'X','b');
scatter3(x_src(xorder),y_src(xorder),valk(porder),150,'X','r');
scatter3(x_bound,y_bound,zeros(length(x_bound),1),100,'o','k');
legend('MV Spec', 'Top peaks', 'Source', 'Sensors');
view(2);
