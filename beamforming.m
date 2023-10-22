clear
close all
clc

num_pos = 10;
noise_levs = 40:-5:0;
iterations = 100;

performance = zeros(num_pos,length(noise_levs),iterations);

for pos=1:num_pos
    
    for noise=1:length(noise_levs)

        pred_perf = zeros(iterations,1);
        
        for iter=1:iterations
            
            disp(pos+"_"+noise_levs(noise)+"_"+iter);
            
            data = load("sig_dataGen\sim_signal_"+pos+"_"+noise_levs(noise)+"_"+iter+"_.mat");
            z_ = data.z_;
            k = data.k;
            m = data.m;
            theta = data.theta;
            Q = data.Q;
            radius= data.radius;
            R = data.R;
            radial_units = data.radial_units;
            canvas_units = data.canvas_units;
            radial_bound = data.radial_bound;
            sim_loc_rad = data.sim_loc_rad;
            Y = data.Y;
            P = data.P;
            Nx = size(Y,1);
            Ny = size(Y,2);

            beamformer = 1;

            [sx_loc,sy_loc] = pol2cart(sim_loc_rad(:,2)*pi/180, sim_loc_rad(:,1)*R/radial_units);

            %The matrix for fourier transformation
            ind1 = repmat(-m:m,Q,1);    % varying M
            ind21 = repmat(theta',1,2*m+1);     % varying theta
            ind22 = repmat(radius',1,2*m+1);    % varying radii

            gamma = besselh(ind1,1,k*ind22).*exp(1i*ind1.*ind21);

            % computing alpha
            alpha = 4/(1i)*pinv(gamma)*z_;

            % the covariance of the fourier coefficients
            Ra = 1/size(alpha,2)*alpha*ctranspose(alpha);

            % initialising the MV spectrum matrix 
            Z = zeros(Nx,Ny);

            RaI = pinv(Ra);

            for i=1:Nx
                for j=1:Ny
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
            valk = val_sorted(1:min(length(pks),size(sim_loc_rad,1)));
            pks_sorted = pks(ind);
            pksk = pks_sorted(1:min(length(pks),size(sim_loc_rad,1)));
            % xmax = xmax(1);
            % ymax = ymax(1);

%             figure;
%             mesh(A,B,C);
%             ylabel('y-position [m]');
%             xlabel('x-position [m]');
%             title("noise: "+noise_levs(noise)+" db");
%             hold on;
%             scatter3(sx_loc,sy_loc,ones(length(sx_loc),1)*valk(1),150,'X','b');
%             scatter3(A(pksk),B(pksk),valk,150,'X','k');
%             scatter3(x_bound,y_bound,zeros(length(x_bound),1),100,'o','k');
%             legend('MV Spec', 'Source', 'Top peaks', 'Boundary');
%             view(2);
            
            performance(pos,noise,iter) = norm(sort(complex(sx_loc,sy_loc),'ComparisonMethod','abs')...
                - sort(complex(A(pksk),B(pksk)),'ComparisonMethod','abs'));
            
        end
        
    end
    
end

figure;
imagesc(1:num_pos,noise_levs,mean(performance,3));
xlabel('pos number');
ylabel('noise db');
title('beamforming');
colorbar;
savefig('performance_mvdr_beamformer');
save 'performance_mvdr_beamformer' performance