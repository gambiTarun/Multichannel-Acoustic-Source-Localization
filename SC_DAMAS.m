clc
clear 
close all

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
            sen = data.Q;
            sen_pos_rad = data.sen_pos_rad;
            sim_loc_rad = data.sim_loc_rad;
            R = data.R;
            radial_units = data.radial_units;
            src_num = size(sim_loc_rad,1);
            radial_bound = data.radial_bound;
            Y = data.Y;
            P = data.P;
            Nx = size(Y,1);
            Ny = size(Y,2);
            
            % prior data
            sim_loc_db_rad = data.sim_loc_db_rad;
            db_size = size(sim_loc_db_rad,1);
                       
            [x_db,y_db] = pol2cart(sim_loc_db_rad(:,2)*pi/180, sim_loc_db_rad(:,1)*R/radial_units);
            [x_sen,y_sen] = pol2cart(sen_pos_rad(:,2)*pi/180, sen_pos_rad(:,1)*R/radial_units);
%             [x_src,y_src] = pol2cart(src_pos_rad(:,2)*pi/180, src_pos_rad(:,1)*R/radial_units);
            [x_loc,y_loc] = pol2cart(sim_loc_rad(:,2)*pi/180, sim_loc_rad(:,1)*R/radial_units);

            A = zeros(sen,db_size);
            for k_ind=1:db_size
                A(:,k_ind) = 1i/4*besselh(0,1,k_ind*vecnorm(([x_sen y_sen] - [x_db(k_ind) y_db(k_ind)]),2,2));
            end
            
            C = 1./sqrt(sum(abs(A).^2,1));
            a = A*diag(C);
            a_til = (A./(abs(A).^2))*diag(1./C);
            A_til = 1/(sen)^2 * abs(a_til'*a).^2;
            
            y = data.z_;           
            Ra = y*y';
            [V,D] = eig(Ra);
            lambda = trace(diag(D) - min(diag(D))*eye(size(D,1)));

            y_til = zeros(db_size,1);
            for n=1:db_size
                y_til(n,:) = 1/sen^2*(a_til(:,n))'*(Ra)*(a_til(:,n));
            end
            
            cvx_begin quiet
                variable x(db_size,1);
                minimize( norm( y_til - A_til*x, 2 ) );
                subject to
                    x >= 0;
                    norm( x, 1 ) <= lambda;
            cvx_end
            
            pks = find(imregionalmax(x));
            [val_sorted, ind] = sort(x(pks), 'descend');
            valk = val_sorted(1:min(length(pks),src_num));
            pks_sorted = pks(ind);
            pksk = pks_sorted(1:src_num);
            
            % Grid Contraints
            canvas_grid = zeros(Nx,Ny);
            [A_grid,B_grid,C_grid] = pol2cart(P,Y,canvas_grid);
            [x_bound,y_bound] = pol2cart(radial_bound(2,:),radial_bound(1,:));
                       
%             figure;
%             mesh(A_grid,B_grid,C_grid);
%             ylabel('x-position [m]');
%             xlabel('y-position [m]');
%             title("noise: "+noise_levs(noise)+" db");
%             hold on;          
%             scatter3(x_loc,y_loc,zeros(length(x_loc),1),100,'^','b');
%             scatter3(x_db(pksk),y_db(pksk),zeros(length(pksk),1),150,'X','k');
%             scatter3(x_db,y_db,zeros(length(x_db),1),100,'.','y');
%             scatter3(x_bound,y_bound,zeros(length(x_bound),1),10,'o','k');
%             legend('MV Spec', 'Source', 'Top peaks', 'Prior', 'Boundary');
%             view(2);
   
            performance(pos,noise,iter) = norm(sort(complex([x_loc],[y_loc]),'ComparisonMethod','abs')...
                - sort(complex([x_db(pksk)],[y_db(pksk)]),'ComparisonMethod','abs'));
            
        end
    end  
end

figure;
imagesc(1:num_pos,noise_levs,mean(performance,3));
xlabel('pos number');
ylabel('noise db');
title('SC-DAMAS');
colorbar;
savefig('performance_sc_damas');
save 'performance_sc_damas' performance