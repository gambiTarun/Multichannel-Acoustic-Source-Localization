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
            y = mean(data.z_,2);
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
            
            epochs = 10;
            
            [x_db,y_db] = pol2cart(sim_loc_db_rad(:,2)*pi/180, sim_loc_db_rad(:,1)*R/radial_units);
            [x_sen,y_sen] = pol2cart(sen_pos_rad(:,2)*pi/180, sen_pos_rad(:,1)*R/radial_units);
%             [x_src,y_src] = pol2cart(src_pos_rad(:,2)*pi/180, src_pos_rad(:,1)*R/radial_units);
            [x_loc,y_loc] = pol2cart(sim_loc_rad(:,2)*pi/180, sim_loc_rad(:,1)*R/radial_units);

            ATF = zeros(sen,db_size);
            for k_ind=1:db_size
                ATF(:,k_ind) = 1i/4*besselh(0,1,k_ind*vecnorm(([x_sen y_sen] - [x_db(k_ind) y_db(k_ind)]),2,2));
            end
            
            Se = zeros(1,src_num);    
            He = zeros(sen,src_num);  
            loss = zeros(1,epochs);
            est_loc = zeros(epochs,src_num);
            for e=1:epochs
                    
                for src=1:src_num
                    Se(:,src) = complex(rand(),rand());                   
                end
                               
                for s=1:src_num
                    He(:,s) = ATF(:,randi(db_size));                   
                end
                
                curr_loss = 0;
                prev_loss = 0;
                
                while(true)
                    
                    for src=1:src_num
                        
                        % opt -> H
                        losses = zeros(1,db_size);
                        for i=1:db_size
                            sum_ = zeros(sen,1);
                            for sum_src=1:src_num
                                if(sum_src==src)
                                    sum_ = sum_ + ATF(:,i);
                                else
                                    sum_ = sum_ + He(:,sum_src);
                                end
                            end
                            losses(i) = norm(y - sum_);
                        end
                        [~,est_loc(e,src)] = min(losses); 
                        He(:,src) = ATF(:,est_loc(e,src));
                                  
                        % opt -> S
                        n = He(:,src)'*(y - (sum(He.*repmat(Se,size(He,1),1),2) - He(:,src)*Se(src)));
                        d = He(:,src)'*He(:,src);
                        Se(src) = n/d;
                        
                    end
                    
                    prev_loss = curr_loss;
                    curr_loss = norm(y - sum(He.*repmat(Se,size(He,1),1),2));
%                     disp(prev_loss+" "+curr_loss);
%                     disp(sim_loc_rad);
%                     disp(sim_loc_db_rad(est_loc,:));
%                     disp("~~~~~~~~~~~~~~~~~~~");
                    if(abs(prev_loss-curr_loss)<1e-10)
                        break;
                    end
                    
                end    
                loss(e) = curr_loss;
                
            end
            [~,min_ind] = min(loss);
            [x_est, y_est] = pol2cart(sim_loc_db_rad(est_loc(min_ind,:),2)*pi/180,...
                sim_loc_db_rad(est_loc(min_ind,:),1)*R/radial_units);

%             % Grid Constraints
%             canvas_grid = zeros(Nx,Ny);
%             [A_grid,B_grid,C_grid] = pol2cart(P,Y,canvas_grid);
%             [x_bound,y_bound] = pol2cart(radial_bound(2,:),radial_bound(1,:));
%                        
%             figure;
%             mesh(A_grid,B_grid,C_grid);
%             ylabel('x-position [m]');
%             xlabel('y-position [m]');
%             title("noise: "+noise_levs(noise)+" db");
%             hold on;          
%             scatter3(x_loc,y_loc,zeros(src_num,1),100,'^','b');
%             scatter3(x_est,y_est,zeros(src_num,1),150,'X','k');
%             scatter3(x_db,y_db,zeros(length(x_db),1),100,'.','y');
%             scatter3(x_bound,y_bound,zeros(length(x_bound),1),10,'o','k');
%             legend('MV Spec', 'Source', 'Top peaks', 'Prior', 'Boundary');
%             view(2);
   
            performance(pos,noise,iter) = norm(sort(complex(x_loc,y_loc),'ComparisonMethod','abs')...
                - sort(complex(x_est,y_est),'ComparisonMethod','abs'));
            
        end
    end
end

figure;
imagesc(1:num_pos,noise_levs,mean(performance,3));
xlabel('pos number');
ylabel('noise db');
title('proposed');
colorbar;
savefig('performance_proposed');
save 'performance_proposed' performance