clear 
close all
clc

noise_levs = 0:2:20;
source_radius = 0:0.5:8;
iterations = 100;
Q = 80;

sen_pos_rad = zeros(Q,2);
for i=1:Q
    sen_pos_rad(i,:) = [8 i*360/Q];
end

Y = linspace(0, 0.15, 200);
P = linspace(0, 2*pi, 180);
[P_grid,Y_grid] = meshgrid(P,Y);
[A,B] = pol2cart(P_grid,Y_grid);

for radi=1:length(source_radius)
     
    fig = figure;
%     sgtitle("radius: "+source_radius(radi));
    disp("radius: "+source_radius(radi));
    for n=1:length(noise_levs)
        
        [x_loc,y_loc] = pol2cart(pi,source_radius(radi)*0.15/8);
        
        subplot(3,4,n);
        scatter(x_loc,y_loc,200,'X','r');
        ylabel('y-position [m]');
        xlabel('x-position [m]');
        title("snr: "+noise_levs(n)+"db");
        axis([-0.15 0.15 -0.15 0.15]);
        hold on;
        grid on;
        
        for i=1:Q
            [x,y] = pol2cart(sen_pos_rad(i,2)*pi/180,sen_pos_rad(i,1)*0.15/8);
            scatter(x,y,50,'o','b');
            hold on
        end
        
        for iter=1:iterations
            
            x = load("accVnoiseVradius/MV_iter"+iter+"_rad"+source_radius(radi)+"u_snr"+noise_levs(n)+"db.mat",'C');
            
            C = x.C;
            
            [xmax,ymax] = find(abs(C) == max(max(abs(C))));
            xmax = xmax(1);
            ymax = ymax(1);

%             figure(6);
%             mesh(A,B,abs(C));
%             ylabel('x-position [m]');
%             xlabel('y-position [m]');
%             zlabel('Z');
            
            scatter(A(xmax,ymax),B(xmax,ymax),10,'*','k');
            hold on;
            
        end
%         break;
    end
    savefig(fig,char("accVnoiseVradius/plots/peaks_"+source_radius(radi)+"u.fig"));
    saveas(fig,char("accVnoiseVradius/plots/peaks_"+source_radius(radi)+"u.jpeg"));
    close all;
%     break;
end
