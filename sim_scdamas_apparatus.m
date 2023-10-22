clear
close all
clc

% =========================================================================
% APPARATUS
% =========================================================================

snr_db = 0;
side = 2*25.4 * 1e-2;
dx = 2.54 * 1e-2;
dy = dx;
sen = 33;
R = 19.71/2 * 1e-2;
depth = 1.52;
windows = 10000;

speed = 40;
freq = 15000;

k = 2*pi*freq/speed;

% =========================================================================
% GRIDSPACE
% =========================================================================

Nx = side/dx+1;
Ny = Nx;

X_grid = -(side)/2:dx:(side)/2;
Y_grid = -(side)/2:dy:(side)/2;
[X,Y] = meshgrid(X_grid,Y_grid);

% =========================================================================
% SENSORS 
% =========================================================================

sen_pos = zeros(sen,3);
% for s=1:sen
%     if(s==1)
%         [x,y] = pol2cart(0, 0);
%     elseif(s<=9)
%         [x,y] = pol2cart(pi/4*mod(s-1,8), R/8*1);
%     elseif(s<=17)
%         [x,y] = pol2cart(pi/4*mod(s-1,8)+pi/8, R/8*2);
%     elseif(s<=25)
%         [x,y] = pol2cart(pi/4*mod(s-1,8), R/8*4);
%     elseif(s<=33)
%         [x,y] = pol2cart(pi/4*mod(s-1,8)+pi/8, R/8*8);
%     end
%     sen_pos(s,:) = [x y 0];
% end

for s=1:sen
    [x,y] = pol2cart(2*pi/sen*s, side);
    sen_pos(s,:) = [x y 0];
end
depth = 0;

% =========================================================================
% SOURCE POSITIONS
% =========================================================================

% src_pos = [-side/4 0 depth; 0 -side/4 depth; 0 0 depth; 0 side/4 depth; side/4 0 depth;];
num_src = 8;
rng(42);
src_pos = [randi([-floor(Nx/2),floor(Nx/2)],num_src,1)*dx randi([-floor(Nx/2),floor(Nx/2)],num_src,1)*dy repmat(depth,num_src,1)];

grid_pos = [reshape(X,[],1) reshape(-Y,[],1) repmat(depth,(Nx)*(Ny),1)];
grids = size(grid_pos,1);
num_src = size(src_pos,1);

for snr=snr_db
    for src=num_src:-1:1
        
        % =========================================================================
        % RECORDINGS 
        % =========================================================================

        A = zeros(sen,grids);
        for s=1:sen
            r = vecnorm(grid_pos - sen_pos(s,:),2,2);
            A(s,:) = 1./(4*pi*r).*exp(1i*k*r);
%             A(s,:) = 1i/4*besselh(0,1,k*r);
        end
        
        S = zeros(grids,windows);
        for s=1:src
            [~,src_ind] = min(vecnorm((src_pos(s,:) - grid_pos),2,2));
            var = 1;
            S(src_ind,:) = sqrt(var/2)*(randn(1,windows)+1i*randn(1,windows));
        end

        y = awgn(A*S,snr); 

        C = 1./sqrt(sum(abs(A).^2,1));
        a = A*diag(C);
        a_til = (A./(abs(A).^2))*diag(1./C);
        A_til = 1/(sen)^2 * abs(a_til'*a).^2;
        
        % dr = a*a' - diag(diag(a*a'));       % Diagonal Removal
        % A_til = 1/(sen^2-sen)*a_til'*(dr)*a_til;

        Ra = 1/windows*(y*y');
        [V,D] = eig(Ra);
        lambda = trace(diag(D) - min(diag(D))*eye(size(D,1)));
        % Ra = Ra-diag(diag(Ra));         % Diagonal Removal

        y_til = zeros(grids,1);
        for n=1:grids
            y_til(n,:) = 1/sen^2*(a_til(:,n))'*(Ra)*(a_til(:,n));
        end

%         % =========================================================================
%         % DAMAS 
%         % =========================================================================
% 
%         disp('DAMAS : Gauss-Siedel');
%         disp('computing...');
% 
%         D = diag(diag(A_til));
%         L = tril(A_til) - D;
%         U = triu(A_til) - D;
% 
%         e = max(eig(-inv(D+L)*(U)));
%         disp("absolute of largest eigenvalue: "+abs(e));
%         if abs(e) >= 1
%             disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1'); 
%             disp ('this process is not convergent.');
%         end
% 
%         X_gs(:,1) = rand(grids,1);
%         thresh = 1e-5*ones(grids,1);
%         iter = 1;
%         err = 1e8*rand(grids,1);
%         while sum(abs(err)>=thresh) ~= zeros(grids,1)
%             X_gs(:,iter+1) = max(0,(D+L)\y_til - (D+L)\(U)*X_gs(:,iter)); 
%             err = X_gs(:,iter+1) - X_gs(:,iter);
%             iter = iter+1;
%         end
% 
%         disp("final error:  "+sum(abs(err)));
%         disp("num_iterations:  "+iter);
% 
%         x_damas = reshape(abs(X_gs(:,iter)),Nx,Ny);
% 
%         disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');
% % 
        % =========================================================================
        % SC DAMAS 
        % =========================================================================

        disp('SC-DAMAS : SEDUMI');
        disp('computing...');

        cvx_begin quiet
            variable x_scdamas(grids,1);
            minimize( norm( y_til - A_til*x_scdamas, 2 ) );
            subject to
                x_scdamas >= 0;
                norm( x_scdamas, 1 ) <= abs(lambda);
        cvx_end

        x_scdamas = reshape(abs(x_scdamas),Nx,Ny);

        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');

        % =========================================================================
        % CMF 
        % =========================================================================

        disp('CMF');
        disp('computing...');
        x_cmf = zeros(grids,1);
        err = 1e8;
        best_sig = 0;

        for sig=0:0.5:10
        %     disp("sigma = "+sig);
            cvx_begin quiet
                variables x(grids,1);
                minimize( norm( Ra - A*diag(x)*A' - (sig^2)*eye(sen), 'fro' ) );
                subject to
                    x_cmf >= 0;
                    sum(x_cmf) <= abs(lambda);
                    sig^2 >= 0;
            cvx_end
            perf = norm( Ra - A*diag(x)*A' - (sig^2)*eye(sen), 'fro' );
            if perf<err
                err = perf;
                x_cmf = x;
                best_sig = sig;
            end
        end

        x_cmf = reshape(abs(x_cmf),Nx,Ny);

        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');

        % =========================================================================
        % DAS
        % =========================================================================

        disp('DAS');
        disp('computing...');

        % initialising the MV spectrum matrix 
        x_das = zeros(Nx,Ny);
        Rm = zeros(Nx,Ny,sen);

        for i=1:Nx
            for j=1:Ny
                [~,min_ind] = min(vecnorm(grid_pos-[X(i,j) Y(i,j) depth],2,2));
                Rm(i,j,:) = vecnorm(grid_pos(min_ind,:) - sen_pos,2,2);
%                 c = 1./sen*exp(-1i*k*reshape(Rm(i,j,:),[],1));
                c = a_til(:,min_ind);
                x_das(i,j) = real(c'*Ra*c); 
            end
        end
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');

        % =========================================================================
        % MVDR
        % =========================================================================

        disp('MVDR');
        disp('computing...');

        % initialising the MV spectrum matrix 
        x_mvdr = zeros(Nx,Ny);
        RaI = pinv(Ra);

        for i=1:Nx
            for j=1:Ny
                [~,min_ind] = min(vecnorm(grid_pos-[X(i,j) Y(i,j) depth],2,2));
%                 c = 1./sen*exp(-1i*k*reshape(Rm(i,j,:),[],1));
                c = A(:,min_ind);
                x_mvdr(i,j) = real(1/(c'*RaI*c)); 
            end
        end

        % min_v = min(min([x_damas; x_scdamas; x_cmf; x_das; x_mvdr]));
        % max_v = max(max([x_damas; x_scdamas; x_cmf; x_das; x_mvdr]));
        % 
        % x_damas = (x_damas - min_v) / max_v;
        % x_scdamas = (x_scdamas - min_v) / max_v;
        % x_das = (x_das - min_v) / max_v;
        % x_cmf = (x_cmf - min_v) / max_v;
        % x_mvdr = (x_mvdr - min_v) / max_v;
        
%         x = abs(x_damas);
%         pks = find(imregionalmax(x));
%         [val_sorted, ind] = sort(x(pks), 'descend');
%         valk = val_sorted(1:min(length(pks),src));
%         pks_sorted = pks(ind);
%         pksk = pks_sorted(1:src);
%         fig=figure;
%         imagesc(X_grid,Y_grid,flip(x_damas,1));
%         set(gca,'YDir','normal');
%         xticks([-(side-dx)/2 dx/2 (side-dx)/2]);
%         yticks([-(side-dx)/2 dy/2 (side-dx)/2]);
%         grid on;
%         title('DAMAS');
%         colorbar;
%         colormap(flipud(bone));
%         hold on;
%         scatter(src_pos(:,1),src_pos(:,2),100,'o','filled','g');
%         scatter(X(pksk),-Y(pksk),30,'o','filled','b');
%         legend('Sources', 'Top peaks'); 
%         saveas(fig,"report figures/set1_damas_8_n"+snr,"epsc");

        x = abs(x_scdamas);
        pks = find(imregionalmax(x));
        [val_sorted, ind] = sort(x(pks), 'descend');
        valk = val_sorted(1:min(length(pks),src));
        pks_sorted = pks(ind);
        pksk = pks_sorted(1:src);
        fig=figure;
        imagesc(X_grid,Y_grid,flip(x_scdamas,1));
        set(gca,'YDir','normal');
        xticks([-(side-dx)/2 dx/2 (side-dx)/2]);
        yticks([-(side-dx)/2 dy/2 (side-dx)/2]);
        grid on;
        title('SC DAMAS');
        colorbar;
        colormap(flipud(bone));
        hold on;
        scatter(src_pos(:,1),src_pos(:,2),100,'o','filled','g');
        scatter(X(pksk),-Y(pksk),30,'o','filled','b');
        legend('Sources', 'Top peaks'); 
        saveas(fig,"report figures/set1_scdamas_8_n"+snr,"epsc");

        x = abs(x_cmf);
        pks = find(imregionalmax(x));
        [val_sorted, ind] = sort(x(pks), 'descend');
        valk = val_sorted(1:min(length(pks),src));
        pks_sorted = pks(ind);
        pksk = pks_sorted(1:min(size(pks_sorted,1),src));
        fig=figure;
        imagesc(X_grid,Y_grid,flip(x_cmf,1));
        set(gca,'YDir','normal');
        xticks([-(side-dx)/2 dx/2 (side-dx)/2]);
        yticks([-(side-dx)/2 dy/2 (side-dx)/2]);
        grid on;
        title("CMF (sig = "+best_sig+")");
        colorbar;
        colormap(flipud(bone));
        hold on;
        scatter(src_pos(:,1),src_pos(:,2),100,'o','filled','g');
        scatter(X(pksk),-Y(pksk),30,'o','filled','b');
        legend('Sources', 'Top peaks'); 
        saveas(fig,"report figures/set1_cmf_8_n"+snr,"epsc");
        
        x = abs(x_das);
        pks = find(imregionalmax(x));
        [val_sorted, ind] = sort(x(pks), 'descend');
        valk = val_sorted(1:min(length(pks),src));
        pks_sorted = pks(ind);
        pksk = pks_sorted(1:src);
        fig=figure;
        imagesc(X_grid,-Y_grid,flip(abs(x_das),1));
        set(gca,'YDir','normal');
        xticks([-(side-dx)/2 dx/2 (side-dx)/2]);
        yticks([-(side-dx)/2 dy/2 (side-dx)/2]);
        grid on
        title('DAS');
        colorbar;
        colormap(flipud(bone));
        hold on;
        scatter(src_pos(:,1),src_pos(:,2),100,'o','filled','g');
        scatter(X(pksk),Y(pksk),30,'o','filled','b');
        legend('Sources', 'Top peaks'); 
        saveas(fig,"report figures/set1_das_8_n"+snr,"epsc");

        x = abs(x_mvdr);
        pks = find(imregionalmax(x));
        [val_sorted, ind] = sort(x(pks), 'descend');
        valk = val_sorted(1:min(length(pks),src));
        pks_sorted = pks(ind);
        pksk = pks_sorted(1:src);
        fig=figure;
        imagesc(X_grid,-Y_grid,flip(abs(x_mvdr),1));
        set(gca,'YDir','normal');
        xticks([-(side-dx)/2 dx/2 (side-dx)/2]);
        yticks([-(side-dx)/2 dy/2 (side-dx)/2]);
        grid on;
        title('MVDR');
        colorbar;
        colormap(flipud(bone));
        hold on;
        scatter(src_pos(:,1),src_pos(:,2),100,'o','filled','g');
        scatter(X(pksk),Y(pksk),30,'o','filled','b');
        legend('Sources', 'Top peaks'); 
        saveas(fig,"report figures/set1_mvdr_8_n"+snr,"epsc");

        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');
        
%         saveas(fig,"setup1/Location_Noise_G3/noise_"+snr+"_src_"+src+".jpg");
        break 
    end
    break
end
        
x = abs(x_mvdr);
pks = find(imregionalmax(x));
[val_sorted, ind] = sort(x(pks), 'descend');
valk = val_sorted(1:min(length(pks),src));
pks_sorted = pks(ind);
pksk = pks_sorted(1:src);

x_loc = src_pos(:,1);
y_loc = src_pos(:,2);
[~,xorder] = sort(x_loc);
[~,porder] = sort(X(pksk));
px_loc = X(pksk);
py_loc = Y(pksk);

fig=figure;
mesh(X,Y,x_mvdr+depth);
ylabel('x-position [m]');
xlabel('y-position [m]');
title('Setup');
hold on;
scatter3(sen_pos(:,1),sen_pos(:,2),sen_pos(:,3),'.','k');
scatter3(x_loc(xorder),y_loc(xorder),valk(porder)+depth,150,'X','r');
% scatter3(px_loc(porder),py_loc(porder),valk(porder)+depth,150,'X','b');
legend('Grid', 'Microphone', 'Sources');
view(2);
saveas(fig,"report figures/set1_8","epsc");
