clc
clear 
close all

i = load('SOM-Acoustics-2(Normal)\hrtf_Ni50_data.mat');
A = i.E_s;
sen = size(A,1);
src = size(A,2);
azim = 0:180/50:180;

C = zeros(src);
a = zeros(sen,src);
a_til = zeros(sen,src);
for k=1:src
    C(k) = 1;  % norm(a_bar(:,k));
    a(:,k) = A(:,k)/C(k);
    a_til(:,k) = C(k)*A(:,k).*((abs(A(:,k))).^2);
end

A_til = zeros(src,src);
for k=1:src
    for n=1:src
        A_til(k,n) = 1/(sen)^2 * abs(a_til(:,k)'*a(:,n))^2;
    end
end

performance = zeros(src,src);
Y = zeros(src,src,src);
for loc1=1:src
    for loc2=loc1+15:src
        
        I = 100;
        s = zeros(src,I);
        s(loc1) = 1.1553+1.0236i;
        s(loc2) = 1.8602+1.1911i;
        x_til = mean(abs(s).^2,2);
     
        y = awgn(A*s,20);           
        R = y*y';
        [V,D] = eig(R);
        lambda = trace(diag(D) - min(diag(D))*eye(size(D,1)));

        y_til = zeros(src,1);
        for n=1:src
            y_til(n,:) = 1/sen^2*(a_til(:,n))'*(R)*(a_til(:,n));
        end
            
        cvx_begin quiet
            variable x(src,1);
            minimize( norm( y_til - A_til*x, 2 ) );
            subject to
                x >= 0;
                norm( x, 1 ) <= lambda;
        cvx_end

        plot(abs(x_til));
        hold on;
        plot(x)
        hold off;

        pks = find(imregionalmax(x));
        [~, ind] = sort(x(pks), 'descend');
        pks_sorted = pks(ind);
        pksk = pks_sorted(1:2);

%         disp(azim([loc1,loc2]));
        disp(azim(pksk));
        
        performance(loc1,loc2) = norm(sort([loc1;loc2])-sort(pksk));
        pause(0.3);
%         break;
    end
    
end

imagesc(performance);