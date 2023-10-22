clear
close all
clc

noise_levs = 0:2:20;
source_radius = 0:0.5:8;

fig = figure;
for radi=1:length(source_radius)
    
    noiseP = zeros(1,length(noise_levs));
    
    for n=1:length(noise_levs)
    
        x = load("accVnoiseVradius/dist_rad"+source_radius(radi)+"u_snr"+noise_levs(n)+"db.mat",'distance');
        
        noiseP(n) = mean(x.distance);
        
    end
    
    subplot(5,4,radi);
    plot(noise_levs,noiseP);
    xlabel('noise');
    ylabel('perf.(dist.(m))');
    title("rad: "+source_radius(radi)+"u");
    axis([0 20 0 0.1442]);
    
end
saveas(fig,char("accVnoiseVradius/plots/snrVperf.fig")); 

fig = figure;
for n=1:length(noise_levs)
    
    radP = zeros(1,length(source_radius));
    
    for radi=1:length(source_radius)
    
        x = load("accVnoiseVradius/dist_rad"+source_radius(radi)+"u_snr"+noise_levs(n)+"db.mat",'distance');
        
        radP(radi) = mean(x.distance);
        
    end
    
    subplot(3,4,n);
    plot(source_radius,radP);
    xlabel('rad(1u=1.87cm)');
    ylabel('perf.(dist.(m))');
    title("snr: "+noise_levs(n)+"db");
    axis([0 8 0 0.1442]);
        
end
saveas(fig,char("accVnoiseVradius/plots/radVperf.fig"));
