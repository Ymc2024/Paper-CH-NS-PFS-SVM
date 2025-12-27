clear;
clc;

dirname = 'ex03_2_SFdataSVM_bdf2';

datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

T = 10;
Tsplit = [0  4  T];
tsave = [0.1 1 ];

time = [];
for kk = 1:length(tsave)
    time = [time, Tsplit(kk):tsave(kk):Tsplit(kk+1)];
end

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

time = [0,0.1,0.2,0.3];
for k = 1:length(time)
    t = time(k)
    ssp = [datadir '/phi_t=' num2str(t) '.txt'];
    ssr = [datadir '/rho_t=' num2str(t) '.txt'];
    phi = load(ssp); rho = load(ssr); 


    figure(1)
    % subplot(1,2,1)
    s=pcolor(X,Y,phi);
%     s=mesh(X,Y,phi);
    s.FaceColor='interp';
    s.EdgeColor='interp';
%     view(0,90);
%     colormap default
    colormap jet;
%     colorbar;
%     colormap viridis;
    axis square;
    axis tight;
    axis equal;
%     caxis([-1 1])
    axis off;
    title(['t=',num2str(t)]);
    drawnow;

    figname = ['../../09figures/',dirname,'_phi_t=' num2str(t) '.png'];
    % print(figname,'-dpng', '-r300')

    figure(2)
    % subplot(1,2,2)
    s=pcolor(X,Y,rho);
%     s=mesh(X,Y,phi);
    s.FaceColor='interp';
    s.EdgeColor='interp';
%     view(0,90);
%     colormap default
    colormap jet;
%     colorbar;
%     colormap viridis;
    axis square;
    axis tight;
    axis equal;
%     caxis([-1 1])
    axis off;
    title(['t=',num2str(t)]);
    
    drawnow;

    figname = ['../../09figures/',dirname,'_rho_t=' num2str(t) '.png'];
    % print(figname,'-dpng', '-r300') 
    pause(1)
end



