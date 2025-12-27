clear;
clc;


dirname = 'ex05_1_SFdataSVM_bdf2';


datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

T = 50;
Tsplit = [0  T];
tsave = [1 ];  % ex03_1_SFdata
Tsplit = [0.1 4 10 T];
tsave = [0.1 1 5];
% Tsplit = [0 10 50];
% tsave = [1 5 10]; 

time = [];
for kk = 1:length(tsave)
    time = [time, Tsplit(kk):tsave(kk):Tsplit(kk+1)];
end

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

time = [1.8,2.5,5,10];
time = [1,1.5,2.5,50];
% time =20;
for k = 1:length(time)
    t = time(k)
    ss = [datadir '/phi_t=' num2str(t) '.txt'];
    ssr = [datadir '/rho_t=' num2str(t) '.txt'];

    phi = load(ss); rho = load(ssr); 

    ssu = [datadir '/u_t=' num2str(t) '.txt'];
    ssv = [datadir '/v_t=' num2str(t) '.txt'];
    ssp = [datadir '/p_t=' num2str(t) '.txt'];
    u = load(ssu); v = load(ssv); p = load(ssp);

    figure(1)
    % subplot(1,2,1)
    s=pcolor(X,Y,phi);
%     s=mesh(X,Y,phi);
    s.FaceColor='interp';
    s.EdgeColor='interp';
    colormap jet;
    hold on

    temp = 4;
    quiver(X(1:temp:end,1:temp:end),Y(1:temp:end,1:temp:end),...
            u(1:temp:end,1:temp:end),v(1:temp:end,1:temp:end),...
            1,'Color','w','linewidth',1.5)
    % axis square;
    box on;

    hold off
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

    s.FaceColor='interp';
    s.EdgeColor='interp';

    colormap jet;

    hold on

    temp = 4;
    quiver(X(1:temp:end,1:temp:end),Y(1:temp:end,1:temp:end),...
            u(1:temp:end,1:temp:end),v(1:temp:end,1:temp:end),...
            1,'Color','w','linewidth',1.5)
    % axis square;
    box on;

    hold off
    axis square;
    axis tight;
    axis equal;
%     caxis([-1 1])
    axis off;
    title(['t=',num2str(t)]);
   
    
    drawnow;
    figname = ['../../09figures/',dirname,'_rho_t=' num2str(t) '.png'];
    % print(figname,'-dpng', '-r300')
    
% 
    pause(0.1)
end



