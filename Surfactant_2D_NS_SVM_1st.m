function [phi,rho,u,v,p,alpha] = Surfactant_2D_NS_SVM_1st(pde,domain,Nx,Ny,time,option)

global dt kx ky k2 k4 hx hy S1 S2...
    Mrho Mphi xi  epsilon eta theta    ...
    alpha1 alpha2 Sigma upsilon beta...
    gamma1 gamma2 Lambda_1 Lambda_2

S1 = pde.S1;
S2 = pde.S2;
Mrho = pde.Mrho;
Mphi = pde.Mphi;
xi = pde.xi;
gamma1 = pde.gamma1;
gamma2 = pde.gamma2;
epsilon = pde.epsilon;
eta   = pde.eta;
theta = pde.theta;
alpha1=pde.alpha1;
alpha2=pde.alpha2;
Sigma=pde.Sigma;
upsilon = pde.upsilon;
% beta = pde.beta;
Lambda_1 = pde.Lambda_1;
Lambda_2 = pde.Lambda_2;

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-10;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end
if ~isfield(option,'plotflag')
    option.plotflag = 0;   
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;  
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;  
end
if ~isfield(option,'printflag')
    option.printflag = 0;   
end
if ~isfield(option,'vtkflag')
    option.printflag = 0;   
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;   
end

if 1 == option.energyflag
    figname_mass = [pde.name,num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,num2str(time.dt),'_energy.txt'];    
    out1 = fopen(figname_mass,'w');
    out2 = fopen(figname_energy,'w');
end

%% options for fsolve
opts = optimoptions('fsolve','Display','off',...
                       'StepTolerance', option.tol, ...
                       'FunctionTolerance',option.tolit,...
                       'MaxIterations',option.maxit);

%%
T  = time.T;
t  = time.t0;
dt = time.dt;
if isfield(time, 'Tsplit')
    Tsplit = time.Tsplit;
else
    Tsplit = [];
end
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name 'SVM_1st/data'];

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;
% x  = domain.left   + hx*(1:Nx);
% y  = domain.bottom + hy*(1:Ny);
x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);   
k2x = k_x.^2; k2y = k_y.^2;
[kxx, kyy] = ndgrid(k2x,k2y);
k2 = kxx + kyy;
[kx,ky]=ndgrid(k_x,k_y);
k4 = k2.^2;

[xx,yy] = ndgrid(x,y);
phi0 = pde.init_phi(xx,yy);
rho0 = pde.init_rho(xx,yy);
u0 = pde.init_u(xx,yy);
v0 = pde.init_v(xx,yy);
p0 = pde.init_p(xx,yy);  


%% plot initial value
nfigure =1;
if 1 == option.saveflag
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
    end
    ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
    ss2 = [dir_data '/rho_t=' num2str(t) '.txt'];
    ssu = [dir_data '/u_t=' num2str(t) '.txt'];
    ssv = [dir_data '/v_t=' num2str(t) '.txt'];
    ssp = [dir_data '/p_t=' num2str(t) '.txt'];
    writematrix(phi0,ss1,'Delimiter',' ');
    writematrix(rho0,ss2,'Delimiter',' ');
    writematrix(u0,ssu,'Delimiter',' ');
    writematrix(v0,ssv,'Delimiter',' ');
    writematrix(p0,ssp,'Delimiter',' ');
    writematrix(xx,[dir_data '/X.txt'],'Delimiter',' ');
    writematrix(yy,[dir_data '/Y.txt'],'Delimiter',' ');
end
if 1 == option.plotflag
    if 1 == option.saveflag
        subplot(1,2,1)
        showsolution_2D(nfigure,xx,yy,phi0,t,dir_fig);
        subplot(1,2,2)
        showsolution_2D(nfigure,xx,yy,rho0,t,dir_fig);
    else
        subplot(1,2,1)
        showsolution_2D(nfigure,xx,yy,phi0,t);
        subplot(1,2,2)
        showsolution_2D(nfigure,xx,yy,rho0,t);
    end
end

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

if 1 == option.energyflag
    calculate_energy(out1,out2,hx,hy,t,u0,v0,phi0,rho0);
end

iter_nsave = 1;
for nt = 1:nplot
    t = t+dt;   
    phi_star = phi0; 
    rho_star = rho0;
    u_star = u0; v_star = v0; p_star = p0;
    
    if isfield(pde,'exact_phi') && isfield(pde,'exact_rho') && isfield(pde,'exact_u') && isfield(pde,'exact_v')
        rhs_mu =  - Lambda_1 *gamma1.*lap_diff(pde.exact_phi(xx,yy,t)) + E_derivative_phi(pde.exact_phi(xx,yy,t),pde.exact_rho(xx,yy,t));
        tmp1 = pde.phi_t(xx,yy,t) + pde.advection1(xx,yy,t) - Mphi .* lap_diff(rhs_mu);
        rhs_omega =  - Lambda_2 *gamma2.*lap_diff(pde.exact_rho(xx,yy,t)) + E_derivative_rho(pde.exact_phi(xx,yy,t),pde.exact_rho(xx,yy,t));
        tmp2 = pde.rho_t(xx,yy,t) + pde.advection2(xx,yy,t) - Mrho .* lap_diff(rhs_omega);
        tmpu = pde.u_t(xx,yy,t) + pde.convective1(xx,yy,t)...
             + pde.px(xx,yy,t) - upsilon.* pde.lap_u(xx,yy,t) ...
             + pde.exact_phi(xx,yy,t).*diff_x(rhs_mu) + pde.exact_rho(xx,yy,t).*diff_x(rhs_omega);
        tmpv = pde.v_t(xx,yy,t) + pde.convective2(xx,yy,t)...
             + pde.py(xx,yy,t) - upsilon.* pde.lap_v(xx,yy,t) ...
             + pde.exact_phi(xx,yy,t).*diff_y(rhs_mu) + pde.exact_rho(xx,yy,t).*diff_y(rhs_omega);
    else
        tmp1 = 0;
        tmp2 = 0;
        tmpu = 0;
        tmpv = 0;
    end
    
    H1 = phi0./(Mphi*dt) - (diff_x(u_star.*phi_star) + diff_y(v_star.*phi_star))./Mphi + lap_diff(E_derivative_phi(phi_star,rho_star)) - S1/epsilon.^2*lap_diff(phi_star);
    H2 = rho0./(dt*Mrho) - (diff_x(u_star.*rho_star) + diff_y(v_star.*rho_star))./Mrho + lap_diff(E_derivative_rho(phi_star,rho_star)) - S2/eta.^2*lap_diff(rho_star);
    
    phi_star_new = inv_A(H1 + tmp1./Mphi);
    rho_star_new = inv_B(H2 + tmp2./Mrho);


    mu_phi_star = -Lambda_1 *gamma1*lap_diff(phi_star_new) + E_derivative_phi(phi_star,rho_star) + S1/epsilon.^2*(phi_star_new - phi0);
    mu_rho_star = -Lambda_2 *gamma2*lap_diff(rho_star_new) + E_derivative_rho(phi_star,rho_star) + S2/eta.^2*(rho_star_new - rho0);
    
    H3 = u0/dt - (u0.*diff_x(u0)+v0.*diff_y(u0)) - diff_x(p0) - phi_star_new.*diff_x(mu_phi_star) - rho_star_new.*diff_x(mu_rho_star);
    H4 = v0/dt - (u0.*diff_x(v0)+v0.*diff_y(v0)) - diff_y(p0) - phi_star_new.*diff_y(mu_phi_star) - rho_star_new.*diff_y(mu_rho_star);

    utilde_star = inv_C(H3 + tmpu);
    vtilde_star = inv_C(H4 + tmpv);
    
    tmp = (diff_x(utilde_star) + diff_y(vtilde_star))./dt + lap_diff(p0);
    p_star_new_hat = fft2(tmp)./k2;
    p_star_new_hat(1,1) = 0;
    p_star_new = real(ifft2(p_star_new_hat));
    u_star_new = utilde_star - dt.*diff_x(p_star_new - p0);
    v_star_new = vtilde_star - dt.*diff_y(p_star_new - p0);

    % Step 2
    

    h11 = phi0/(dt*Mphi) - (diff_x(u_star_new.*phi_star_new)+diff_y(v_star_new.*phi_star_new))./Mphi + lap_diff(E_derivative_phi(phi_star_new,rho_star_new)) - S1/epsilon.^2*lap_diff(phi_star);
    % h12 = Mphi.*lap_diff(mu_phi_star)./Mphi;
    if 1 == strcmp('SVM1', option.type)
        h12 = Mphi.*lap_diff(mu_phi_star)./Mphi;
    elseif 1 == strcmp('SVM2', option.type)
        h12 = Mphi.*lap_diff(phi_star_new)./Mphi;
    end

    h21 = rho0/(dt*Mrho) - (diff_x(u_star_new.*rho_star_new)+diff_y(v_star_new.*rho_star_new))./Mrho + lap_diff(E_derivative_rho(phi_star_new,rho_star_new)) - S2/eta.^2*lap_diff(rho_star);
    % h22 = Mrho.*lap_diff(mu_rho_star)./Mrho;
    if 1 == strcmp('SVM1', option.type)
        h22 = Mrho.*lap_diff(mu_rho_star)./Mrho;
    elseif 1 == strcmp('SVM2', option.type)
        h22 = Mrho.*lap_diff(rho_star_new)./Mrho;
    end
    
    phi_new1 = inv_A(h11+tmp1./Mphi);
    phi_new2 = inv_A(h12);    
    rho_new1 = inv_B(h21+tmp2./Mrho);
    rho_new2 = inv_B(h22);  
    
    h31_u = u0/dt - (u_star_new.*diff_x(u_star_new)+v_star_new.*diff_y(u_star_new)) - diff_x(p_star_new) - phi_star_new.*diff_x(mu_phi_star) - rho_star_new.*diff_x(mu_rho_star);
    h31_v = v0/dt - (u_star_new.*diff_x(v_star_new)+v_star_new.*diff_y(v_star_new)) - diff_y(p_star_new) - phi_star_new.*diff_y(mu_phi_star) - rho_star_new.*diff_y(mu_rho_star);
    
    utilde_new1 = inv_C(h31_u + tmpu);
    vtilde_new1 = inv_C(h31_v + tmpv);
    
    h32_u = u_star_new;
    h32_v = v_star_new;
    
    utilde_new2 = inv_C(h32_u);
    vtilde_new2 = inv_C(h32_v);
    
    tmp = (diff_x(utilde_new1) + diff_y(vtilde_new1))./dt + lap_diff(p_star_new);
    p_new_hat = fft2(tmp)./k2;
    p_new_hat(1,1) = 0; 
    p_new = real(ifft2(p_new_hat));

    u_new1 = utilde_new1 - dt.*diff_x(p_new - p_star_new);
    v_new1 = vtilde_new1 - dt.*diff_y(p_new - p_star_new);
    
    u_new2 = utilde_new2;
    v_new2 = vtilde_new2;
    
    % Step 3
    alpha = 0;
    alpha = fsolve(@(alpha) non_fun_alpha(alpha,u_new1,u_new2,u0,v_new1,v_new2,v0,phi_new1,phi_new2,phi0,...
        rho_new1,rho_new2,rho0,mu_phi_star,mu_rho_star,u_star_new,v_star_new,tmp1,tmp2,tmpu,tmpv),alpha,opts);
    
    phi = phi_new1 + alpha * phi_new2;
    rho = rho_new1 + alpha * rho_new2;
    u = u_new1 + alpha * u_new2;
    v = v_new1 + alpha * v_new2;
    p = p_new;
    
%% update phi0 and u0
    phi0 = phi;
    rho0 = rho;
    u0 = u;
    v0 = v;
    p0 = p;
    
    if 1 == option.energyflag
        calculate_energy(out1,out2,hx,hy,t,u,v,phi,rho);
    end
    
    % print and plot
    if  0 == mod(nt,nsave(iter_nsave))
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('alpha=%.3f, ',alpha);
            fprintf('epsilon=%.3f,t=%.4f/%.3f, dt=%.2e, N=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,timeElapsed);
        end

        if 1 == option.saveflag
            if round(t) == Tsplit(iter_nsave)
                iter_nsave = iter_nsave + 1;
            end
            if ~exist(dir_data,'dir')
                mkdir(dir_data);
            end
            ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
            ss2 = [dir_data '/rho_t=' num2str(t) '.txt'];
            ssu = [dir_data '/u_t=' num2str(t) '.txt'];
            ssv = [dir_data '/v_t=' num2str(t) '.txt'];
            ssp = [dir_data '/p_t=' num2str(t) '.txt'];
            writematrix(phi0,ss1,'Delimiter',' ');
            writematrix(rho0,ss2,'Delimiter',' ');
            writematrix(u0,ssu,'Delimiter',' ');
            writematrix(v0,ssv,'Delimiter',' ');
            writematrix(p0,ssp,'Delimiter',' ');

            timefile = [dir_data '/iteration_times.txt'];
            fid = fopen(timefile, 'a');
            fprintf(fid, '%.6f, %.6f, %d, %.2e, %.6f\n', t, alpha, Nx, dt, timeElapsed);
            fclose(fid);
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                subplot(1,2,1)
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
                subplot(1,2,2)
                showsolution_2D(nfigure,xx,yy,rho,t,dir_fig);
            else
                subplot(1,2,1)
                showsolution_2D(nfigure,xx,yy,phi,t);
                subplot(1,2,2)
                showsolution_2D(nfigure,xx,yy,rho,t);                
            end
        end
    end
end

if 1 == option.savefinal
    % Save data
    filename=[option.scheme,'_phi_e',num2str(epsilon),'Mphi',num2str(Mphi),'Mrho',num2str(Mrho),'Nx',num2str(Nx),'Ny',num2str(Ny),'dt',num2str(dt),'.mat'];
    save(filename,'xx','yy','hx','hy','Nx','Ny','dt','T','phi','rho','u','v','p','alpha'); 
end
if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function r = fun_inner(f,g)
global hx hy
    r1 = fft2(f.*g);
    r = r1(1,1)*hx*hy;
end

function r = inv_A(phi)
global dt Mphi k2 S1 epsilon gamma1 Lambda_1 
    r = real(ifft2(fft2(phi)./(1/(dt*Mphi) + Lambda_1 *gamma1.* k2.^2 - S1 / epsilon.^2 * k2)));
end

function r = inv_B(rho)
global dt Mrho k2 S2 eta gamma2 Lambda_2
    r = real(ifft2(fft2(rho)./(1/(dt*Mrho) + Lambda_2 *gamma2 * k2.^2 - S2 / eta.^2 * k2)));
end

function r = inv_C(phi)
global dt upsilon k2
    r = real(ifft2(fft2(phi)./(1/dt - upsilon.* k2)));
end

function lap = lap_diff(phi)
global k2
    lap=real(ifft2((k2.*fft2(phi))));
end

function lap = inv_lap(phi)
global k2
    lap=real(ifft2((fft2(phi)./k2)));
end

function lap = diff_x(phi)
global kx
    lap=real(ifft2((kx.*fft2(phi))));
end

function lap = diff_y(phi)
global ky
    lap=real(ifft2((ky.*fft2(phi))));
end


function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2;
end

function r = E_derivative_phi(phi,rho)
global epsilon    alpha1 alpha2 xi theta Lambda_1
    grad_phi_rho = diff_x((theta * rho - xi * grad_square(phi)) .* diff_x(phi)) ...
                 + diff_y((theta * rho - xi * grad_square(phi)) .* diff_y(phi));   
    r =  Lambda_1/(epsilon^2) * f(phi) + grad_phi_rho;
end

function r = E_derivative_rho(phi,rho)
global  beta  alpha1 alpha2 theta eta Lambda_2 
    r =  Lambda_2/eta^2 * G_derivative(rho) - theta/2 * grad_square(phi);
end
function r = F(phi) 
    r =  1/4 * (phi.^2-1).^2;
end

function r = f(phi)
    r = (phi.^2-1) .* phi;
end

function r = E_no_linear(phi,rho)
global epsilon beta   alpha1 alpha2 theta xi eta Lambda_1 Lambda_2
    r =   Lambda_1/ ( epsilon.^2) .* F(phi) ...
        +  Lambda_2/eta.^2 *G(rho) -theta./2 *rho .* grad_square(phi) ...
        + xi/4 * grad_square(phi).^2;
end

function r = G(rho)
global Sigma
    r = 3.62*rho.^4 - 7.25*rho.^3 + 7.3*rho.^2 - 3.68*rho;
end

function r = G_derivative(rho)
global Sigma
    r = 3.62*4*rho.^3 - 7.25*3*rho.^2 + 7.3*2*rho.^1 - 3.68;
end

function r = Energy(u,v,phi,rho)
global gamma1 gamma2  Lambda_1 Lambda_2
%     r = fun_inner(1, gamma1 .* grad_square(phi)/2 +  gamma2.*grad_square(rho)/2 + E_no_linear(phi,rho));
    r = fun_inner(1,1/2 * (u.^2 + v.^2) + Lambda_1* gamma1 .* grad_square(phi)/2 +  Lambda_2*gamma2.*grad_square(rho)/2 + E_no_linear(phi,rho));
end

function L = augmented_lagrangian(alpha, lambda, kappa, u_new1, u_new2, u0, v_new1, v_new2, v0, ...
                                 phi_new1, phi_new2, phi0, rho_new1, rho_new2, rho0, ...
                                 mu_phi_star, mu_rho_star, u_star, v_star, rhs1, rhs2, rhs3, rhs4)
    % 计算约束 N(alpha)
    N = non_fun_alpha(alpha, u_new1, u_new2, u0, v_new1, v_new2, v0, ...
                      phi_new1, phi_new2, phi0, rho_new1, rho_new2, rho0, ...
                      mu_phi_star, mu_rho_star, u_star, v_star, rhs1, rhs2, rhs3, rhs4);
    % 增广拉格朗日函数
    L = alpha^2 + lambda * N + (kappa / 2) * N^2;
end

function r = non_fun_alpha(alpha,u_new1,u_new2,u0,v_new1,v_new2,v0,phi_new1,phi_new2,phi0,rho_new1,rho_new2,rho0,mu_phi_star,mu_rho_star,u_star,v_star,rhs1,rhs2,rhs3,rhs4)
global dt Mphi Mrho upsilon
    left = Energy(u_new1+ alpha.*u_new2, v_new1 + alpha.*v_new2, phi_new1 + alpha.*phi_new2, rho_new1+alpha.*rho_new2)-Energy(u0,v0,phi0,rho0);
    right = fun_inner(mu_phi_star,Mphi.*lap_diff(mu_phi_star)) + fun_inner(rhs1,mu_phi_star)...
           +fun_inner(mu_rho_star,Mrho.*lap_diff(mu_rho_star)) + fun_inner(rhs2,mu_rho_star)...
           +fun_inner(u_star,upsilon.*lap_diff(u_star)) + fun_inner(rhs3,u_star)...
           +fun_inner(v_star,upsilon.*lap_diff(v_star)) + fun_inner(rhs4,v_star);
    r = left - dt*right;
end

function [] = calculate_energy(out1,out2,hx,hy,t,u,v,phi,rho)
energy0 = Energy(u,v,phi,rho);
massphi   = fft2(phi);
massphi    = hx*hy*massphi(1,1);
massrho    = fft2(rho);
massrho    = hx*hy*massrho(1,1);
fprintf(out1,'%14.6e  %f   %f\n',t,massphi,massrho);
fprintf(out2,'%14.6e  %f   %f\n',t,energy0,energy0);
end


