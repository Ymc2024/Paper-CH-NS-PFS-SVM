%% test order of convergence in time

close all;
clear; clc;

% add path
addpath('../','-begin');

% Parameters
para.Mphi =1e-4;
para.Mrho = 1e-4;
para.epsilon = 0.04;
para.eta = 0.01;
para.gamma1 = 1;
para.gamma2 = 1;
para.gamma = 1;
para.theta = 0.2;
para.theta1 = 0.02;
para.theta2 = 0.02;
para.xi = 1e-3;
para.S1=10;
para.S2=10;
para.Lambda_1=1;
para.Lambda_2=1;
para.alpha=1;
para.alpha1=2;
para.alpha2=-2;
para.beta = 4;
para.Sigma=0.001;
para.upsilon=1;
para.rhos=1;
para.C0=1e10;

% Time: dt T
T = 50;
t0 = 0;
% tsave = 10 * T;
Tsplit = [5 10 T];
tsave = [0.1 1 5];

%% pde
dt_ref = 1e-4;
 pde = ex05_1_SFdata(para);
%pde = ex05_2_SFdata(para);
%% Space: Domain and N
domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;
Nx = 128; Ny = Nx;

% scheme = 's1st';   % First-order scheme
% scheme = 's2cn';   % Second-order CN scheme, no need to write
scheme = 's2bdf';  % Second-order BDF scheme

%% option
option.scheme = scheme;
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 0;
option.energyflag = 0;
option.tol = 1e-14;
option.tolit = 1e-12;
option.maxit = 2000;
option.type = 'SVM1';

%% Run:
time = struct('T',T,'t0',t0,'dt',dt_ref,'Tsplit',Tsplit,'tsave',tsave);
if 1 == strcmp(scheme,'s1st')
    Surfactant_2D_NS_SVM_1st(pde,domain,Nx,Ny,time,option);
    % Surfactant_2D_NS_SVM_1st_modified(pde,domain,Nx,Ny,time,option);
elseif 1 == strcmp(scheme,'s2cn')
    Surfactant_2D_NS_SVM_2cn(pde,domain,Nx,Ny,time,option);
    % Surfactant_2D_NS_SVM_2cn_modified(pde,domain,Nx,Ny,time,option);
elseif 1 == strcmp(scheme,'s2bdf')
    Surfactant_2D_NS_SVM_2bdf(pde,domain,Nx,Ny,time,option);
     % Surfactant_2D_NS_SVM_2bdf_modified(pde,domain,Nx,Ny,time,option);
end






