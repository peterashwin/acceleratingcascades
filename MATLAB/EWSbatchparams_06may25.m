% EWSbatch.m

%% Default parameters
%% parameter shift from lambdam to lambdap
par.lambdam=1.0; %% lower limit lambda
par.lambdap=3.0; %% upper limit lambda
%par.lambdap=1.0; %% upper limit lambda
%% add lambdah*bump(a*s+b)
par.a=1.0;  %
par.b=0.0; %
par.lambdah=0.0; %
%% default rate
par.r=0.025; % rate

%% choose model
% =1 for 1D fold, =2 for 1D fold plus stable, =3 for fifth order
par.model=1;
% model params
par.model_xi12=0;
par.model_xi21=0;
par.model_xi22=0;
% Initial condition
par.x0=-1.7;

%% noise
par.eta=0.05;

%% Simulation Setup
%
par.Name = 'EWS default';
par.StartTime = -105;
par.ShowStartTime = -100;
par.EndTime = 20;
par.ShowEndTime = 20;
par.dt = 0.01;

%% Integrator parameters (Heun)
par.sde_opts.dt=par.dt;
par.sde_opts.thin=5;
par.sde_opts.eta=par.eta;

%% Test project forward time
par.StartFit = [-60 -30];
par.StopFit = [-40 -10];

%% repeats
par.repeats=1;

%% EWS parameters
par.nw=4000; % estimator window length
par.alphasqmax=30; % largest alpha^2 to plot

%% RNG seed
par.seed=2024;

%% default constants for coupling function
par.c0=0;
par.c1=0;
par.c2=2.5;
par.c3=2;
par.c4=0.5;

%% default number
pno=0;
par.no=pno;

%% default
pno=pno+1;
pars(pno)=par;
pars(pno).no=pno;
pars(pno).repeats=5;

%% bump term on earlier start to fit
pno=pno+1;
pars(pno)=par;
pars(pno).no=pno;
pars(pno).a=1.0;
pars(pno).b=0;
pars(pno).lambdah=0.9;
pars(pno).lambdap=1.0;
pars(pno).repeats=5;

%% fold+stable earlier fit
pno=pno+1;
par.no=pno;
pars(pno)=par;
pars(pno).no=pno;
pars(pno).model=2;
% model params
pars(pno).model_xi12=0.2;
pars(pno).model_xi21=0.2;
pars(pno).model_xi22=-sqrt(5);
% Initial condition
pars(pno).x0=[-1.7; 1];
pars(pno).StartFit=[-50 -15];
pars(pno).StopFit=[-30 5];
pars(pno).repeats=5;


% %% fifth order earlier fit
% pno=pno+1;
% par.no=pno;
% pars(pno)=par;
% pars(pno).no=pno;
% pars(pno).model=3;
% % Initial condition
% pars(pno).x0=1.7;
% pars(pno).r=0.002;
% pars(pno).eta=0.001;
% pars(pno).StartTime = -205;
% pars(pno).ShowStartTime = -200;
% pars(pno).StartFit=-60;
% pars(pno).StopFit=-40;
% pars(pno).ShowEndTime = 300;
% pars(pno).EndTime = 300;
% pars(pno).alphasqmax=20;
% pars(pno).repeats=5;
% 
% %% fifth order later fit
% pno=pno+1;
% par.no=pno;
% pars(pno)=par;
% pars(pno).no=pno;
% pars(pno).model=3;
% pars(pno).r=0.002;
% pars(pno).eta=0.001;
% pars(pno).StartTime = -205;
% pars(pno).ShowStartTime = -200;
% pars(pno).x0=1.7;
% pars(pno).StartFit=50;
% pars(pno).StopFit=70;
% pars(pno).ShowEndTime = 300;
% pars(pno).EndTime = 300;
% pars(pno).alphasqmax=20;
% pars(pno).repeats=5;
% % 


