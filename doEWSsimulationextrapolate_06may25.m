%% doEWSsimulationextrapolate
% This script runs a numerical simulation of fold tipping to illustrate EWS
% methodology, including extrapolation
%   Peter Ashwin
%   Apr 2025
%
%% Start with a clean slate
close all
clear all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');


%% read in parameters to pars
run('EWSbatchparams_06may25.m')

%% run through parameter sets
for ip=1:length(pars)
    %for ip=7:8
    %% start timer
    tic;
    %% run through each of the parameter values
    par=pars(ip);
    %% set random seed
    rng(par.seed);
    dt=par.sde_opts.dt;

    %% choose model
    if par.model==1
        fun=@(t,y,par) fold_model(t,y,par);
        jac=@(t,y,par) jac_fold_model(t,y,par);
    elseif par.model==2
        fun=@(t,y,par) fold_model_stable(t,y,par);
        jac=@(t,y,par) jac_fold_model_stable(t,y,par);
    elseif par.model==3
        fun=@(t,y,par) fifth_model(t,y,par);
        jac=@(t,y,par) jac_fifth_model(t,y,par);
    else
        echo('PANIC: MODEL NOT KNOWN!!');
        keyboard;
    end

    %% Plot timeeries
    f1=figure;
    clf;
    f1.PaperUnits='centimeters';
    f1.Units='centimeters';

    f1.PaperSize=[10 12];
    f1.InnerPosition=[10 1 10 12];



    for rep=1:par.repeats
        sprintf('Params %d, run %d of %d',i,rep,par.repeats)

        %% timeseries fig
        figure(f1);

        %% Initial state vector
        par.xinit = [par.x0];

        %% Call the function that runs the actual numerical simulation
        [var] = EWSsimulator(fun,par);
        vars(rep)=var;

        tt=var.t;
        xx=var.x;

        % plot only first component if multidimensional!
        xs=xx(1,:);


        %% plot timeseries with tipping region


        p0=subplot(4,1,1);
        hold on
        lambdas=fLambda(tt,par);
        if rep==1
            plot(tt,lambdas,"Color","k");
            text(0.03,0.8,'(a)','Units','normalized');
        end
        if par.model==1 || par.model==2
            plot([par.ShowStartTime par.ShowEndTime],[2 2],'r:')
        elseif par.model==3
            plot([par.ShowStartTime par.ShowEndTime],[2.1066 2.1066],'k:')
            plot([par.ShowStartTime par.ShowEndTime],[2.471 2.471],'r:')
        else
            echo('PANIC: MODEL NOT KNOWN!!');
            keyboard;
        end

        xlim([par.ShowStartTime par.ShowEndTime]);
        ylabel('$\Lambda(rt)$');

        %% plot timeseries with tipping region
        p1=subplot(4,1,2);
        hold on

        if rep==1
            plot(tt,xs,"Color",[1 0 0]);
            text(0.03,0.8,'(b)','Units','normalized');
        else
            plot(tt,xs,"Color",[1 0 0 0.2]);
        end
        xlim([par.ShowStartTime par.ShowEndTime]);
        ylim([-2.2 2.2]);
        %legend('$x(t)$');%,'$x^*(t)$');
        ylabel('$x$');
        p1.XTickLabel=[];


        %% perform analysis of (possibly truncated) EWS
        tts=tt(tt<=par.ShowEndTime);
        xxs=xx(:,tt<=par.ShowEndTime);

        %% EWS
        [xSD,xAR1]=EWSindicators(tts,xxs,par,jac);

        %% AR1 estimator variance (valid if \alpha dt is small)
        nwin=par.nw;
        xsigma = sqrt(2*xAR1*dt/nwin);

        %%
        %p2=subplot(4,1,3);
        %hold on;
        timewin= par.ShowStartTime+par.nw*par.dt;
        pgon=polyshape([par.ShowStartTime par.ShowStartTime timewin timewin],[0 1 1 0]);

        %%

        p3=subplot(4,1,3);

        hold on;
        if rep==1
            % show estimator window
            plot(pgon,"FaceColor",[1 0 0],"FaceAlpha",0.1,"EdgeAlpha",0.0);
            % plot plus and minus 2 estimated SD
            plot(tts,xAR1,"Color",[1 0 0]);
            plot(tts,xAR1+2*xsigma,':');
            plot(tts,xAR1-2*xsigma,':');
            text(0.03,0.8,'(c)','Units','normalized');
        else
            plot(tts,xAR1,"Color",[1 0 0 0.1]);
        end

        xlim([par.ShowStartTime par.ShowEndTime]);
        ylim([min(xAR1)*0.9 1.005]);
        %xlabel('t [yr]')
        p3.XTickLabel=[];
        ylabel('$\hat{\rho}$');
        %
        hold off;

        p4=subplot(4,1,4);

        xalp= -log(xAR1)/dt;

        % plot plus and minus 2 estimated SD
        xalpminus = -log(xAR1-2*xsigma)/dt;
        xalpplus = -log(xAR1+2*xsigma)/dt;

        hold on;
        if rep==1
            plot(tts,xalp.^2,"Color",[1 0 0]);
            plot(tts,xalpplus.^2,':',"Color",[1 0 0]);
            plot(tts,xalpminus.^2,':',"Color",[1 0 0]);

            % plot alpha using linearization at xs
            [jacxs,jacv]=jac(tt,xx,par);

            plot(tt,jacxs.^2,"Color",[0 0 0]);
        else
            plot(tts,xalp.^2,"Color",[1 0 0 0.1]);
        end

        %% fit within time window [par.StartFit,par.StopFit]


        % only fit the first repeat - others are there for ensemble
        if rep==1
            % show estimator window
            pgon=polyshape([par.ShowStartTime par.ShowStartTime timewin timewin],[0 par.alphasqmax par.alphasqmax 0]);
            plot(pgon,"FaceColor",[1 0 0],"FaceAlpha",0.1,"EdgeAlpha",0.0);
            nfits=length(par.StopFit);
            for k=1:nfits;
                thisStartFit=par.StartFit(k);
                thisStopFit=par.StopFit(k);
                ttsw=tts(and(tts>=thisStartFit,tts<=thisStopFit));
                xalpwsq=xalp(and(tts>=thisStartFit,tts<=thisStopFit)).^2;
                % show fitting window
                pgon2=polyshape([thisStartFit thisStartFit thisStopFit thisStopFit],[0 par.alphasqmax par.alphasqmax 0]);
                plot(pgon2,"FaceColor",[0 1 0],"FaceAlpha",0.2*(k/nfits),"EdgeAlpha",0.0);
                % fit according to model
                modelfun = @(b,x) b(1)+b(2)*x;
                beta0 = [5; 0];

                mdl = fitnlm(ttsw,xalpwsq,modelfun,beta0);
                % report model
                mdl

                % compute extrapoation
                betaest=mdl.Coefficients.Estimate;
                ttext=tt(and(tt>thisStartFit,tt<par.ShowEndTime));
                xalpsqext=modelfun(betaest,ttext);
                pp=plot(ttext,xalpsqext,"LineWidth",2,"Color",[0 1 0 0.5*(k/nfits)]);

                % find and mark first zero to within tolerance
                tol=0;
                ttip=find(xalpsqext<=tol,1);
                if(~isempty(ttip))
                    plot(ttext(ttip),0,"*","Color","green");
                end
            end
            text(0.03,0.8,'(d)','Units','normalized');
        end
        xlim([par.ShowStartTime par.ShowEndTime]);
        ylim([0 par.alphasqmax]);
        xlabel('$t$')
        ylabel('$\hat{\alpha}^2$');
        %
        hold off;
    end

    fontsize(gcf,scale=1.5);

    %% Save data, pars and figures
    spath = './Data/';
    fpath = './Figures/';

    name = par.Name;
    fname=sprintf('par-%i-%s',par.no,string(datetime("today")));

    sfile_name = sprintf('%s%s%s',spath,fname,'.m');
    save( sfile_name, 'par', 'var');

    % save timeseries figure
    ffile_name = sprintf('%s%s-1%s',fpath,fname,'.pdf');
    %exportgraphics(f1,ffile_name,'ContentType','vector');
    exportgraphics(f1,ffile_name,'ContentType','image',"Resolution",300);

end

%% Stop timer
toc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FUNCTIONS           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do simulation
function [var] = EWSsimulator(fun,par)

x0=par.x0;
ff=@(t,y) fun(t,y,par);


% compute transient
tspan = [par.StartTime par.ShowStartTime];
[t,x] = heun(ff , tspan, x0, par.sde_opts);

% compute trajectory from end of transient
tspan = [par.ShowStartTime par.EndTime];
x0=x(:,end);
[t,x] = heun(ff , tspan, x0, par.sde_opts);

var.t = t;
var.x = x;

end

%% 1D fold
function dxdt = fold_model(t, x, par)

temp=3*x-x.^3+fLambda(t,par);
dxdt=temp;

end

%
function [jac,jacv]=jac_fold_model(t, x, par)

temp=3-3*x.^2;

jac=temp;
jacv=1;

end

%% Fifth order fold
function dxdt = fifth_model(t, x, par)

temp=-(x+1).*x.*(x-1).*(x-1.4).*(x-1.8)+2-fLambda(t,par);
dxdt=temp;

end

%
function [jac,jacv]=jac_fifth_model(t, x, par)

temp=-5.*x.^4 + 12.8*x.^3 - 4.56*x.^2 - 6.40*x + 2.52;

jac=temp;
jacv=1;

end




%% 1D fold + stable
function dxxdt = fold_model_stable(t, xx, par)

x=xx(1,:);
y=xx(2,:);
lambda=fLambda(t,par);
temp=3*x-x.^3+par.model_xi12*y+lambda;
dxdt=temp;
dydt=par.model_xi21*x+par.model_xi22*y;
dxxdt=[dxdt;dydt];

end

%% return (sorted) eigenvalues and left vectors of jacobian
function [jac,jacev]=jac_fold_model_stable(t, xx, par)

x=xx(1,:);
y=xx(2,:);
lambda=fLambda(t,par);
ty=y+lambda;
for i=1:size(xx,2)
    jj=[3-3*x(i)^2, par.model_xi12 ; ...
        par.model_xi21, par.model_xi22];
    [tempv,tempe]=eigs(jj',2,'largestreal');
    % most positive eigenvalue first
    jac(1,i)=tempe(1,1);
    jac(2,i)=tempe(2,2);
    jacev(1,i,:)=tempv(1,:);
    jacev(2,i,:)=tempv(2,:);
end

end

%%
function temp=fLambda(t,par)
% NB quadratic dependence in tanh!

s=par.r*t;
temp=par.lambdam+(par.lambdap-par.lambdam)*(tanh(s)+1)/2+par.lambdah*bump(par.a*s+par.b);

end


%% coupling function (output)
function temp=fnM(x,par)

temp=par.c0+par.c1*x+par.c2*bump(par.c3*(x-par.c4));

end

%% Heun integrator
function [Tout,Yout]=heun(ode,tspan,y0,options)
%   ode must be a function of Y, t and params @ode(t,y,params)
%   ode must be a function handle (ie @ode)
%
%   here use additive noise, ie g(t,Y)=eta
%
%   params must contain:
%   params.dt timestep split into
%   param.thin substeps
%   params.eta
%
%   tspan=[Tstart Tend]

n=size(y0,1);
dt=options.dt; % dt for output
eta=options.eta;

%only store every 'thin' points
thin=options.thin;
% substep
sdt=dt/thin;
jthin=0;

t0=tspan(1);

Tout=t0;
Yout=y0;

tp=t0;
Yp=y0;

tn=t0;

while tn<tspan(end)
    %update t
    tn=tp+sdt;

    %compute noise
    Wn=sqrt(sdt)*randn(n,1);

    %RK2/heun steps
    m1=feval(ode,tp,Yp);
    Yt=Yp+m1*sdt+eta*Wn;
    m2=feval(ode,tn,Yt);
    Yn=Yp+0.5*(m1+m2)*sdt+eta*Wn;

    %update previous data
    tp=tn;
    Yp=Yn;

    jthin=jthin+1;

    if jthin==thin
        %add to end of array
        Tout=[Tout,tn];
        Yout=[Yout,Yn];
        jthin=0;
    end

end

end


%% Calculate SD and AR1 indicators
function [xSD,xAR1] = EWSindicators(tt,xx,par,jac)
%
np=length(tt);

%% moving average window length
nw=par.nw;

% moving SD
xSD=nan(size(tt));
% moving AR
xAR1=nan(size(tt));
% moving variance
xvar=nan(size(tt));

% delay for AR
if nw>=np
    sprinft('WINDOW LONGER THAN TIMESERIES!!');
    return;
end

projxx=zeros(size(tt));
for i=nw+1:np
    %% find leading ev at end of window
    [jacxs,jacxv] = jac(tt(i),xx(:,i),par);
    leadingev=reshape(jacxv(1,end,:),1,size(par.x0,1));
    %% project onto leading ev at end of window
    for j=i-nw:i
        %leadingev=[0 1];
        projxx(j)=leadingev*xx(:,j);
    end
    %% linear detrend
    tempdet=detrend(projxx(i-nw:i),1);
    temp=tempdet(2:end);
    temps=tempdet(1:end-1);
    % correlation coefficient
    cc=corrcoef(temps,temp);
    % var
    xvar(i)=sum(temp.^2);
    % SD
    xSD(i)=sqrt(xvar(i)/nw);
    % AR1
    % xAR1(i)=sum(temps.*temp)/xvar(i);
    xAR1(i)=cc(1,2);
end

end

function temp=bump(x)
xsq=(abs(x)<1).*x.^2;
temp=(abs(x)<1).*exp(xsq./(xsq-1));
end