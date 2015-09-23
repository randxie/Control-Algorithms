function [xout,uout]=Controller(sys,para)
%Controller designed for Nonlinear control homework
%   para: contain h,StartTime,EndTime,TimeVec,NumOfStates,CtlType,InitState
%   author: Rand Xie
%   date:   2015 Apr 18
global PFfun %prefilter values

%% Predefined Parameters
CtlType=para.CtlType;
StartTime=para.StartTime;
h=para.h;
EndTime=para.EndTime;
TimeVec=para.TimeVec;

%% Select controller
switch CtlType
%     case 'Test'
%         x0=sys.x0;
%         xout=RK4(TimeVec,x0,h,@DynEqn);
%         xout=xout';
    case 'FBL'
        PFfun=Prefilter(TimeVec,sys);
        xout=CalState(TimeVec,sys,para);
        [zeta,uout]=CalInnerState(TimeVec,xout);
    case 'SMC'
        PFfun=Prefilter(TimeVec,sys);
        xout=CalState(TimeVec,sys,para);
        [zeta,uout]=CalInnerState(TimeVec,xout);
end
end

%% Calculate zetas and inputs
function [zeta,u]=CalInnerState(TimeVec,x)
global sys zetaTrans HasMeasurementNoise
HasMeasurementNoise=false; %stop generating noise when calculating inner states
zeta=zeros(length(TimeVec),sys.NumOfStates);
u=zeros(length(TimeVec),sys.NumOfInputs);
for i=1:length(TimeVec)
    r=x(i,1);     dr=x(i,2);        s=x(i,3);     ds=x(i,4);
    theta=x(i,5); dtheta=x(i,6);    phi=x(i,7);   dphi=x(i,8);
    zetatmp=zetaTrans(r,dr,s,ds,theta,dtheta,phi,dphi);
    zeta(i,:)=zetatmp';
    u(i,:)=CalInput(TimeVec(i),zetatmp,x(i,:));
end
end

%% Caluclate state changes using RK4
function [xout]=CalState(TimeVec,sys,para)
x0=sys.x0;
if(strcmp(para.CtlType,'FBL'))
    xout=RK4(TimeVec,x0,para.h,@FBL);
    xout=xout';
elseif (strcmp(para.CtlType,'SMC'))
    xout=RK4(TimeVec,x0,para.h,@SMC);
    xout=xout';
elseif (strcmp(para.CtlType,'Optimal'))
    xout=RK4(TimeVec,x0,para.h,@Optimal);
    xout=xout';
end
end

%% Feedback Linearized controller
function [xdot]=FBL(t,x)
global  HasMeasurementNoise
if(HasMeasurementNoise)
    x=x+CalMeasurementNoise();
end
r=x(1);     dr=x(2);        s=x(3);     ds=x(4);
theta=x(5); dtheta=x(6);    phi=x(7);   dphi=x(8);
global zetaTrans
zeta=zetaTrans(r,dr,s,ds,theta,dtheta,phi,dphi);
u=CalInput(t,zeta,x);
xdot=DynEqn(x,u);
end

%% Sliding mode controller
function [xdot]=SMC(t,x)
global  HasMeasurementNoise
if(HasMeasurementNoise)
    x=x+CalMeasurementNoise();
end
r=x(1);     dr=x(2);        s=x(3);     ds=x(4);
theta=x(5); dtheta=x(6);    phi=x(7);   dphi=x(8);
global zetaTrans
zeta=zetaTrans(r,dr,s,ds,theta,dtheta,phi,dphi);
u=CalInput(t,zeta,x);
xdot=DynEqn(x,u);
end

%% Noise Generator for measurement
function [CalNoise]=CalMeasurementNoise()
CalNoise=randn(8,1)/10;
end

%% Calculate input
function [uout]=CalInput(t,zeta,x)
global sys PFfun para Dr Cr
r=x(1);     dr=x(2);        s=x(3);     ds=x(4);
theta=x(5); dtheta=x(6);    phi=x(7);   dphi=x(8);
yd=PFfun(t);
if(strcmp(para.CtlType,'FBL'))
    uout=(Dr(r,dr,s,ds,theta,dtheta,phi,dphi))\...
        (-[sys.LrFBL,zeros(1,4);zeros(1,4),sys.LsFBL]*(zeta-yd)-Cr(r,dr,s,ds,theta,dtheta,phi,dphi));
elseif (strcmp(para.CtlType,'SMC')) %hybrid controller when Dr is badly scaled
    us=-sys.ita*sign([sys.LrSMC,zeros(1,4);zeros(1,4),sys.LsSMC]*(zeta-yd));
    if(isnan(abs(det(Dr(r,dr,s,ds,theta,dtheta,phi,dphi)))))
        uout=(us-Cr(r,dr,s,ds,theta,dtheta,phi,dphi));
    else
        uout=(Dr(r,dr,s,ds,theta,dtheta,phi,dphi))\(us-Cr(r,dr,s,ds,theta,dtheta,phi,dphi));
    end
end
end

%% Prefilter for FBL and SMC
function [PFfun]=Prefilter(TimeVec,sys)
global PF
pf=sys.pf;
Ayr=pf.Ayr;  Byr=pf.Byr;    Kyr=pf.Kyr;
Pyr=pf.Pyr;  yT0=pf.yT0;
Aarg=[(Ayr-Byr*Kyr),zeros(4,4);zeros(4,4),(Ayr-Byr*Kyr)];
Barg=[Byr*Pyr zeros(4,1);zeros(4,1) Byr*Pyr];
PreF=@(t,x) Aarg*x+Barg*(Map(t));
[tout,PF]=ode45(PreF,TimeVec,yT0);
h=TimeVec(2)-TimeVec(1);
PFfun=@(t) (PF(1+floor(t/h),:)+(PF(1+ceil(t/h),:)-PF(1+floor(t/h),:))*(t/h-floor(t/h)))';
end

%% Dynamic equation for the system
function [xdot]=DynEqn(x,u)
global fx gx fxHat gxHat SysUncertainty HasStateNoise
r=x(1);     dr=x(2);        s=x(3);     ds=x(4);
theta=x(5); dtheta=x(6);    phi=x(7);   dphi=x(8);
if(~SysUncertainty)
    fnum=fx(r,dr,s,ds,theta,dtheta,phi,dphi);
    gnum=gx(r,dr,s,ds,theta,dtheta,phi,dphi);
else
    fnum=fxHat(r,dr,s,ds,theta,dtheta,phi,dphi);
    gnum=gxHat(r,dr,s,ds,theta,dtheta,phi,dphi);
end
xdot=(fnum)+(gnum)*u;
if(HasStateNoise)
    xdot=xdot+CalStateNoise();
end
end

function [CalNoise]=CalStateNoise()
CalNoise=randn(8,1)/10;
end