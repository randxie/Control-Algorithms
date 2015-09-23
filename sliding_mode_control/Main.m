% Main.m
clc
clear
close all

%% Define global variables
global fx gx zetaTrans Cr Dr % function used to calculate f(x), g(x), T(x)
global para sys              % system struct 
global DrawMap MapPoints MaxTime 
% with model uncertainty
global zetaTransHat fxHat gxHat CrHat DrHat HasSysUncertainty
global HasStateNoise HasMeasurementNoise

%% Adjustable parameters
DrawMap=false;              %whether take input from user
Animation=true;             %whether animate trajectory
MaxTime=30;                 %end time
HasSysUncertainty=true;     %add system uncertainty or not
para.CtlType='FBL';         %control type, FBL or SMC
para.h=0.001;               %step size for RK4
HasStateNoise=true;         %whether add noise into states equation
HasMeasurementNoise=true;   %whether measurement noise exists

%% Get map, user can draw map by mouse
if(DrawMap)
    h1=figure;
    axis([-0.1 0.1 -0.1 0.1]);
    xlabel('r');
    ylabel('s');
    title('Draw anything you like');
    MapPoints = get_pencil_curve(h1);
    disp('Finish drawing, Starting calculation');
end

%% Symbolic calculation (those have "Hat" are with system uncertainty
[zetaTrans,fx,gx,Cr,Dr]=CalNormForm(false);
[zetaTransHat,fxHat,gxHat,CrHat,DrHat]=CalNormForm(HasSysUncertainty);

%% Define parameters
para.StartTime=0;
para.EndTime=MaxTime;
para.TimeVec=para.StartTime:para.h:para.EndTime;

%% Define constants and intial values
NumOfStates=8;
NumOfInputs=2;
% state order: r, r_dot, s, s_dot, theta, theta_dot, fi, fi_dot
MapInitPoint=Map(0);
x0=[MapInitPoint(1);0;MapInitPoint(2);0;0;0;0;0]; 

%% Prefilter for r and s has similar structure
% Prefilter structure: y is estimated states, r(t) is reference
% y_dot=(Ayr-Byr*Kyr)y+Byr*Pyr*r(t)
% 
Ayr=[0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     0 0 0 0];       
Byr=[0 0 0 1]';
Poleyr=[-12+12j,-12-12j,-10+10j,-10-10j];   %pole selected
Kyr=place(Ayr,Byr,Poleyr);                  %pole place ment

tmp=([1 0 0 0]*inv(Ayr-Byr*Kyr)*Byr);
Pyr=-1./(tmp);                              %constant Pyr

% Store into struct pf
pf.Ayr=Ayr;
pf.Byr=Byr;
pf.Kyr=Kyr;
pf.Pyr=Pyr;

% Calculate initial value for prefilter
r=x0(1);     dr=x0(2);        s=x0(3);     ds=x0(4);
theta=x0(5); dtheta=x0(6);    phi=x0(7);   dphi=x0(8);
yT0=zetaTrans(r,dr,s,ds,theta,dtheta,phi,dphi);
pf.yT0=yT0;

%% Design sliding surface for SMC
As=[0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        0 0 0 0]; 
Bs=[0 0 0 1]';
PoleSMC=[-15+15j,-15-15j,-5-5j,-5+5j];                %pole for SMC controller
PoleFBL=[-5+5j,-5-5j,-15-15j,-15+15j];      %pole for FBL controller
KsFBL=place(As,Bs,PoleFBL);
KsSMC=place(As,Bs,PoleSMC);
LrSMC=KsSMC/KsSMC(1);
LsSMC=KsSMC/KsSMC(1);
ita=20000*[1 0;0 1];

%% Define system structure
sys.pf=pf;
sys.NumOfStates=NumOfStates;
sys.NumOfInputs=NumOfInputs;
sys.x0=x0;
sys.zetaTrans=zetaTrans;
sys.LrFBL=KsFBL;
sys.LsFBL=KsFBL;
sys.LrSMC=LrSMC;
sys.LsSMC=LsSMC;
sys.ita=ita;

tic
[States,Input]=Controller(sys,para);
toc

%% For plotting
TimeVec=para.TimeVec;
ydesired=zeros(2,length(TimeVec));
for i=1:length(TimeVec)
    ydesired(:,i)=Map(para.TimeVec(i));
end

figure(1);
set(0,'DefaultAxesFontSize',12.5);
subplot(4,2,1);
plot(TimeVec,States(:,1));
xlabel('Time(s)');
ylabel('r');

subplot(4,2,2);
plot(TimeVec,States(:,2));
xlabel('Time(s)');
ylabel('$\dot{r}$','Interpreter','latex');

subplot(4,2,3);
plot(TimeVec,States(:,3));
xlabel('Time(s)');
ylabel('s','Interpreter','latex');

subplot(4,2,4);
plot(TimeVec,States(:,4));
xlabel('Time(s)');
ylabel('$\dot{s}$','Interpreter','latex');

subplot(4,2,5);
plot(TimeVec,States(:,5));
xlabel('Time(s)');
ylabel('${\theta}$','Interpreter','latex');

subplot(4,2,6);
plot(TimeVec,States(:,6));
xlabel('Time(s)');
ylabel('$\dot{\theta}$','Interpreter','latex');

subplot(4,2,7);
plot(TimeVec,States(:,7));
xlabel('Time(s)');
ylabel('${\phi}$','Interpreter','latex');

subplot(4,2,8);
plot(TimeVec,States(:,8));
xlabel('Time(s)');
ylabel('$\dot{\phi}$','Interpreter','latex');

figure(2);
subplot(2,1,1);
plot(TimeVec,Input(:,1),'b','LineWidth',0.5);
ylabel('Input $u_1$','Interpreter','latex');
subplot(2,1,2);
plot(TimeVec,Input(:,2),'b','LineWidth',0.5);
xlabel('Time(s)');
ylabel('Input $u_2$','Interpreter','latex');

global PF
if(~Animation)
    figure(3);
    plot(States(1,1),States(1,3),'rp');
    hold on;
    plot(States(:,1),States(:,3),'r--','LineWidth',2);
    hold on;
    plot(ydesired(1,:),ydesired(2,:),'g','LineWidth',1);
    hold on;
    % the prefileter states are y1,y1dot,y1ddot,y1dddot,y2....
    plot(PF(:,1),PF(:,5),'b--','LineWidth',1);
    hold on;
    plot(States(1,1),States(1,3),'ro');
    xlabel('r');
    ylabel('s','Interpreter','latex');
    title([para.CtlType 'Controller']);
    axis equal

else
    figure(3);
    plot(ydesired(1,:),ydesired(2,:),'b--','LineWidth',0.5);
    hold on;
    ha = animatedline('Color','r','LineWidth',1);
    axis([-0.1 0.1 -0.1 0.1])
    for k = 1:15:length(TimeVec)
        addpoints(ha,States(k,1),States(k,3));
        drawnow
    end
    axis equal;
end