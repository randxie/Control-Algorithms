%% Test OptCtrSolver
clc;
clear;
para.StartTime=0;
para.h=0.01;
para.EndTime=20;
para.TimeVec=(para.StartTime:para.h:para.EndTime)';
para.NumOfStates=2;
para.CtlType='LQT';
para.InitState=[0;0];

sys.A=[0 1;0 0];
sys.B=[0;1];
sys.C=[1 0;0 sqrt(10)];
sys.D=0;
sys.Q=[1 0;0 1];    
sys.R=1;
sys.r=@(t) [1;0];  %control x1 to certain point and x2 to 0
sys.P=0; 

[Time,States]=OptCtrSolver(sys,para);
uout=zeros(length(Time),1);
for uiter=1:length(Time)
    uout(uiter)=u(Time(uiter),States(uiter,:)');
end


subplot(3,1,1);
plot(Time,States(:,1));
% axis([0 Time(end) 0 2]);
xlabel('Time(s)');
ylabel('x1');

subplot(3,1,2);
plot(Time,States(:,2));
% axis([0 Time(end) 0 0.5]);
xlabel('Time(s)');
ylabel('x2');

subplot(3,1,3);
plot(Time,uout);
%axis([0 Time(end) 0 0.5]);
xlabel('Time(s)');
ylabel('Input u(t)');
