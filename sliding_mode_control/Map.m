% Define any curve you like here!
function [yout]=Map(t)
global DrawMap MapPoints MaxTime
if(~DrawMap)
    if(t<5)
        tmp=t;
        yout=[0.05;tmp/100];
    elseif (t>=5)&(t<10)
        tmp=t-5;
        yout=[0.05-tmp/100;0.05];
    elseif (t>=10) & (t<15)
        tmp=t-10;
        yout=[0;0.05-tmp/100];
    else
        yout=[0;0];
    end
%     yout=[0.05*cos(t);0.05*sin(t)];
else
    t=(t/MaxTime*(length(MapPoints)-1));
    yout=(MapPoints(1+floor(t),:)+...
            (MapPoints(1+ceil(t),:)-MapPoints(1+floor(t),:))*(t-floor(t)))';
end

end