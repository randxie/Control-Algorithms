function [y_res]=RK4(t,y0,h,f)
y_res=zeros(numel(y0),numel(t));
y_res(:,1)=y0;
for n=2:numel(t)
y_res(:,n)=RK4_i(t(n-1),y_res(:,n-1),h,f) ;
end
end

function [yn_new]=RK4_i(tn,yn,h,f)
k1=f(tn,yn) ;
k2=f(tn+h/2,yn+1/2*k1*h) ;
k3=f(tn+h/2,yn+1/2*k2*h) ;
k4=f(tn+h,yn+k3*h) ;
yn_new=(yn+h/6*(k1+2*k2+2*k3+k4));
end