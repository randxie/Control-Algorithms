% NL Project Feedback Linearization (normal form) 
% original states [r dr s ds theta dtheta phi dphi]'; 
% to norm form
function [zetaTrans,fx,gx,Cr,Dr]=CalNormForm(SysUncertainty)
syms r dr s ds theta dtheta phi dphi real 
syms b m gty Jpx Jpy Jpz B Jf  real

f = [dr;
    -5*b/(7*m)*dr-5/7*gty*sin(theta);
    ds;
    -5*b/(7*m)*ds+5/7*gty*cos(theta)*sin(phi);
    dtheta;
    1/(Jpy*cos(phi)^2+Jpz*sin(phi)^2)*(-B*dtheta-m*gty*(s*sin(theta)*sin(phi)+r*cos(theta)));
    dphi;
    1/(Jf+Jpx)*(-B*dphi+m*gty*s*cos(theta)*cos(phi))];

g = [0 0 0 0 0 1/(Jpy*cos(phi)^2+Jpz*sin(phi)^2)  0  0;
     0 0 0 0 0 0                                  0  1/(Jf+Jpx)]'; 
 
h = [r, s]';
x = [r dr s ds theta dtheta phi dphi]';
zeta = CalLie(f,h,g,x);

%% load system data
bc=0.01;
mc=0.00427; %kg
gtyc=9.8;
Jpxc=0.00108;
Jpyc=0.00108;
Jpzc=0.00216;
Bc=0.2;
Jfc=0.033;
if(SysUncertainty)
%     bc=bc-10/100*bc;
    bc=0.05;
    Jfc=Jfc-Jfc*0.5;
end

zetaTrans=subs(zeta,{'b','m','gty','Jpx','Jpy','Jpz','B','Jf'},...
                {bc, mc, gtyc, Jpxc, Jpyc, Jpzc, Bc, Jfc});
zetaTrans=matlabFunction(zetaTrans, 'vars', [r dr s ds theta dtheta phi dphi]);           
fx=subs(f,{'b','m','gty','Jpx','Jpy','Jpz','B','Jf'},...
                {bc, mc, gtyc, Jpxc, Jpyc, Jpzc, Bc, Jfc});
            
fx=matlabFunction(fx, 'vars', [r dr s ds theta dtheta phi dphi]);

gx=subs(g,{'b','m','gty','Jpx','Jpy','Jpz','B','Jf'},...
                {bc, mc, gtyc, Jpxc, Jpyc, Jpzc, Bc, Jfc});
            
gx=matlabFunction(gx, 'vars', [r dr s ds theta dtheta phi dphi]);

lf4h1tem = liederivative(f,h(1),x,4);
lf4h1 = lf4h1tem(5);
lf3h1 = lf4h1tem(4);
lf4h2tem = liederivative(f,h(2),x,4);
lf4h2 = lf4h2tem(5);
lf3h2 = lf4h2tem(4);
Cr = [lf4h1; lf4h2];
Cr=subs(Cr,{'b','m','gty','Jpx','Jpy','Jpz','B','Jf'},...
                {bc, mc, gtyc, Jpxc, Jpyc, Jpzc, Bc, Jfc});
            
Cr=matlabFunction(Cr, 'vars', [r dr s ds theta dtheta phi dphi]);    

lg1lf3h1tem =  liederivative(g(:,1),lf3h1,x,1);
lg1lf3h1 = lg1lf3h1tem(2);
lg2lf3h1tem =  liederivative(g(:,2),lf3h1,x,1);
lg2lf3h1 = lg2lf3h1tem(2);
lg1lf3h2tem =  liederivative(g(:,1),lf3h2,x,1);
lg1lf3h2 = lg1lf3h2tem(2);
lg2lf3h2tem =  liederivative(g(:,2),lf3h2,x,1);
lg2lf3h2 = lg2lf3h2tem(2);
Dr = [lg1lf3h1,lg2lf3h1;lg1lf3h2, lg2lf3h2];
Dr=subs(Dr,{'b','m','gty','Jpx','Jpy','Jpz','B','Jf'},...
                {bc, mc, gtyc, Jpxc, Jpyc, Jpzc, Bc, Jfc});

Dr=matlabFunction(Dr, 'vars', [r dr s ds theta dtheta phi dphi]);   
end
