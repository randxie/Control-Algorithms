%Calculate Norm Form from initial states
%f is n*1,  h is p*1   g is n*m
function zeta= CalLie(f,h,g,x)
zeta = sym([]);
p = length(h);
m = size(g,2);
for i=1: p
    j=1;
    lfh=h(i);
     while(j~=0)
        for k =1:m
            lied=liederivative(g(:,k),lfh,x,1);
            l_g_j_hi(k) = lied(2);
        end
         if nnz(l_g_j_hi)==0 
            j=j+1; 
            lfhtem = liederivative(f,lfh,x,1); 
            lfh = lfhtem(2);
        else 
             order = j; 
             j=0;
        end
     end
     l = length(zeta);
     zeta((l+1):(order+l)) = liederivative(f,h(i),x,order-1); 
    
end
zeta = zeta'; 
end      
        
        


