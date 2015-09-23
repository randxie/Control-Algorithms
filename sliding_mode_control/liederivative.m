function l_fnh = liederivative(f,h,x,n)
% LIEDERIVATIVE  (Nth order iterated Lie derivative of f and h)
%
% Usage:
% l_fnh = liederivative(f,h,x,n)
%
% Input:
% f = symbolic vector field of length l1 (first operand)
% h = symbolic vector field of length 1 (second operand)
% x = symbolic vector of variables of length l1
% n = order of the Lie Derivative
%
% Output:
% l_fnh =  [ h    l_fh  l_f^2(h)   l_f^3(h).... ]
%        n=   0      1         2         3
%             
%
%
% Note: define the symbolic variables as 
%       x1 =sym(x1,'real')
%       in order to get rid of the conjugates
%       (only if they are reals).
%
% Author: Atakan Varol
% Date: 03.23.2006

%l_fnh = sym(zeros(n+1,1));
l_fnh =sym([]);
l_fnh(1) = h;

if n>0
    for t = 2:n+1
        l_fnh(t) = jacobian(l_fnh(t-1),x)*f;
    end
end

l_fnh = expand(l_fnh);

% End of code