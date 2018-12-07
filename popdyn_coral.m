function dy = popdyn_coral(t,y,pr)
% Author: Karlo Hock, University of Queensland.
% Equations from Fabina et al. 2015, Ecological Applications, 25(6), 1534–1545

%Coral-algal ODEs for a given set of parameter values

a = pr(1);
s = pr(2);
n = pr(3);
g = pr(4);
b = pr(5);
m = pr(6);
h = pr(7);
z = pr(8);
o = pr(9);
r = pr(10);
C = y(1);
M = y(2);


dy = zeros(2,1);
R = -(C*r*(C + M - 1))/(a + n + C*r + M*s);%pre-solved coral recruitment equation dR = r*A*(1-R-A-M) - a*R - s*R*M - n*R to obtain R
dy(1) = a*R + g*C*(1-R-C-M) - b*C*M - m*C;
dy(2) = s*M*(1-R-C-M) + s*R*M + b*C*M - h*M - (z*M*o*C)/(1+o*C);


