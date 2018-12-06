%coral_Fabina ODE's

function dy = popdyn_fabina_2015(t,y,pr)


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
R = y(1);
A = y(2);
M = y(3);


dy = zeros(3,1);
dy(1) = r*A*(1-R-A-M) - a*R - s*R*M - n*R;
dy(2) = a*R + g*A*(1-R-A-M) - b*A*M - m*A;
dy(3) = s*M*(1-R-A-M) + s*R*M + b*A*M - h*M - (z*M*o*A)/(1+o*A);


