function [stable_eq, unstable_eq] = stability_coral(pr)
% Author: Karlo Hock, University of Queensland.
% Equations from Fabina et al. 2015, Ecological Applications, 25(6), 1534–1545

%Perfom stability analysis for coral-algal ODEs for a given set of parameter values

syms C M;
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

stable_eq=[];%container for stable equilibria values
unstable_eq=[];%container for unstable equilibria values
R = -(C*r*(C + M - 1))/(a + n + C*r + M*s);%pre-solved coral recruitment equation dR = r*A*(1-R-A-M) - a*R - s*R*M - n*R to obtain R
dA = a*R + g*C*(1-R-C-M) - b*C*M - m*C;%adult coral
dM = s*M*(1-R-C-M) + s*R*M + b*C*M - h*M - (z*M*o*C)/(1+o*C);%macroalgae
J = jacobian([dA, dM ],[C M]);%obtain Jacobian matrix for the ODE system
ss1 = vpasolve([0 == a*R + g*C*(1-R-C-M) - b*C*M - m*C,0==s*M*(1-R-C-M) + s*R*M + b*C*M - h*M - (z*M*o*C)/(1+o*C)],[C,M]);%obtain all numerical solutions for the ODEs
for jj=1:length(ss1.C)%sequentially plot equilibria for a given range of parameter values
    if isreal(ss1.C(jj)) && isreal(ss1.M(jj))%only concerned with real equilibria
        if all([ss1.C(jj) ss1.M(jj)]>=0 & [ss1.C(jj) ss1.M(jj)]<=1)%only concerned with nontrivial solutions, i.e. if proportion od coral and macroalgal cover is less than 1
            jaceigv=eig(double(subs(J,[C M],[ss1.C(jj) ss1.M(jj)])));%substitute equilibrai into Jacobian to obtain eigenvalues
            if any(jaceigv>=0)%eigenvalues >0 means equlilibrium is a saddle
                unstable_eq=double(vertcat(unstable_eq,[ss1.M(jj) ss1.C(jj)]));%if unstable
            else
                stable_eq=double(vertcat(stable_eq,[ss1.M(jj) ss1.C(jj)]));%if stable
            end
        end
    end
end


end