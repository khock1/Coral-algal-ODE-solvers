% Author: Karlo Hock, University of Queensland.
% Equations from Fabina et al. 2015, Ecological Applications, 25(6), 1534–1545

%Perfom stability analysis for coral-algal ODEs over a range of parameter values
%Test a range of grazing pressures to examine how stability structure changes
%Plot hysteresis plot with different colours for stable vs unstable equilibria

clear;
syms C M;
% parameter values ---------
r = 0.05;%coral recruitment
a = 0.2;%recruit maturation
g = 0.1;%adult growth
n = 0.8;%recruit mortality
m = 0.03;%adult mortality
s = 0.4;%macoralgal recruitment/recruit overgrowth
h = 0.2;%baseline macroalgal mortality
z = 0.4;%supplemental macroalgal mortality, from herbivores
b = 0.4;%adult overgrowth
o = 4;%herbivore habitat provisioning

parameter_range= 0:0.02:0.4;

figure;hold
for i=1:length(parameter_range)
    z = parameter_range(i);%subsititute z here to perform stability analysis with a different paramter
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
                    scatter(parameter_range(i),ss1.C(jj),30,'filled','MarkerFaceColor','m');%if unstable
                else
                    scatter(parameter_range(i),ss1.C(jj),30,'filled','MarkerFaceColor','r');%if stable
                end
            end
        end
    end
end
ylabel('Adult Coral');
xlabel('External Supply');
