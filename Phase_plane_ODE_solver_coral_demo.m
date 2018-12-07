% Author: Karlo Hock, University of Queensland.
% Equations from Fabina et al. 2015, Ecological Applications, 25(6), 1534�1545

% Solve ODEs for coral-macroalgal dynamics, plot pahse plane portrait of the solutions
% Solve for nullclines, then plot nullclines

clear;
%constrain parameter space to realistic values, ensuring that covers together cannot exceed 1
x=0:0.05:1;
y=0:0.05:1;
spc=[0 0];
axlim=0.8;%limit to parameter space to be xexplored, has to be <=1
for i=1:size(x,2)
    for j=1:size(y,2)
        if x(i)+y(j)<=axlim
            spc=vertcat(spc,[x(i) y(j)]);
        end
    end
end
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

parameters = [a,s,n,g,b,m,h,z,o,r];
figure;
hold on;

%solve and draw solution for parameter space
for ii = 1:length(spc)
    options = odeset('NonNegative',1);
    [T,Y] = ode45(@popdyn_coral,[0 200],[spc(ii,2),spc(ii,1)],options, parameters );%solve ODEs for parameter space
    C = Y(:,1); % adult coral
    M = Y(:,2); % macroalgae
    phaseplane = plot(M,C,'Color','k'); % plot solutions
end

title('Phase plane - Coral-macroalgal dynamics', 'FontSize',11)
xlabel('Macroalgae');
ylabel('Adult Coral');
%perform stability analysis, draw stable and unstable equilbiria on a phase plane
[stable_eq, unstable_eq] = stability_coral(parameters);
for jj=1:size(unstable_eq,1)
    scatter(unstable_eq(jj,1),unstable_eq(jj,2),140,'d','m','filled');
end
for jj=1:size(stable_eq,1)
    scatter(stable_eq(jj,1),stable_eq(jj,2),140,'d','c','filled');%if stable
end
ylabel('Adult Coral');
xlabel('Macroalgae');
axis([0 axlim 0 axlim]);
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 coral_phaseplane.eps

%solve for coral and algal nullclines, draw them and draw the axes
figure; hold;
syms Cs Ms;
Rs = -(Cs*r*(Cs + Ms - 1))/(a + n + Cs*r + Ms*s);%pre-solved coral recruitment equation dR = r*A*(1-R-A-M) - a*R - s*R*M - n*R to obtain R
coral_nullcline=solve(0 == a.*Rs + g.*Cs.*(1-Rs-Cs-Ms) - b.*Cs.*Ms - m.*Cs,Cs);
algal_nullcline=solve(0 ==  s.*Ms.*(1-Rs-Cs-Ms) + s.*Rs.*Ms + b.*Cs*Ms - h.*Ms - (z.*Ms*o.*Cs)/(1+o.*Cs),Ms);
plh=fplot(coral_nullcline(1),[0 1]);set(plh,'color','red','linewidth',2);
plh=fplot(coral_nullcline(2),[0 1]);set(plh,'color','red','linewidth',2);
plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','blue','linewidth',2);
plh=fplot(algal_nullcline(2),[0 1]);set(plh,'color','blue','linewidth',2);
title('Nullclines - Coral-macroalgal dynamics', 'FontSize',11)
xlabel('Macroalgae');
ylabel('Adult Coral');
axis([0 axlim 0 axlim])
%draw axes
xL = xlim;
yL = ylim;
line([0 0], yL, 'Color','k');  %x-axis
line(xL, [0 0], 'Color','k');  %y-axis

