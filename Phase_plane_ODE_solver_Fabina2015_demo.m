% Solve coral-macroalgal dynamics usign equations from Fabina et al. 2015paper

clear;
%constrain parameter space
x=0:0.05:1;
y=0:0.05:1;
spc=[0 0];
for i=1:size(x,2)
    for j=1:size(y,2)
        if x(i)+y(j)<=0.8%<= 1 for sure
            spc=vertcat(spc,[x(i) y(j)]);
        end
    end
end
spc(:,3)=repmat(0.01,size(spc,1),1);
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
p = 1;



figure;
hold on;

%solve and draw solution for parameter space
for ii = 1:length(spc)
    options = odeset('NonNegative',1);
    [T,Y] = ode45(@popdyn_fabina_2015,[0 200],[spc(ii,3),spc(ii,2),spc(ii,1)],options,[a,s,n,g,b,m,h,z,o,r] );
    Rs = Y(:,1); % coral recruits
    As = Y(:,2); % coral adults
    Ms = Y(:,3); % macroalgae
    phaseplane = plot(Ms,As,'Color','k'); % plot solutions
end

title('Phase plane - Coral-macroalgal dynamics', 'FontSize',11)
xlabel('Macroalgae');
ylabel('Adult Coral');
axis([0 1 0 1])
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 Fabina_2015_phaseplane_ermid.eps

%solve for nullclines for constant recruitment
figure; hold;
syms As Ms;
Rs=0.01;
S1=solve(0 == a.*Rs + g.*As.*(1-Rs-As-Ms) - b.*As.*Ms - m.*As,As);
S2=solve(0 ==  s.*Ms.*(1-Rs-As-Ms) + s.*Rs.*Ms + b.*As*Ms - h.*Ms - (z.*Ms*o.*As)/(1+o.*As),Ms);
plh=fplot(S1(1),[0 1]);set(plh,'color','red','linewidth',2);
plh=fplot(S1(2),[0 1]);set(plh,'color','red','linewidth',2);
plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','blue','linewidth',2);
plh=fplot(S2(2),[0 1]);set(plh,'color','blue','linewidth',2);
title('Nullclines - Coral-macroalgal dynamics', 'FontSize',11)
xlabel('Macroalgae');
ylabel('Adult Coral');
axis([-0.1 1 -0.1 1])
