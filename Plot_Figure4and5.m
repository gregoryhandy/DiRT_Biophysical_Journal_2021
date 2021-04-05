%%
% Reproduces Figures 4 and 5 from
% Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
% Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
%
% Written by Gregory Handy and Sean Lawley, 04/05/2021
%%
close all; close all; clc;

%% choose parameters

kb=0; % breakup rate
epsilon=1e-3; % dimensionless receptor radius
m=1e3; % number of receptors
% u0=600; % extracellular concentration
D=1e3; % diffusion coefficient
R=1; % radius

kcs=[1e1,1e2,1e3,1e4]; % turnover (catalysis) rates


%% make data for plots

u0small=6*1e-1;
u0big=6*1e5;

u0s=logspace(log10(u0small),log10(u0big),1e2); % ms to plot.

Js=zeros(length(u0s),length(kcs));
Jmms=zeros(length(u0s),length(kcs));
Jbp=zeros(length(u0s),1);

for i=1:length(u0s)
    u0=u0s(i);
    Jmax=4*pi*D*R*u0;
    kappa=epsilon*m/pi;
    Jbp(i)=kappa/(kappa+1)*Jmax;
    for j=1:length(kcs)
        chic=kcs(j)*m/Jmax;
        chib=kb*m/Jmax;
        a=(1/2).*(1+chic+kappa.^(-1).*(chib+chic)+(-1).*((-4).*chic+(1+chic+kappa.^(-1).*(chib+ ...
            chic)).^2).^(1/2));
        Js(i,j)=a*Jmax;
        kc=kcs(j);
        Vmax=m*kc;
        Km=Vmax*u0/Jbp(i);
        Jmms(i,j)=Vmax*u0/(Km+u0);
    end
end


%% plotting
% figure('units','inch','position',[20,10,16,12]);
figure(5); clf;
hold on

lw=3; %linewidth

u0sMOLAR=u0s/600;

p1=plot(u0sMOLAR,Js(:,1),'r','LineWidth',lw);
plot(u0sMOLAR,kcs(1)*m*ones(length(u0s),1),'r--','LineWidth',lw/2)

p2=plot(u0sMOLAR,Js(:,2),'b','LineWidth',lw);
plot(u0sMOLAR,kcs(2)*m*ones(length(u0s),1),'b--','LineWidth',lw/2)

p3=plot(u0sMOLAR,Js(:,3),'k','LineWidth',lw);
plot(u0sMOLAR,kcs(3)*m*ones(length(u0s),1),'k--','LineWidth',lw/2)

p4=plot(u0sMOLAR,Js(:,4),'m','LineWidth',lw);
plot(u0sMOLAR,kcs(4)*m*ones(length(u0s),1),'m--','LineWidth',lw/2)

p5=plot(u0sMOLAR,Jbp,'color',[0 .75 0],'LineWidth',lw);

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ll=legend([p1,p2,p3,p4,p5],...
    '$k_{\mbox{c}}=10\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{2}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{3}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{4}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=\infty$ ($J_{\mbox{bp}}$)',...
    'Location','northwest');
xl=xlabel('$u_{0}$ (extracellular conc.) [$\mu$M]');
yl=ylabel('$J_{*}$ [s$^{-1}$]');
set(ll,'Interpreter','Latex')
set(xl,'Interpreter','Latex')
set(yl,'Interpreter','Latex')
set(gca,'fontsize',16)
ylim([1e3,1e9])
xlim([1e-3,1e3])

%%

figure(4); clf;
hold on
p1=plot(u0sMOLAR,Js(:,1)./(kcs(1)*m),'r','LineWidth',lw);
p11=plot(u0sMOLAR,Jmms(:,1)./(kcs(1)*m),'r:','LineWidth',lw);
p2=plot(u0sMOLAR,Js(:,2)./(kcs(2)*m),'b','LineWidth',lw);
plot(u0sMOLAR,Jmms(:,2)./(kcs(2)*m),'b:','LineWidth',lw);
p3=plot(u0sMOLAR,Js(:,3)./(kcs(3)*m),'k','LineWidth',lw);
plot(u0sMOLAR,Jmms(:,3)./(kcs(3)*m),'k:','LineWidth',lw);
p4=plot(u0sMOLAR,Js(:,4)./(kcs(4)*m),'m','LineWidth',lw);
plot(u0sMOLAR,Jmms(:,4)./(kcs(4)*m),'m:','LineWidth',lw);
set(gca, 'XScale', 'log')
ll=legend([p1,p11,p2,p3,p4],...
    '$k_{\mbox{c}}=10\,$s$^{-1}$, $J_{*}$',...
    '$k_{\mbox{c}}=10\,$s$^{-1}$, $J_{\mbox{mm}}$',...
    '$k_{\mbox{c}}=10^{2}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{3}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{4}\,$s$^{-1}$',...
    'Location','southeast');
xl=xlabel('$u_{0}$ (extracellular conc.) [$\mu$M]');
yl=ylabel('$J/V_{\mbox{max}}$');
set(ll,'Interpreter','Latex')
set(xl,'Interpreter','Latex')
set(yl,'Interpreter','Latex')
set(gca,'fontsize',16)
ylim([0,1.01])
xlim([1e-3,1e3])


