clear; close all; clc;

%% choose parameters

kb=0; % breakup rate
epsilon=1e-3; % dimensionless receptor radius
% m=1e3; % number of receptors
u0=600; % extracellular concentration
D=1e3; % diffusion coefficient
R=1; % radius
Jmax=4*pi*D*R*u0;
kcs=[1e1,1e2,1e3,1e4]; % turnover (catalysis) rates

%% make data for plots
msmall=1e2;
mbig=1e5;
Ns=logspace(log10(msmall),log10(mbig),1e2); % Ns to plot.

as=zeros(length(Ns),length(kcs));
abp=zeros(length(Ns),1);
for i=1:length(Ns)
    m=Ns(i);
    kappa=epsilon*m/pi;
    abp(i)=kappa/(kappa+1);
    for j=1:length(kcs)
        chic=kcs(j)*m/Jmax;
        chib=kb*m/Jmax;
        a=(1/2).*(1+chic+kappa.^(-1).*(chib+chic)+(-1).*((-4).*chic+(1+chic+kappa.^(-1).*(chib+ ...
            chic)).^2).^(1/2));
        as(i,j)=a;
    end
end


%% plotting
figure(3); subplot(2,1,1)
hold all

lw=3; %linewidth

plot(Ns,as(:,1),'r','LineWidth',lw)
plot(Ns,as(:,2),'b','LineWidth',lw)
plot(Ns,as(:,3),'k','LineWidth',lw)
plot(Ns,as(:,4),'m','LineWidth',lw)
plot(Ns,abp,'color',[0 .75 0],'LineWidth',lw)
set(gca, 'XScale', 'log')
ll=legend(...
    '$k_{\mbox{c}}=10\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{2}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{3}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=10^{4}\,$s$^{-1}$',...
    '$k_{\mbox{c}}=\infty$ ($J_{\mbox{bp}}$)',...
    'Location','northwest');
xl=xlabel('$N$ (number of receptors)');
yl=ylabel('$J_{*}/J_{\mbox{max}}$');
set(ll,'Interpreter','Latex')
set(xl,'Interpreter','Latex')
set(yl,'Interpreter','Latex')
set(gca,'fontsize',16)

%% choose parameters

kb=0; % breakup rate
epsilon=1e-3; % dimensionless receptor radius
% m=1e3; % number of receptors
u0=600; % extracellular concentration
D=1e3; % diffusion coefficient
R=1; % radius
Jmax=4*pi*D*R*u0;



%% make data for plots

kcsmall=1e1;
kcbig=1e4;
kcs=logspace(log10(kcsmall),log10(kcbig),1e2); % kc's to plot.

N12s=zeros(length(kcs),1);
abp=zeros(length(kcs),1);

for i=1:length(kcs)
    kc=kcs(i);
    N12s(i)=pi*(kb+kc+2*D*R*u0*epsilon)/(kc*epsilon);
end



%% plotting
% figure('units','inch','position',[20,10,8,6]);
subplot(2,1,2)
hold all

lw=3; %linewidth

plot(kcs,N12s,'k','LineWidth',lw)
plot(kcs,pi/epsilon*ones(length(kcs),1),'b--','LineWidth',lw)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ll=legend(...
    '$k_{\mbox{c}}<\infty$',...
    '$k_{\mbox{c}}=\infty$ ($J_{\mbox{bp}}$)',...
    'Location','northeast');
xl=xlabel('$k_{\mbox{c}}$ (receptor turnover rate) [s$^{-1}$]');
yl=ylabel('$N_{\mbox{half}}$');
set(ll,'Interpreter','Latex')
set(xl,'Interpreter','Latex')
set(yl,'Interpreter','Latex')
set(gca,'fontsize',16)

