%% Data processing for TAC
clear all
close all
global p dt
%% Simulation data
p = 5;
dt = 0.5;
load('DataTAC/plots_data')
FiguresPositionsSC(q_his,qd,Ncon)
FiguresUeffMSE(this,ueff_his,mse_hisQ,Ncon,N)
fout = FiguresUeffMSETAC();
load('DataTAC/simdata')
ss_mse = mean(DATA2,1);
ss_ueff = mean(DATA4,1);
ss_mse = ss_mse/ss_mse(1);
ss_ueff = ss_ueff/ss_ueff(1);

%% Functions
function FiguresPositionsSC(q_his,qd,Ncon)
global p
f_size = 18;
f_size_ld = 18;
num_size = 14;
lwidth = 1.5;
lwidthd = 1.5;
title_size = 17;
% positions
ConNames = ["CV-STEM controller","Nominal controller","PID controller",...
            "$\mathcal{H}_{\infty}$ controller"];
FileNames = ["sc_cvstem","sc_nominal","sc_pid","sc_hinf"];
Nt = 5500;
tplot = 0:100:Nt;
qdplot = zeros(3,Nt/100);
for i = 1:Nt/100
    qdplot(:,i) = qd(tplot(i));
end
figure
for c = 1:Ncon
    subplot(2,2,c)
    plot(qdplot(1,:),qdplot(2,:),'k--','LineWidth',lwidthd)
    hold on
    for j = 1:p
        plot3(q_his(3*(j-1)+1,1:end,c),q_his(3*(j-1)+2,1:end,c),q_his(3*j,1:end,c),'LineWidth',lwidth)
        hold on
    end
    set(gca,'FontSize',num_size);
    xlabel({'$x$ (km)'},'Interpreter','latex','FontSize',f_size)
    ylabel({'$y$ (km)'},'Interpreter','latex','FontSize',f_size)
    zlabel({'$z$ (km)'},'Interpreter','latex','FontSize',f_size)
    title({ConNames(c)},'Interpreter','latex','FontSize',title_size)
    if c == Ncon
        legend({'$q_d(t)$'},'Interpreter','latex','FontSize',f_size_ld,'Location','southeast')
    end
    grid on
    xlim([-3 3])
    ylim([-3 3])
end
saveas(gcf,fullfile('figures','lvlh'),'epsc')
saveas(gcf,'figures/lvlh.png')
end
% Figures of control effort and MSE
function FiguresUeffMSE(this,ueff_his,mse_his,Ncon,N)
% control effort
f_size = 25;
f_size_ax = 35;
num_size = 25;
lwidth = 5.0;
LineShapes = ["-","-.","--",":"];
figure
for c = 1:Ncon
    plot(this(1:N),(ueff_his(:,c)),LineShapes(c),'LineWidth',lwidth)
    hold on
end
set(gca,'FontSize',num_size);
xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size_ax)
ylabel({'${\int_{0}^{t}\|u(\tau)\|^2d\tau}$'},'Interpreter','latex','FontSize',f_size_ax)
title({'Control Effort'},'Interpreter','latex','FontSize',f_size_ax)
fname1 = "sc_ueff";
saveas(gcf,fullfile('figures',fname1),'epsc')
saveas(gcf,"figures/"+fname1+".png")
% mean squared error
samples = 150;
coeff = ones(1,samples)/samples;
figure
for c = 1:Ncon
    avg = filter(coeff,1,(mse_his(:,c)));
    plot(this(1:N),log10(avg),LineShapes(c),'LineWidth',lwidth)
    hold on
end
ylim([-3 2.5])
set(gca,'FontSize',num_size);
xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size_ax)
ylabel({'$\log_{10}(\|x-x_d\|^2)$'},'Interpreter','latex','FontSize',f_size_ax)
legend({'CV-STEM','Nominal Controller','PID','$\mathcal{H}_{\infty}$'},'Interpreter','latex','FontSize',f_size)
title({'Tracking Error'},'Interpreter','latex','FontSize',f_size_ax)
fname2 = 'sc_mse';
saveas(gcf,fullfile('figures',fname2),'epsc')
saveas(gcf,"figures/"+fname2+".png")
end

function [fout] = FiguresUeffMSETAC()
global dt
% control effort
f_size = 22;
lwidth = 2.0;
num_size = 14;
mksize = 6;
f_size_led = 15;
LineShapes = ["-","-.","--",":"];
matlabG = [0.4660 0.6740 0.1880];
xo = linspace(130,520,1000);
yo = ones(1,1000);
load('DataTAC/simdata')
ss_mse = mean(DATA2,1);
ss_ueff = mean(DATA4,1);
load('DataTAC/simdata_CVSTEMs')
ss_mse2 = mean(DATA2,1);
ss_ueff2 = mean(DATA4,1);
figure
subplot(2,1,1)
for c = 1:length(ss_mse)
    semilogy(xo,(ss_mse(c)/ss_mse(1))*yo,LineShapes(c),'LineWidth',lwidth)
    hold on
end
Ufreqs = [300,400,500,600,700,800,900,1000];
pg = semilogy(Ufreqs*dt,(ss_mse2(2:end)/ss_mse(1)),'-o','MarkerSize',mksize,...
         'MarkerEdgeColor',matlabG,'MarkerFaceColor',matlabG);
hold on
fout = ss_mse2(2:end)/ss_mse(1);
set(gca,'FontSize',num_size);
grid on
ylabel({'$\displaystyle \lim_{t\to 500}\|x-x_d\|^2$'},'Interpreter','latex','FontSize',f_size)
xlim([100,550])
ylim([0.4,10^2])
legend(pg,{'CV-STEM with different sampling periods'},'Interpreter',...
    'latex','FontSize',f_size_led,'Location','best')
subplot(2,1,2)
for c = 1:length(ss_mse)
    plot(xo,ss_ueff(c)/ss_ueff(1)*yo,LineShapes(c),'LineWidth',lwidth)
    hold on
end
plot(Ufreqs*dt,ss_ueff2(2:end)/ss_ueff(1),'-o','MarkerSize',mksize,...
    'MarkerEdgeColor',matlabG,'MarkerFaceColor',matlabG);
set(gca,'FontSize',num_size);
grid on
xlabel({'${\rm Sampling~period}~\Delta t$ (s)'},'Interpreter','latex','FontSize',f_size)
ylabel({'${\int_{0}^{500}\|u(t)\|^2dt}$'},'Interpreter','latex','FontSize',f_size)
xlim([100,550])
ylim([0.73,1.33])
lgd = legend({'CV-STEM','Nominal Controller','PID','$\mathcal{H}_{\infty}$'},...
    'Interpreter','latex','FontSize',f_size_led,'Location','best');
lgd.NumColumns = 4;
fname = 'sc_mse_ueff_0';
saveas(gcf,fullfile('figures',fname),'epsc')
saveas(gcf,"figures/"+fname+".png")
end