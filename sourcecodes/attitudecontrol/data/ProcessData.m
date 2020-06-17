%% Data processing for TAC
clear all
close all
syms t
global dt qd dqd_dt dqd_dt_2
%% Desired trajectories
% Desired Trajectory q1d and its derivatives
q1d_sym = 0.3*sin(2*pi*0.1*t);
dq1d_dt_sym = diff(q1d_sym);
dq1d_dt_2_sym = diff(dq1d_dt_sym);
q1d = matlabFunction(q1d_sym);
dq1d_dt = matlabFunction(dq1d_dt_sym);
dq1d_dt_2 = matlabFunction(dq1d_dt_2_sym);
% Desired Trajectory q2d and its derivatives
q2d_sym = 0.2*sin(2*pi*0.2*t+pi/6);
dq2d_dt_sym = diff(q2d_sym);
dq2d_dt_2_sym = diff(dq2d_dt_sym);
q2d = matlabFunction(q2d_sym);
dq2d_dt = matlabFunction(dq2d_dt_sym);
dq2d_dt_2 = matlabFunction(dq2d_dt_2_sym);
% Desired Trajectory q3d and its derivatives
q3d_sym = 0*t;
dq3d_dt_sym = diff(q3d_sym);
dq3d_dt_2_sym = diff(dq3d_dt_sym);
q3d = @(t) 0*t;
dq3d_dt = @(t) 0*t;
dq3d_dt_2 = @(t) 0*t;
% Desired Trajectory
qd = @(t) [q1d(t);q2d(t);q3d(t)];
dqd_dt = @(t) [dq1d_dt(t);dq2d_dt(t);dq3d_dt(t)];
dqd_dt_2 = @(t) [dq1d_dt_2(t);dq2d_dt_2(t);dq3d_dt_2(t)];
%% Simulation data
dt = 0.1;
load('DataTAC/plots_data0')
FiguresMRP(t_his,X_his,Ncon,N)
FiguresUeffMSE(t_his,ueff_his,mse_hisQ,Ncon,N)
FiguresUeffMSETAC()
load('DataTAC/simdata')
ss_mse = mean(DATA2,1);
ss_ueff = mean(DATA4,1);
ss_mse = ss_mse/ss_mse(1);
ss_ueff = ss_ueff/ss_ueff(1);
%% Functions
% Figures of modified Rodrigues parameter
function FiguresMRP(t_his,X_his,Ncon,N)
global dt qd dqd_dt
tf = N*dt;
qdhis = qd(t_his);
dqdthis = dqd_dt(t_his);
f_size_label = 20;
f_size_legend = 15;
lwidth = 2.0;
num_size = 12;
figure
LineShapes = ["-","-.","--",":"];
samples = 150;
coeff = ones(1,samples)/samples;
TPval = 0.7;
for is = 1:3
    subplot(3,2,is*2-1)
    for c = 1:Ncon
        avg = filter(coeff,1,(X_his(is,:,c)-qdhis(is,:)));
        pc = plot(t_his,avg,LineShapes(c),'LineWidth',lwidth);
        pc.Color(4) = TPval;
        hold on
    end
    set(gca,'FontSize',num_size);
    ylab = "$q_"+num2str(is)+"-q_{"+num2str(is)+"d}$";
    ylabel({ylab},'Interpreter','latex','FontSize',f_size_label)
    if is == 3
        xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size_label)
    end
    if is == 1
        lgd = legend({'CV-STEM','Nominal Controller','PID','$\mathcal{H}_{\infty}$'}, ...
        'Interpreter','latex','FontSize',f_size_legend,'Location','best');
        lgd.NumColumns = 4;
    end
    xlim([0 tf])
    grid on
    subplot(3,2,is*2)
    for c = 1:Ncon
        avg = filter(coeff,1,(X_his(is+3,:,c)-dqdthis(is,:)));
        pc = plot(t_his,avg,LineShapes(c),'LineWidth',lwidth);
        pc.Color(4) = TPval;
        hold on
    end
    set(gca,'FontSize',num_size);
    ylab = "$\dot{q}_"+num2str(is)+"-\dot{q}_{"+num2str(is)+"d}$";
    ylabel({ylab},'Interpreter','latex','FontSize',f_size_label)
    xlim([0 tf])
    if is == 3
        xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size_label)
    end
    grid on
end
saveas(gcf,fullfile('figures',"test"),'epsc')
saveas(gcf,"figures/test.png")
end
% Figures of control effort and MSE
function FiguresUeffMSE(t_his,ueff_his,mse_his,Ncon,N)
% control effort
f_size = 35;
f_size_led = 25;
num_size = 25;
lwidth = 5.0;
LineShapes = ["-","-.","--",":"];
TPval = 0.7;
figure
for i = 1:Ncon
    pc = plot(t_his(1:N),(ueff_his(:,i)),LineShapes(i),'LineWidth',lwidth);
    pc.Color(4) = TPval;
    hold on
end
set(gca,'FontSize',num_size);
xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size)
ylabel({'${\int_{0}^{t}\|u(\tau)\|^2d\tau}$'},'Interpreter','latex','FontSize',f_size)
title({'Control Effort'},'Interpreter','latex','FontSize',f_size)
fname1 = 'ueff';
saveas(gcf,fullfile('figures',fname1),'epsc')
saveas(gcf,"figures/"+fname1+".png")
% mean squared error
samples = 150;
coeff = ones(1,samples)/samples;
figure
for i = 1:Ncon
    avg = filter(coeff,1,(mse_his(:,i)));
    pc = plot(t_his(1:N),(avg),LineShapes(i),'LineWidth',lwidth);
    pc.Color(4) = TPval;
    hold on
end
set(gca,'FontSize',num_size);
xlabel({'${\rm time}$ (s)'},'Interpreter','latex','FontSize',f_size)
ylabel({'$\|x(t)-x_d(t)\|^2$'},'Interpreter','latex','FontSize',f_size)
title({'Tracking Error'},'Interpreter','latex','FontSize',f_size)
ylim([0 0.5])
legend({'CV-STEM','Nominal Controller','PID','$\mathcal{H}_{\infty}$'},'Interpreter','latex','FontSize',f_size_led)
fname2 = 'mse';
saveas(gcf,fullfile('figures',fname2),'epsc')
saveas(gcf,"figures/"+fname2+".png")
end
% Figures of control effort and MSE (for TAC)
function FiguresUeffMSETAC()
global dt
% control effort
f_size = 22;
lwidth = 2.0;
num_size = 13;
mksize = 6;
f_size_led = 15;
LineShapes = ["-","-.","--",":"];
matlabG = [0.4660 0.6740 0.1880];
xo = linspace(8,52,1000);
yo = ones(1,1000);
load('DataTAC/simdata')
ss_mse = mean(DATA2,1);
ss_ueff = mean(DATA4,1);
load('DataTAC/simdata_CVSTEMs')
ss_mse2 = mean(DATA2,1);
ss_ueff2 = mean(DATA4,1);
ss_ueff
ss_ueff2
figure
subplot(2,1,1)
for c = 1:length(ss_mse)
    plot(xo,ss_mse(c)/ss_mse(1)*yo,LineShapes(c),'LineWidth',lwidth)
    hold on
end
Ufreqs = [100,150,200,250,300,350,400,450,500];
pg = plot(Ufreqs*dt,ss_mse2/ss_mse(1),'-o','MarkerSize',mksize,...
         'MarkerEdgeColor',matlabG,'MarkerFaceColor',matlabG);
hold on
set(gca,'FontSize',num_size);
grid on
ylabel({'$\displaystyle\lim_{t\to 50}\|x-x_d\|^2$'},'Interpreter','latex','FontSize',f_size)
xlim([5,55])
ylim([0.3,4.8])
legend(pg,{'CV-STEM with different sampling periods'},'Interpreter',...
    'latex','FontSize',f_size_led,'Location','best')
subplot(2,1,2)
for c = 1:length(ss_mse)
    plot(xo,ss_ueff(c)/ss_ueff(1)*yo,LineShapes(c),'LineWidth',lwidth)
    hold on
end
plot(Ufreqs*dt,ss_ueff2/ss_ueff(1),'-o','MarkerSize',mksize,...
    'MarkerEdgeColor',matlabG,'MarkerFaceColor',matlabG);
set(gca,'FontSize',num_size);
grid on
xlabel({'${\rm Sampling~period}~\Delta t$ (s)'},'Interpreter','latex','FontSize',f_size)
ylabel({'${\int_{0}^{50}\|u(t)\|^2dt}$'},'Interpreter','latex','FontSize',f_size)
xlim([5,55])
ylim([0.92,1.6])
lgd = legend({'CV-STEM','Nominal Controller','PID','$\mathcal{H}_{\infty}$'},...
    'Interpreter','latex','FontSize',f_size_led,'Location','best');
lgd.NumColumns = 4;
fname = "mse_ueff_0";
saveas(gcf,fullfile('figures',fname),'epsc')
saveas(gcf,"figures/"+fname+".png")
end