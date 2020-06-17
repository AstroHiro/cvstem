%% Spacecraft Attitude Control using CV-STEMs with various update frequencies
tic
clear all
close all
syms t
global dt qd dqd_dt dqd_dt_2 od d
global Jtot Kr Kp Ki Kd Lambda R alp umax_tot gamhinf
rng(226) % Let Nsim = 1 or 60
%% Constants
Nsim = 60;
n = 6; % #states
m = 3; % #inputs
d = 1; % dimension of Wiener process
dt = 0.1;
R = 1*eye(3);
alp = 1e-3;
umax_tot = 700;
J11 = 150; J12 = 0; J13 = -100; 
J22 = 270; J23 = 0; J33 = 300;
J21 = J12; J32 = J23; J31 = J13;
Jtot = [J11 J12 J13;J21 J22 J23;J31 J32 J33];
Lambda = 1*eye(3);
Kr = 100*eye(3);
Kp = 1300*eye(3);
Ki = 300*eye(3);
Kd = 1300*eye(3);
N = 500;
% Initial Values
t0 = 0;
q0 = [0.9;-0.9;0.7];
dqdt0 = [0.6;0.7;-0.5];
X0 = [q0;dqdt0];
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
% Desired Input
qd_sym = [q1d_sym;q2d_sym;q3d_sym];
dqd_dt_sym = [dq1d_dt_sym;dq2d_dt_sym;dq3d_dt_sym];
Zqd_sym = 1/2*(eye(3)*((1-qd_sym.'*qd_sym)/2)+qd_sym*qd_sym.'+SkewMat(qd_sym));
od_sym = Zqd_sym\dqd_dt_sym;
dod_dt_sym = diff(od_sym);
ud_sym = Jtot*dod_dt_sym-SkewMat(Jtot*od_sym)*od_sym;
od = matlabFunction(od_sym);
ud = matlabFunction(ud_sym);

%% Simulation using CV-STEMs with various update frequencies
Ufreqs = [100,150,200,250,300,350,400,450,500];
Ncon = length(Ufreqs);
DATA1 = zeros(Nsim,Ncon);
DATA2 = zeros(Nsim,Ncon);
DATA3 = zeros(Nsim,Ncon);
DATA4 = zeros(Nsim,Ncon);
DATA5 = zeros(N,Ncon,Nsim);
for nsim = 1:Nsim
    % Initialization
    mse0 = 0;
    ueff0 = 0;
    t = t0;
    Xmat = repmat(X0,1,Ncon);
    ueff_mat = ueff0*ones(1,Ncon);
    t_his = zeros(1,N+1);
    X_his = zeros(n,N+1,Ncon);
    u_his = zeros(m,N,Ncon);
    mse_hisQ = zeros(N,Ncon);
    ueff_his = zeros(N,Ncon);
    t_his(:,1) = t;
    for c = 1:Ncon
        X_his(:,1,c) = X0; 
    end
    [K,Pprev,Qprev,cvx_status] = CVSTEM0(t0,X0);
    gam = 1e6;
    Pprevs = zeros(3,3,Ncon);
    Qprevs = zeros(3,3,Ncon);
    Kss = zeros(3,3,Ncon);
    gams = zeros(1,Ncon);
    for c = 1:Ncon
        Pprevs(:,:,c) = Pprev;
        Qprevs(:,:,c) = Qprev;
        Kss(:,:,c) = eye(3);
        gams(:,c) = gam;
    end
    % Spacecraft Attitude Control with various CV-STEMs
    for k = 1:N
        tex = "SIM2 Current step: k = "+num2str(k)+", Current sim: nsim = "+num2str(nsim)+"\n";
        if rem(k-1,10) == 0
            fprintf(tex)
        end
        dW = mvnrnd(zeros(d,1),dt*eye(d))'; % Same noise for all dynamics
        for c = 1:Ncon
            X = Xmat(:,c);
            G = sde_diffusion(t,X);
            ueff = ueff_mat(c);
            [M,C] = get_MC(X);
            Pprev = Pprevs(:,:,c);
            Qprev = Qprevs(:,:,c);
            Ks = Kss(:,:,c);
            gam = gams(c);
            Ufreq = Ufreqs(c);
            [u,Pprev,Qprev,gam,Ks] = AttitudeCVSTEM(t,X,Pprev,Qprev,gam,G,Ks,k,Ufreq);
            Pprevs(:,:,c) = Pprev;
            Qprevs(:,:,c) = Qprev;
            Kss(:,:,c) = Ks;
            gams(:,c) = gam;
            drift  = @(X,u) AttitudeSC(X,u);
            Xd = [qd(t);dqd_dt(t)];
            mse_hisQ(k,c) = norm(X-Xd)^2;
            X = rk4(dt,X,u,drift)+[zeros(3,1);G]*dW;
            ueff = ueff+norm(u)^2*dt;
            X_his(:,k+1,c) = X;
            u_his(:,k,c) = u;
            ueff_his(k,c) = ueff;
            Xmat(:,c) = X;
            ueff_mat(c) = ueff;
        end
        t = t+dt;
        t_his(:,k+1) = t;
    end
    % Summarize simulation data
    if Nsim == 1
        save('data/plots_data.mat','t_his','ueff_his','mse_hisQ','Ncon','N')
    end
    samples = 150;
    peak_ns = zeros(1,Ncon);
    end_ns = zeros(1,Ncon);
    peak_us = zeros(1,Ncon);
    end_ueffs = zeros(1,Ncon);
    for c = 1:Ncon
        peak_ns(c) = max(mse_hisQ(:,c));
        end_ns(c) = sum(mse_hisQ((end-(samples-1)):end,c))/samples;
        peak_us(c) = max(sum(u_his(:,:,c).*u_his(:,:,c)));
        end_ueffs(c) = ueff_his(end,c);
    end
    DATA1(nsim,:) = peak_ns;
    DATA2(nsim,:) = end_ns;
    DATA3(nsim,:) = peak_us;
    DATA4(nsim,:) = end_ueffs;
    DATA5(:,:,nsim) = mse_hisQ;
end
if Nsim == 60
    save('data/simdata_CVSTEMs.mat','DATA1','DATA2','DATA3','DATA4','DATA5')
end
sim_time = toc

%% Functions
% CV-STEM for initialization
function [K,P,Q,cvx_status] = CVSTEM0(t,X)
global Kr R alp umax_tot
s = get_s(t,X);
[M,C] = get_MC(X);
q = X(1:3);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
A = -M\(C+Kr);
B = inv(M);
Bm = inv(Zq');
N = length(X);
n = N/2;
I = eye(n);
cvx_begin sdp quiet
    variable Q(n,n) semidefinite
    variables tau
    variable gam nonnegative
    variable nu nonnegative
    variable chi nonnegative
    minimize (tau)
    subject to
        I <= Q <= chi*I;
        A*Q+Q*A'-nu*B*inv(R)*B'+gam*I <= 0;
        [Q*Bm*inv(R)*B'+B*inv(R)*Bm'*Q+gam*I+nu*B*inv(R)*B'-2*alp*Q Q;Q 1/2/alp*inv(M)*nu] >= 0;
        [tau-n*chi chi;chi nu/trace(M)] >= 0;
        nu*norm(inv(R)*B')*norm(s) <= lambda_min(Q)*umax_tot;
cvx_end
Pinv = Q/nu;
P = inv(Pinv);
K = inv(R)*B'/Pinv;
end
% Compute composite state s
function s = get_s(t,X)
global qd dqd_dt Lambda
q = X(1:3);
dqdt = X(4:6);
s = (dqdt-dqd_dt(t)+Lambda*(q-qd(t)));
end
% Nominal controller that guarantees exponential convergence
function u = ExpController(t,X)
global Jtot qd dqd_dt dqd_dt_2 Kr Lambda
q = X(1:3);
dqdt = X(4:6);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
if isnan(Zq)
    error('Not a number')
end
omega = Zq\dqdt;
dZqdt =1/2*(eye(3)*(-q'*dqdt)+dqdt*q'+q*dqdt'+SkewMat(dqdt));
omega_r = Zq\(dqd_dt(t)+Lambda*(qd(t)-q));
omega_r_dot = Zq\(-dZqdt*omega_r+(dqd_dt_2(t)+Lambda*(dqd_dt(t)-dqdt)));
u = Jtot*omega_r_dot-SkewMat(Jtot*omega)*omega_r-Kr*(omega-omega_r);
end
% PID controller
function u = pid_controller(t,X)
global dt qd dqd_dt Kp Ki Kd
persistent int_e
if isempty(int_e)
    int_e = zeros(3,1);
end
q = X(1:3);
dqdt = X(4:6);
int_e = int_e+(q-qd(t))*dt;
if norm(int_e) >= 10.0
    int_e = zeros(3,1);
end
u = -Kp*(q-qd(t))-Ki*int_e-Kd*(dqdt-dqd_dt(t));
end
% Diffusion term of SDE 
function G = sde_diffusion(t,X)
global d
n = length(X)/2;
G = ones(n,d)*0.2;
end
% Find matrices M and C in Lagrangian dynamics
function [M,C] = get_MC(X)
% Q = [q;qdot]
% q = Modified Rodrigues Parameters
global Jtot
q = X(1:3);
dqdt = X(4:6);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
dZqdt =1/2*(eye(3)*(-q'*dqdt)+dqdt*q'+q*dqdt'+SkewMat(dqdt));
M = Zq'\Jtot/Zq;
C = -Zq'\Jtot/Zq*dZqdt/Zq-Zq'\SkewMat(Jtot*(Zq\dqdt))/Zq;
end
% Relaxed CVSTEM algorithm
function [K,Pnext,Qnext,gam_next,cvx_status] = RelaxedCVSTEM(t,X,P,Q,gam,G)
global dt Kr R alp umax_tot
persistent ss_error
if isempty(ss_error)
        ss_error = -Inf;
end
s = get_s(t,X);
[M,C] = get_MC(X);
mink = min(eig(Kr));
lmax = max(eig(M+P));
lmin = min(eig(M+P));
c = trace(inv(M)*G*G'*inv(M));
ss_error_next = c*trace(M+P)/(2*mink/lmax+2*alp)/lmin;
[f1,f2,f3] = check_feasibility(t,X,P,gam);
if norm(P) == 0
    f1 = Inf;
    f2 = Inf;
end
if (ss_error_next-ss_error >= 0) || (f1 >= 0) || (f2 >= 0) || (f3 >= 0)
    ss_error = ss_error_next;
    q = X(1:3);
    Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
    Kz = Zq'\Kr/Zq;
    A = -M\(C+Kz);
    Bm = inv(Zq');
    B = inv(M)*Bm;
    N = length(X);
    n = N/2;
    I = eye(n);
    Qprev = Q;
    cvx_begin sdp quiet
        variable Q(n,n) semidefinite
        variables tau
        variable gam nonnegative
        variable nu nonnegative
        variable chi nonnegative
        minimize (tau)
        subject to
            I <= Q <= chi*I;
            dQdt = (Q-Qprev)/dt;
            -dQdt+A*Q+Q*A'-nu*B*inv(R)*B'+gam*I <= 0;
            [Q*Bm*inv(R)*B'+B*inv(R)*Bm'*Q+gam*I+nu*B*inv(R)*B'-2*alp*Q Q;Q 1/2/alp*inv(M)*nu] >= 0;
            [tau-n*chi chi;chi nu/trace(M)] >= 0;
            nu*norm(inv(R)*B')*norm(s) <= lambda_min(Q)*umax_tot
    cvx_end
    Qnext = Q;
    Pinv = Q/nu;
    Pnext= inv(Pinv);
    K = inv(R)*B'/Pinv;
    gam_next = gam/nu;
else
    B = inv(M);
    K = inv(R)*B'*P;
    Pnext = P;
    Qnext = Q;
    gam_next = gam;
    cvx_status = 'NotPerformed';
end
end
% Check feasibility of CV-STEM constraints
function [f1,f2,f3] = check_feasibility(t,X,P,gam)
global Kr R alp umax_tot
s = get_s(t,X);
[M,C] = get_MC(X);
q = X(1:3);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
A = -M\(C+Kr);
B = inv(M);
Bm = inv(Zq');
f1 = s'*(P*A+A'*P-P*B*inv(R)*B'*P+gam*P^2)*s;
f2 = -s'*(Bm*inv(R)*B'*P+P*B*inv(R)*Bm'+gam*P^2+P*B*inv(R)*B'*P-2*alp*(M+P))*s;
f3 = norm(-inv(R)*B'*P*s)-umax_tot;
end
% CVSTEM for spacecraft attitude control
function [u,Pnext,Qnext,gam,Ks] = AttitudeCVSTEM(t,X,P,Q,gam,G,Ks,k,Ufreq)
s = get_s(t,X);
un = ExpController(t,X);
if rem(k-1,Ufreq) == 0
    [Ks,P,Q,gam,cvx_status] = RelaxedCVSTEM(t,X,P,Q,gam,G);
end
if isnan(Ks)
    us = [0;0;0];
    Pnext = zeros(3,3);
    Qnext = zeros(3,3);
else
    us = -Ks*s;
    Pnext = P;
    Qnext = Q;
end
u = un+us;
end
% H infinity control
function [u,P,Q,cvx_status] = Hinfinity(X,Xd,Qprev)
global dt R gamhinf
gam = gamhinf;
A = AmatSDC(X);
B = BmatSDC(X);
n = length(X);
I = eye(n);
cvx_begin sdp quiet
    variable Q(n,n) semidefinite
    variable g_inv nonnegative
    minimize (g_inv)
    subject to
        dQdt = (Q-Qprev)/dt;
        [-dQdt+A*Q+Q*A'-B*inv(R)*B' Q;Q -I*g_inv] <= 0;
        g_inv >= 1/gam;
cvx_end
if g_inv >= 1*1e-6
    flag = 1;
    cvx_begin sdp quiet
        variable Q(n,n) semidefinite
        minimize (g_inv)
        subject to
            [A*Q+Q*A'-B*inv(R)*B' Q;Q -I*g_inv] <= 0;
            g_inv >= 1/gam;
    cvx_end
else
    flag = 0;
end
Pinv = Q;
P = inv(Pinv);
K = inv(R)*B'/Pinv;
u = -K*(X-Xd);
end
% A matrix of SDC formulation
function A = AmatSDC(X)
global Jtot
q = X(1:3);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
if isnan(Zq)
    error('Not a number')
end
dqdt = X(4:6);
om = Zq\dqdt;
dZqdt =1/2*(eye(3)*(-q'*dqdt)+dqdt*q'+q*dqdt'+SkewMat(dqdt));
A11 = zeros(3,3);
A12 = eye(3);
A21 = zeros(3,3);
A22 = (dZqdt+Zq*(Jtot\SkewMat(Jtot*om)))/Zq;
A = [A11 A12;A21 A22];
end
% B matrix of SDC formulation 
function B = BmatSDC(X)
global Jtot
q = X(1:3);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
if isnan(Zq)
    error('Not a number')
end
B = [zeros(3);Zq/Jtot];
end
% Spacecraft attitude dynamics 
function dXdt = AttitudeSC(X,u)
[M,C] = get_MC(X);
q = X(1:3);
dqdt = X(4:6);
Zq = 1/2*(eye(3)*((1-q'*q)/2)+q*q'+SkewMat(q));
tau = Zq'\u; % Disturbance torque
dqdt_2 = M\(-C*dqdt+tau);
dXdt = [dqdt;dqdt_2];
end
% Cross-product in matrix form
function skew = SkewMat(X)
x1 = X(1); x2 = X(2); x3 = X(3);
skew = [  0 -x3  x2;
         x3   0 -x1;
        -x2  x1   0];
end
% 4th order Runge-Kutta
function Xnext = rk4(dt,X,u,f)
N = 1;
dt_runge = dt/N;
t = 0;
for i = 1:N
    k1 = f(X,u);
    k2 = f(X+dt_runge*k1/2,u);
    k3 = f(X+dt_runge*k2/2,u);
    k4 = f(X+dt_runge*k3,u);
    X = X+dt_runge/6*(k1+2*k2+2*k3+k4);
    t = t+dt_runge;
end
Xnext = X;
end