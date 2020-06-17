%% CV-STEM for Spacecraft Phase Synchronization with Fixed Desired Trajectory
tStart = tic;
clear all
close all
global mu kJ2 mj psie0 psiz0 p Tall Rfd d
global dt Kp Ki Kd Lamdd L R alp umax_tot
rng(19940808) % Let Nsim = 1 or 60 for the TAC plots
%% Constants
mu = 3.9862e5; % [km^3/s^2]
re = 6378.15; % [km]
kJ2 = 2.633e10; % [km^5/s^2]
mj = 100;
% Cheaf initial conditions
height = 400;
r0 = re+height;
vx0 = 0;
h0 = sqrt(mu*r0); % r0*sqrt(mu/r0)
Om0 = pi/4;
i0 = pi/4;
th0 = pi/4;
X0 = [r0;vx0;h0;Om0;i0;th0];
% Leader and followers desired trajectory
xe = 2;
ye = 2;
ze = 0.;
psie0 = deg2rad(0.573);
psiz0 = deg2rad(11.46);
p = 5; % # Spacecraft
n = sqrt(mu/r0^3); % mean motion
Xd0 = zeros(3,p);
qd = @(t) [xe*sin(n*t+psie0);ye*cos(n*t+psie0);ze*sin(n*t+psiz0)];
qd_dot = @(t) n*[xe*cos(n*t+psie0);-ye*sin(n*t+psie0);ze*cos(n*t+psiz0)];
qd_ddot = @(t) n^2*[-xe*sin(n*t+psie0);-ye*cos(n*t+psie0);-ze*sin(n*t+psiz0)];
[Rfd,Tall] = GetTransMats(xe,ye,ze);
Qd = @(t) [qd(t);qd_dot(t);qd_ddot(t)];
% Deputies initialization
umax = 0.8;
gam = 0.4;
rc = 0.1;
delta = 0.000001;
k1 = 1;
Lamdd = 1*eye(3);
d0 = 0.2;
v0 = 0.4;
q0_all = zeros(3*p,1);
q0_dot_all = zeros(3*p,1);
for j = 1:p
    xj0 = -d0+2*d0*rand;
    yj0 = -d0+2*d0*rand;
    zj0 = -d0+2*d0*rand;
    q0_all(3*(j-1)+1:3*j) = [xj0;yj0;zj0];
end
%% Simulation with CV-STEM, ExpControl, PID, and Hinf
t0 = 0;
tf = 500;
N = 1000;
dt = (tf-t0)/N;
Ncon = 4;
ueff0 = 0;
K1 = eye(3)*5;
K2 = eye(3)*2;
L = GetLaplacian(K1,-K2);
R = eye(3*p);
alp = 0.001;
d = 1;
umax_tot = 1; % 10 worked
Kp = eye(3*p)*7;
Ki = eye(3*p)*0;
Kd = eye(3*p)*11;
Nsim = 60;
DATA1 = zeros(Nsim,Ncon);
DATA2 = zeros(Nsim,Ncon);
DATA3 = zeros(Nsim,Ncon);
DATA4 = zeros(Nsim,Ncon);
DATA5 = zeros(N,Ncon,Nsim);
for nsim = 1:Nsim
    % Initialization
    t = t0;
    Xmat = repmat(X0,1,Ncon);
    q_all_mat = repmat(q0_all,1,Ncon);
    q_dot_all_mat = repmat(q0_dot_all,1,Ncon);
    ueff_mat = ueff0*ones(1,Ncon);
    t_his = zeros(1,N+1);
    this = zeros(1,N+1);
    Xhis = zeros(6,N+1,Ncon);
    q_his = zeros(3*p,N+1,Ncon);
    q_dot_his = zeros(3*p,N+1,Ncon);
    Uhis = zeros(3*p,N,Ncon);
    mse_hisQ = zeros(N,Ncon);
    ueff_his = zeros(N,Ncon);
    this(1) = t0;
    for c = 1:Ncon
        Xhis(:,1,c) = X0;
        q_his(:,1,c) = q0_all;
        q_dot_his(:,1,c) = q0_dot_all;
    end
    Qprev = eye(3*p);
    Qhinfp = zeros(6*p,6*p);
    Phinf = 1;
    gamhinf = 20;
    Ks = eye(3*p);
    % Spacecraft Phase Synchronization with CV-STEM, ExpControl, PID, and Hinf
    for k = 1:N
        tex = "SIM1 Current step: k = "+num2str(k)+", Current sim: nsim = "+num2str(nsim)+"\n";
        if rem(k-1,10) == 0
            fprintf(tex)
        end
        dW = mvnrnd(zeros(d,1),dt*eye(d))'; % Same noise for all dynamics
        for c = 1:Ncon
            X = Xmat(:,c);
            q_all = q_all_mat(:,c);
            q_dot_all = q_dot_all_mat(:,c);
            ueff = ueff_mat(c);
            [s_d,qr_dot_d,qr_ddot_d] = Get_s_d(q_all,q_dot_all,Qd(t));
            G = sde_diffusion(t,X,q_all,q_dot_all);
            tau_n = ExpControl(X,q_all,q_dot_all,Qd(t));
            if c == 1
                Ufreq = 1;
                tic
                [tau,Qprev,Ks] = PhaseSyncCVSTEM(tau_n,Ks,k,Ufreq,X,q_all,q_dot_all,s_d,Qprev);
                tCVSTEM = toc
            elseif c == 2
                tau = tau_n;
            elseif c == 3
                tau = pid_controller(t,q_all,q_dot_all,Qd(t));
            elseif c == 4
                [tau,Qhinfp] = Hinfinity(X,q_all,q_dot_all,Qd(t),Qhinfp,gamhinf);
            end
            [XX,qq,qqd] = GetXX(X,q_all,q_dot_all,Qd(t));
            ueff = ueff+norm(tau)^2*dt;
            fun = @(XX) RelDynamics(XX,tau);
            XX = rk4(XX,dt,fun)+[zeros(3*p,1);G;zeros(6,1)]*dW;
            Xmat(:,c) = XX(6*p+1:6*p+6);
            q_all_mat(:,c) = XX(1:3*p);
            q_dot_all_mat(:,c) = XX(3*p+1:6*p);
            Xhis(:,k+1,c) = X;
            q_his(:,k+1,c) = q_all;
            q_dot_his(:,k+1,c) = q_dot_all;
            mse_hisQ(k,c) = norm(qq-qqd)^2;
            ueff_his(k,c) = ueff;
            Uhis(:,k,c) = tau;
            ueff_mat(c) = ueff;
        end
        t = t+dt;
        this(k+1) = t;
        if rem(k,100) == 0
            fprintf('k = %d\n',k)
        end
    end
    % Summarize simulation data
    samples = 150;
    peak_ns = zeros(1,Ncon);
    end_ns = zeros(1,Ncon);
    peak_us = zeros(1,Ncon);
    end_ueffs = zeros(1,Ncon);
    for c = 1:Ncon
        peak_ns(c) = max(mse_hisQ(:,c));
        end_ns(c) = sum(mse_hisQ((end-(samples-1)):end,c))/samples;
        peak_us(c) = max(sum(Uhis(:,:,c).*Uhis(:,:,c)));
        end_ueffs(c) = ueff_his(end,c);
    end
    DATA1(nsim,:) = peak_ns;
    DATA2(nsim,:) = end_ns;
    DATA3(nsim,:) = peak_us;
    DATA4(nsim,:) = end_ueffs;
    DATA5(:,:,nsim) = mse_hisQ;
    if Nsim == 1
        save('data/plots_data.mat','this','q_his','Uhis','qd','ueff_his','mse_hisQ','Ncon','N')
    end
end
if Nsim == 60
    save('data/simdata.mat','DATA1','DATA2','DATA3','DATA4','DATA5')
end

%% Simulation using CV-STEMs with various update freqencies
Ufreqs = [200,300,400,500,600,700,800,900,1000];
Ncon = length(Ufreqs);
DATA1 = zeros(Nsim,Ncon);
DATA2 = zeros(Nsim,Ncon);
DATA3 = zeros(Nsim,Ncon);
DATA4 = zeros(Nsim,Ncon);
DATA5 = zeros(N,Ncon,Nsim);
for nsim = 1:Nsim
    % Initialization
    t = t0;
    Xmat = repmat(X0,1,Ncon);
    q_all_mat = repmat(q0_all,1,Ncon);
    q_dot_all_mat = repmat(q0_dot_all,1,Ncon);
    ueff_mat = ueff0*ones(1,Ncon);
    t_his = zeros(1,N+1);
    this = zeros(1,N+1);
    Xhis = zeros(6,N+1,Ncon);
    q_his = zeros(3*p,N+1,Ncon);
    q_dot_his = zeros(3*p,N+1,Ncon);
    Uhis = zeros(3*p,N,Ncon);
    mse_hisQ = zeros(N,Ncon);
    ueff_his = zeros(N,Ncon);
    this(1) = t0;
    Qprev = eye(3*p);
    Qprevs = zeros(3*p,3*p,Ncon);
    Kss = zeros(3*p,3*p,Ncon);
    for c = 1:Ncon
        Xhis(:,1,c) = X0;
        q_his(:,1,c) = q0_all;
        q_dot_his(:,1,c) = q0_dot_all;
        Qprevs(:,:,c) = Qprev;
        Kss(:,:,c) = eye(3*p);
    end
    Qhinfp = zeros(6*p,6*p);
    Phinf = 1;
    gamhinf = 20;
    Ks = eye(3*p);
    % Spacecraft Phase Synchronization with CV-STEM, ExpControl, PID, and Hinf
    for k = 1:N
        tex = "SIM2 Current step: k = "+num2str(k)+", Current sim: nsim = "+num2str(nsim)+"\n";
        if rem(k-1,10) == 0
            fprintf(tex)
        end
        dW = mvnrnd(zeros(d,1),dt*eye(d))'; % Same noise for all dynamics
        for c = 1:Ncon
            X = Xmat(:,c);
            q_all = q_all_mat(:,c);
            q_dot_all = q_dot_all_mat(:,c);
            ueff = ueff_mat(c);
            [s_d,qr_dot_d,qr_ddot_d] = Get_s_d(q_all,q_dot_all,Qd(t));
            G = sde_diffusion(t,X,q_all,q_dot_all);
            tau_n = ExpControl(X,q_all,q_dot_all,Qd(t));
            Qprev = Qprevs(:,:,c);
            Ks = Kss(:,:,c);
            [tau,Qprev,Ks] = PhaseSyncCVSTEM(tau_n,Ks,k,Ufreqs(c),X,q_all,q_dot_all,s_d,Qprev);
            Qprevs(:,:,c) = Qprev;
            Kss(:,:,c) = Ks;
            [XX,qq,qqd] = GetXX(X,q_all,q_dot_all,Qd(t));
            ueff = ueff+norm(tau)^2*dt;
            fun = @(XX) RelDynamics(XX,tau);
            XX = rk4(XX,dt,fun)+[zeros(3*p,1);G;zeros(6,1)]*dW;
            Xmat(:,c) = XX(6*p+1:6*p+6);
            q_all_mat(:,c) = XX(1:3*p);
            q_dot_all_mat(:,c) = XX(3*p+1:6*p);
            Xhis(:,k+1,c) = X;
            q_his(:,k+1,c) = q_all;
            q_dot_his(:,k+1,c) = q_dot_all;
            mse_hisQ(k,c) = norm(qq-qqd)^2;
            ueff_his(k,c) = ueff;
            Uhis(:,k,c) = tau;
            ueff_mat(c) = ueff;
        end
        t = t+dt;
        this(k+1) = t;
        if rem(k,100) == 0
            fprintf('k = %d\n',k)
        end
    end
    % Summarize simulation data
    samples = 150;
    peak_ns = zeros(1,Ncon);
    end_ns = zeros(1,Ncon);
    peak_us = zeros(1,Ncon);
    end_ueffs = zeros(1,Ncon);
    for c = 1:Ncon
        peak_ns(c) = max(mse_hisQ(:,c));
        end_ns(c) = sum(mse_hisQ((end-(samples-1)):end,c))/samples;
        peak_us(c) = max(sum(Uhis(:,:,c).*Uhis(:,:,c)));
        end_ueffs(c) = ueff_his(end,c);
    end
    DATA1(nsim,:) = peak_ns;
    DATA2(nsim,:) = end_ns;
    DATA3(nsim,:) = peak_us;
    DATA4(nsim,:) = end_ueffs;
    DATA5(:,:,nsim) = mse_hisQ;
    if Nsim == 1
        save('data/plots_data2.mat','this','ueff_his','mse_hisQ','Ncon','N')
    end
end
if Nsim == 60
    save('data/simdata_CVSTEMs.mat','DATA1','DATA2','DATA3','DATA4','DATA5')
end
sim_time = toc(tStart)

%% Functions
% Compute transformation matrices
function [Rfd,Tall] = GetTransMats(xe,ye,ze)
global psie0 psiz0 p
% Rfd
[lmin,lmax,PHI] = GetMaxes(xe,ye,ze);
R11 = -xe/lmin*sin(PHI-psie0);
R12 = ye/lmin*cos(PHI-psie0);
R13 = -ze/lmin*sin(PHI-psiz0);
R21 = -xe/lmax*cos(PHI-psie0);
R22 = -ye/lmax*sin(PHI-psie0);
R23 = -ze/lmax*cos(PHI-psiz0);
R31 = -ye*ze/lmin/lmax*cos(psie0-psiz0);
R32 = xe*ze/lmin/lmax*sin(psie0-psiz0);
R33 = xe*ye/lmin/lmax;
R1d = [R11 R12 R13;
       R21 R22 R23;
       R31 R32 R33];
R2d = [1 0 0;0 lmin/lmax 0;0 0 1];
Rfd = R2d*R1d;
% T_j-1
Tall = zeros(3*p,3*p);
idx = 1;
for j = 1:p
    if j > 5/2*idx*(idx+1)
        idx = idx+1;
    end
    sat_num = 5*idx;
    phi = 2*pi/sat_num;
    Tz = [cos((j-1)*phi) -sin((j-1)*phi);sin((j-1)*phi) cos((j-1)*phi)];
    Tz = idx*eye(2)*Tz;
    Tall(3*(j-1)+1:3*j,3*(j-1)+1:3*j) = [Tz zeros(2,1);zeros(1,2) 1];
end
end
% Compute semi-major and -minor axes
function [lmin,lmax,PHI] = GetMaxes(xe,ye,ze)
global psie0 psiz0
numP = (xe^2-ye^2)*sin(2*psie0)+ze^2*sin(2*psiz0);
denP = (xe^2-ye^2)*cos(2*psie0)+ze^2*cos(2*psiz0);
PHI = 1/2*atan2(numP,denP);
xdfun = @(psi) [xe*sin(psi+psie0);ye*cos(psi+psie0);ze*sin(psi+psiz0)];
lfun = @(psi) norm(xdfun(psi));
lmin = lfun(-PHI);
lmax = lfun(3/2*pi-PHI);
end
% Compute Laplacian of Lagrangian system
function L = GetLaplacian(M1,M2)
global p
L = zeros(3*p,3*p);
for j = 1:p
    if j == 1
        L(1:3,1:3) = M1;
        L(1:3,4:6) = M2;
        L(1:3,3*(p-1)+1:3*p) = M2;
    elseif j == p
        L(3*(p-1)+1:3*p,3*(p-1)+1:3*p) = M1;
        L(3*(p-1)+1:3*p,1:3) = M2;
        L(3*(p-1)+1:3*p,3*(p-2)+1:3*(p-1)) = M2;
    else
        L(3*(j-1)+1:3*j,3*(j-1)+1:3*j) = M1;
        L(3*(j-1)+1:3*j,3*j+1:3*(j+1)) = M2;
        L(3*(j-1)+1:3*j,3*(j-2)+1:3*(j-1)) = M2;
    end
end
end
% Compute composite state s_d
function [s_d,qr_dot_d,qr_ddot_d] = Get_s_d(q_all,q_dot_all,Qd_t)
global Lamdd Tall Rfd p
qd = Qd_t(1:3);
qd_dot = Qd_t(4:6);
qd_ddot = Qd_t(7:9);
R_d = kron(eye(p),Rfd);
Lam_all = kron(eye(p),Lamdd);
qd_all = repmat(qd,p,1);
qd_dot_all = repmat(qd_dot,p,1);
qd_ddot_all = repmat(qd_ddot,p,1);
q_d = R_d*q_all;
q_dot_d = R_d*q_dot_all;
qd_d = R_d*qd_all;
qd_dot_d = R_d*qd_dot_all;
qd_ddot_d = R_d*qd_ddot_all;
qr_dot_d = Tall*qd_dot_d-Lam_all*(q_d-Tall*qd_d);
qr_ddot_d = Tall*qd_ddot_d-Lam_all*(q_dot_d-Tall*qd_dot_d);
s_d = q_dot_d-qr_dot_d;
end
% Find matrices M, C, and G in Lagrangian dynamics
function [M,C,G] = GetMCG(X,q_all,q_dot_all)
global p
M = zeros(3*p,3*p);
C = zeros(3*p,3*p);
G = zeros(3*p,1);
for j = 1:p
    qj = q_all(3*(j-1)+1:3*j);
    dqdtj = q_dot_all(3*(j-1)+1:3*j);
    XXj = [qj;dqdtj;X];
    [Mj,Cj,Gj] = GetMCGj(XXj);
    M(3*(j-1)+1:3*j,3*(j-1)+1:3*j) = Mj;
    C(3*(j-1)+1:3*j,3*(j-1)+1:3*j) = Cj;
    G(3*(j-1)+1:3*j) = Gj;
end
end
% Find matrices Mj, Cj, and Gj for each agent
function [Mj,Cj,Gj] = GetMCGj(XXj)
global mu kJ2 mj
Xj = XXj(1:6);
Xcheaf = XXj(7:12);
xj = Xj(1);
yj = Xj(2);
zj = Xj(3);
r = Xcheaf(1);
vx = Xcheaf(2);
h = Xcheaf(3);
i = Xcheaf(5);
th = Xcheaf(6);
dXcdt = PhaseSyncDynamics(Xcheaf);
dOmdt = dXcdt(4);
didt = dXcdt(5);
dthdt = dXcdt(6);
omx = didt*cos(th)+dOmdt*sin(th)*sin(i);
% check1 omy = -didt*sin(th)+dOmdt*cos(th)*sin(i)
omz = dthdt+dOmdt*cos(i);
% check2 omz-h/r^2
alx = -kJ2*sin(2*i)*cos(th)/r^5+3*vx*kJ2*sin(2*i)*sin(th)/r^4/h...
    -8*kJ2^2*sin(i)^3*cos(i)*sin(th)^2*cos(th)/r^6/h^2;
alz = -2*h*vx/r^3-kJ2*sin(i)^2*sin(2*th)/r^5;
rj = sqrt((r+xj)^2+yj^2+zj^2);
rjZ = (r+xj)*sin(i)*sin(th)+yj*sin(i)*cos(th)+zj*cos(i);
zt = 2*kJ2*sin(i)*sin(th)/r^4;
ztj = 2*kJ2*rjZ/rj^5;
et = sqrt(mu/r^3+kJ2/r^5-5*kJ2*sin(i)^2*sin(th)^2/r^5);
etj = sqrt(mu/rj^3+kJ2/rj^5-5*kJ2*rjZ^2/rj^7);
Mj = mj*eye(3);
Cj = 2*mj*[0 -omz 0;omz 0 -omx;0 omx 0];
Gj1 = -(-xj*(etj^2-omz^2)+yj*alz-zj*omx*omz-(ztj-zt)*sin(i)*sin(th)...
      -r*(etj^2-et^2));
Gj2 = -(-xj*alz-yj*(etj^2-omz^2-omx^2)+zj*alx-(ztj-zt)*sin(i)*cos(th));
Gj3 = -(-xj*omx*omz-yj*alx-zj*(etj^2-omx^2)-(ztj-zt)*cos(i));
Gj = mj*[Gj1;Gj2;Gj3];
end
% Diffusion term of SDE
function G = sde_diffusion(t,X,q_all,q_dot_all)
global d
[M,C,G] = GetMCG(X,q_all,q_dot_all);
n = length(q_all);
G = ones(n,d)*1;
G = M\G;
end
% Nominal controller that guarantees exponential convergence
function tau = ExpControl(X,q_all,q_dot_all,Qd_t)
global Rfd mj L p
R_d = kron(eye(p),Rfd);
[M,C,G] = GetMCG(X,q_all,q_dot_all);
M_d = (R_d*M)/R_d;
C_d = (R_d*C)/R_d;
G_d = R_d*G;
[s_d,qr_dot_d,qr_ddot_d] = Get_s_d(q_all,q_dot_all,Qd_t);
Rfdtau = M_d*qr_ddot_d+C_d*qr_dot_d+G_d-L*s_d;
tau = R_d\Rfdtau;
end
% CVSTEM for spacecraft phase synchronization
function [tau,Qnext,Ks] = PhaseSyncCVSTEM(tau_n,Ks,k,Ufreq,X,q_all,q_dot_all,s_d,Q)
global p
if rem(k-1,Ufreq) == 0
    [Ks,P,dPdt,Q,gam,cvx_status] = CVSTEM(X,q_all,q_dot_all,s_d,Q);
end
if isnan(Ks)
    Q = eye(3*p);
    [Ks,P,dPdt,Qnext,gam,cvx_status] = CVSTEM(X,q_all,q_dot_all,s_d,Q);
    us = -Ks*s_d;
else
    us = -Ks*s_d;
    Qnext = Q;
end
tau = tau_n+us;
end
% CV-STEM
function [K,P,dPdt,Q,gam,cvx_status] = CVSTEM(X,q_all,q_dot_all,s_d,Qprev)
global dt Rfd L R alp mj p umax_tot 
[M,C,G] = GetMCG(X,q_all,q_dot_all);
R_d = kron(eye(p),Rfd);
M_d = (R_d*M)/R_d;
C_d = (R_d*C)/R_d;
A = -M_d\(C_d+L);
Bm = R_d;
B = M_d\Bm;
n = length(q_all);
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
        dQdt = (Q-Qprev)/dt;
        -dQdt+A*Q+Q*A'-nu*B*inv(R)*B'+gam*I <= 0;
        [Q*Bm*inv(R)*B'+B*inv(R)'*Bm'*Q+gam*I+nu*B*inv(R)*B'-2*alp*Q Q;Q 1/2/alp*inv(M_d)*nu] >= 0;
        [tau-n*chi chi;chi nu/trace(M_d)] >= 0;
        nu*norm(inv(R)*B')*norm(s_d) <= lambda_min(Q)*umax_tot;
cvx_end
Pinv = Q/nu;
P = inv(Pinv);
dQdt_tilde = (Q-Qprev)/dt;
dPdt = -P*dQdt_tilde*P/nu;
K = R\B'/Pinv;
gam = gam/nu;
end
% Spacecraft phase sync dynamics
function dXdt = PhaseSyncDynamics(X)
global mu kJ2
r = X(1);
vx = X(2);
h = X(3);
Om = X(4);
i = X(5);
th = X(6);
drdt = vx;
dvxdt = -mu/r^2+h^2/r^3-kJ2/r^4*(1-3*sin(i)^2*sin(th)^2);
dhdt = -kJ2*sin(i)^2*sin(2*th)/r^3;
dOmdt = -2*kJ2*cos(i)*sin(th)^2/h/r^3;
didt = -kJ2*sin(2*i)*sin(2*th)/2/h/r^3;
dthdt = h/r^2+2*kJ2*cos(i)^2*sin(th)^2/h/r^3;
dXdt = [drdt;dvxdt;dhdt;dOmdt;didt;dthdt];
end
% PID controller
function u = pid_controller(t,q_all,q_dot_all,Qd_t)
global dt Tall Rfd Kp Ki Kd p
persistent int_e
if isempty(int_e)
    int_e = zeros(3*p,1);
end
qd = Qd_t(1:3);
qd_dot = Qd_t(4:6);
qd_all = repmat(qd,p,1);
qd_dot_all = repmat(qd_dot,p,1);
R_d = kron(eye(p),Rfd);
qd_d = R_d*qd_all;
qd_dot_d = R_d*qd_dot_all;
e = q_all-R_d\(Tall*qd_d);
e_dot = q_dot_all-R_d\(Tall*qd_dot_d);
int_e = int_e+e*dt;
u = -Kp*e-Ki*int_e-Kd*e_dot;
end
% 4th order Runde-Kutta
function [Xnext] = rk4(X,dt,fun)
k1 = fun(X);
k2 = fun(X+dt*k1/2);
k3 = fun(X+dt*k2/2);
k4 = fun(X+dt*k3);
Xnext = X+dt/6*(k1+2*k2+2*k3+k4);
end
% Dynamics in LVLH frame
function dXXdt = RelDynamics(XX,tau,mj,p)
global mj p
q_all = XX(1:3*p);
q_dot_all = XX(3*p+1:6*p);
X = XX(6*p+1:6*p+6);
[M,C,G] = GetMCG(X,q_all,q_dot_all);
q_ddot_all = -M\(C*q_dot_all+G-tau);
dXdt = PhaseSyncDynamics(X);
dXXdt = [q_dot_all;q_ddot_all;dXdt];
end
% H infinity controller
function [tau,Q] = Hinfinity(X,q_all,q_dot_all,Qd_t,Qprev,gam)
global dt Rfd p R Tall
qd_val = Qd_t(1:3);
qd_dot_val = Qd_t(4:6);
R_d = kron(eye(p),Rfd);
qd_all = repmat(qd_val,p,1);
qd_dot_all = repmat(qd_dot_val,p,1);
qd_all_trans = R_d\(Tall*R_d*qd_all);
qd_dot_all_trans = R_d\(Tall*R_d*qd_dot_all);
qq = [q_all;q_dot_all];
qqd = [qd_all_trans;qd_dot_all_trans];           
A = AmatSDC(X,q_all,q_dot_all);
B = BmatSDC(X,q_all,q_dot_all);
n = length(q_all)*2;
I = eye(n);
cvx_begin sdp quiet
    variable Q(n,n) semidefinite
    variable g_inv nonnegative
    minimize (g_inv)
    subject to
        dQdt = (Q-Qprev)/(dt);
        [-dQdt+A*Q+Q*A'-B*inv(R)*B' Q;Q -I*g_inv] <= 0;
        g_inv >= 1/gam;
cvx_end
if g_inv >= 0.7
    flag = 1;
    cvx_begin sdp quiet
        variable Q(n,n) semidefinite
        variable g_inv nonnegative
        minimize (g_inv)
        subject to
            [A*Q+Q*A'-B*inv(R)*B' Q;Q -I*g_inv] <= 0;
            g_inv >= 1/gam;
    cvx_end
else
    flag = 0;
end
Pinv = Q;
K = inv(R)*B'/Pinv;
tau = -K*(qq-qqd);
end
% A matrix of SDC formulation 
function A = AmatSDC(X,q_all,q_dot_all)
global p
[M,C,G] = GetMCG(X,q_all,q_dot_all);
Gmat = GetGmat(X,q_all,q_dot_all);
A11 = zeros(3*p,3*p);
A12 = eye(3*p);
A21 = -M\Gmat;
A22 = -M\C;
A = [A11 A12;A21 A22];
end
% B matrix of SDC formulation 
function B = BmatSDC(X,q_all,q_dot_all)
global p
[M,C,G] = GetMCG(X,q_all,q_dot_all);
B = [zeros(3*p,3*p);inv(M)];
end
% Compute whole states of spacecraft phase synchronization
function [XX,qq,qqd] = GetXX(X,q_all,q_dot_all,Qd_t)
global Rfd Tall p
qd_val = Qd_t(1:3);
qd_dot_val = Qd_t(4:6);
R_d = kron(eye(p),Rfd);
qd_all = repmat(qd_val,p,1);
qd_dot_all = repmat(qd_dot_val,p,1);
qd_all_trans = R_d\(Tall*R_d*qd_all);
qd_dot_all_trans = R_d\(Tall*R_d*qd_dot_all);
qq = [q_all;q_dot_all];
qqd = [qd_all_trans;qd_dot_all_trans];
XX = [q_all;q_dot_all;X];
end
% Compute Gmat for SDC formulation
function Gmat = GetGmat(X,q_all,q_dot_all)
global p
Gmat = zeros(3*p,3*p);
for j = 1:p
    qj = q_all(3*(j-1)+1:3*j);
    dqdtj = q_dot_all(3*(j-1)+1:3*j);
    XXj = [qj;dqdtj;X];
    Gjmat = GetGmatj(XXj);
    Gmat(3*(j-1)+1:3*j,3*(j-1)+1:3*j) = Gjmat;
end
end
% Compute Gmat for each agent
function Gjmat = GetGmatj(XXj)
global mu kJ2 mj
Xj = XXj(1:6);
Xcheaf = XXj(7:12);
xj = Xj(1);
yj = Xj(2);
zj = Xj(3);
r = Xcheaf(1);
vx = Xcheaf(2);
h = Xcheaf(3);
i = Xcheaf(5);
th = Xcheaf(6);
dXcdt = PhaseSyncDynamics(Xcheaf);
dOmdt = dXcdt(4);
didt = dXcdt(5);
dthdt = dXcdt(6);
omx = didt*cos(th)+dOmdt*sin(th)*sin(i);
omz = dthdt+dOmdt*cos(i);
alx = -kJ2*sin(2*i)*cos(th)/r^5+3*vx*kJ2*sin(2*i)*sin(th)/r^4/h...
    -8*kJ2^2*sin(i)^3*cos(i)*sin(th)^2*cos(th)/r^6/h^2;
alz = -2*h*vx/r^3-kJ2*sin(i)^2*sin(2*th)/r^5;
rj = sqrt((r+xj)^2+yj^2+zj^2);
rjZ = (r+xj)*sin(i)*sin(th)+yj*sin(i)*cos(th)+zj*cos(i);
zt = 2*kJ2*sin(i)*sin(th)/r^4;
ztj = 2*kJ2*rjZ/rj^5;
et = sqrt(mu/r^3+kJ2/r^5-5*kJ2*sin(i)^2*sin(th)^2/r^5);
etj = sqrt(mu/rj^3+kJ2/rj^5-5*kJ2*rjZ^2/rj^7);
Remjmat = [(ztj-zt)*sin(i)*sin(th)+r*(etj^2-et^2);(ztj-zt)*sin(i)*cos(th);(ztj-zt)*cos(i)];
Gjmat = mj*[(etj^2-omz^2)+Remjmat(1)/xj -alz omx*omz;alz (etj^2-omz^2-omx^2)+Remjmat(2)/yj -alx;omx*omz alx (etj^2-omx^2)+Remjmat(3)/zj];
end