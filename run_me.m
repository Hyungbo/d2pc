%%% 
%%% Simulation for D2PC (Data-Driven Predictive Control)
%%% including comparisons with DeePC and MPC
%%%
%%% (Parameters are the same as those used in the paper.)
%%%
%%% Nam H. Jo and H. Shim
%%% July 22, 2021
%%%
%%% tested in MATLAB 2020b
%%% with OSQP (https://osqp.org); a QP solver 
%%%
%%% can handle MIMO through MISO approach
%%%
%%% the cost function to be minimized: 
%%% \sum_{k=0}^{N-1} ||y_k - r||^2 + ||u_k||^2
%%%

%%% OSQP or CVX is REQUIRED to solve the control input under constraints.
%%% Specify your solver option around line 76.
%%% If you don't have OSQP or CVX, choose the solver option = 'handful',
%%% but the control problem is solved without constraints
%%% (so that the outcomes are different from the paper).

%%% For more testing, change the model below and 
%%% changethe noise level An by yourself
%%% at the end of this file where the plant is described.

clear

% Select your model; uncomment one of three, 
% or, you can define your own at the end of this file in f_DefinePlant.
model = 'InvPen';   
%model = 'TwoCart';
%model = 'FourTank';

switch model
    case 'InvPen'
        % Get the discrete-time real plant
        P_real = f_DefinePlant('InvPen');
        
        OP.umax = 20;   % size: (nu by 1) Input constraint
        OP.umin = -20;  % size: (nu by 1) Input constraint
        OP.ymax = 100;  % size: (ny by 1) Output constraint
        OP.ymin = -100; % size: (ny by 1) Output constraint
        OP.Q = 1000;    % size: (ny by ny)
        OP.R = 1;       % size: (nu by nu)
        OP.r = 1;       % target set-point of the output (ny by 1)
        OP.N = 20;      % prediction horizon (step) of cost function (from 0 to N-1)
        
    case 'TwoCart'
        P_real = f_DefinePlant('TwoCart');
        
        OP.umax = 2;    % size: (nu by 1) Input constraint
        OP.umin = -2;   % size: (nu by 1) Input constraint
        OP.ymax = Inf;  % size: (ny by 1) Output constraint
        OP.ymin = -Inf; % size: (ny by 1) Output constraint
        OP.Q = 200;     % size: (ny by ny)
        OP.R = 1;       % size: (nu by nu)
        OP.r = 1;       % target set-point of the output (ny by 1)
        OP.N = 20;      % prediction horizon (step) of cost function (from 0 to N-1)
        
    case 'FourTank'
        P_real = f_DefinePlant('FourTank');
        
        OP.umax = [Inf; Inf];   % size: (nu by 1) Input constraint
        OP.umin = -OP.umax;     % size: (nu by 1) Input constraint
        OP.ymax = [Inf; Inf];   % size: (ny by 1) Output constraint
        OP.ymin = -OP.ymax;     % size: (ny by 1) Output constraint
        OP.Q = 3*eye(2);        % size: (ny by ny)
        OP.R = 0.01*eye(2);     % size: (nu by nu)
        OP.r = [0.65; 0.77];    % target set-point of the output (ny by 1)
        OP.N = 30;      % prediction horizon (step) of cost function (from 0 to N-1)
end

OP.solver = 'osqp'; % solver = 'osqp' or 'cvx' or 'handful'
%OP.solver = 'handful';
OP.P_qp = [];   % to define QP problem
OP.q_qp = [];   % to define QP problem
OP.A_qp = [];   % to define QP problem
OP.l_qp = [];   % to define QP problem
OP.u_qp = [];   % to define QP problem
OP.prob = [];   % for use of osqp
OP.ApBp = [];   % for use of d2pc
OP.q_qp0 = [];  % for internal use of DeePC
OP.q_qp1 = [];  % for internal use of DeePC
OP.Uf1 = [];    % for internal use of DeePC

% for Simulation
Simulation.nsim = 60;   % simulation step
Trial.no = 10;          % Number of repeated trials

% Parameters for d2pc
d2pc.NBar = 4;          % Estimated order of the plant
d2pc.Nd = 1;            % # of averaging for Ad and Bd
d2pc.Tsampling = 100;   % length of sampling input (step)
d2pc.gen_sampling_input = @(nu,T) 2*rand(nu,T)-1;   % generating function of sampling input 

% Parameters for DeePC
deepc.Tini = 4;     % Tini in the paper
deepc.Tsampling = 100;   % length of sampling input
deepc.gen_sampling_input = @(nu,T) 2*rand(nu,T)-1;  % sampling input
deepc.q = 1;        % number of datasets
deepc.Nd = 1;       % number of averaging Hankel matrix
deepc.lambda_g = 0;
deepc.lambda_y = 0; % 0 or non-zero determines regularization

switch model
    case 'InvPen'
        d2pc.Tsampling = 17;
        
        deepc.q = 10;
        kPE = 2*deepc.Tini+OP.N;
        TdeepcM = (kPE-1) + ceil(kPE/deepc.q);  % Tdeepc Multi data set
        deepc.Tsampling = TdeepcM;
        
        % Below: for Tables 2 and 3 in the paper
        % d2pc.NBar = 14; d2pc.Nd = 50; d2pc.Tsampling = 55
        % deepc.lambda_g = 500; deepc.lambda_y = 5e5;
    case 'TwoCart'
        d2pc.NBar = 20;   

        deepc.Tini = 15;
        deepc.lambda_g = 500;        deepc.lambda_y = 5e5;

        % Below: for Table 6 in the paper
        % d2pc.Nd = 50;      deepc.Nd = d2pc.Nd;
    case 'FourTank'
        d2pc.Tsampling = 400;  
        d2pc.NBar=30; 

        deepc.Tsampling = d2pc.Tsampling;
        deepc.Tini = 30;
        deepc.lambda_g = 0.1;        deepc.lambda_y = 1000;
        % Below: for Table 9 in the paper
        % d2pc.Nd = 50;        deepc.Nd = d2pc.Nd;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate: Model-based MPC w/o noise to get the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default');     % reset random number generator
Pideal = P_real; Pideal.An = 0;
OP = f_MPC_setup(OP,Pideal);

% Simulate the closed-loop
Xmpc = [];  % for state recording
Umpc = [];  % for input recording
x_cur = Pideal.x0; % current state of the plant
for i = 1 : Simulation.nsim
    
    % measure the current state x_cur & record
    x_cur_measured = Pideal.gen_noise(x_cur, Pideal.An);
    Xmpc = [Xmpc, x_cur];
    
    % compute feedback input u0 & record
    fprintf('=== MPC: step %d out of %d === \n', i, Simulation.nsim);
    [u0, err_code] = f_MPC_solve(OP,Pideal,x_cur_measured);
    if err_code == 1, error('! MPC failure'); end
    Umpc = [Umpc, u0];
    
    % advance the plant to the next x_cur
    x_cur = Pideal.A*x_cur + Pideal.B*u0;
end
Trial.Xmpc = Xmpc;
% temporary plot of simulation result
% v_drawSIM(Xmpc,Umpc,Simulation.nsim,Pideal.Ts)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate: D2PC 
%%% Most of the codes in this block is for displaying.
%%% Key functions are performed in SimulateD2PC function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default');     % reset random number generator

plot_time = (0:Simulation.nsim-1) * P_real.Ts;
figure(1);clf
subplot(2,1,1)
plot(plot_time, P_real.C*Xmpc, 'r--','Linewidth',2);
hold on
subplot(2,1,2)
plot(plot_time, Umpc, 'r--','Linewidth',2);
hold on
MAEtable = [];
for i = 1:Trial.no
   Trial.i = i;
   
   %%% Main function %%%
   [MAE, Xd2pc, Ud2pc, err_code] = SimulateD2PC(P_real,d2pc,OP,Simulation,Trial);
   
   if err_code == 0
        MAEtable = [MAEtable, MAE];
        figure(1)
        subplot(2,1,1)
        plot(plot_time, P_real.C*Xd2pc);
        subplot(2,1,2)
        plot(plot_time, Ud2pc);
        figure(3)
        subplot(221);plot(plot_time, P_real.C*Xd2pc,'b');ylabel('y');hold on
        grid on; xlabel('Time [sec]')
   end
end
FailureRatio = 1-length(MAEtable)/Trial.no;
MAE = mean(MAEtable);
figure(1)
subplot(2,1,1), title('Output')
subplot(2,1,2), title('Input')
text = sprintf('D2PC: NBar = %d, An = %g, Nd = %d, FR = %1.2f, MAE = %1.3f', ...
    d2pc.NBar,P_real.An,d2pc.Nd,FailureRatio,MAE);
xlabel(text)
figure(3); 
subplot(221);title('D2PC (Proposed)');
if FailureRatio==1
    title('D2PC: FailureRatio = 1')
else
    plot(plot_time, P_real.C*Xmpc,'r--','Linewidth',2);hold off
end
MAEd2pc = MAE   % to display


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate: DeePC
%%% Most of the codes in this block is for displaying.
%%% Key functions are performed in SimulateDeePC function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(OP.solver, 'handful')
    disp('DeePC is not simulated without a solver, and the code stops.')
    return
end

rng('default');     % reset random number generator

plot_time = (0:Simulation.nsim-1) * P_real.Ts;
figure(2);clf
subplot(2,1,1)
plot(plot_time, P_real.C*Xmpc, 'r--','Linewidth',2);
hold on
subplot(2,1,2)
plot(plot_time, Umpc, 'r--','Linewidth',2);
hold on
MAEtable = [];
for i = 1:Trial.no
   Trial.i = i;
    
   %%% Main function %%%
   [MAE, Xdeepc, Udeepc, err_code] = SimulateDeePC(P_real,deepc,OP,Simulation,Trial);
   
   if err_code == 0
        MAEtable = [MAEtable, MAE];
        figure(2)
        subplot(2,1,1)
        plot(plot_time, P_real.C*Xdeepc);
        subplot(2,1,2)
        plot(plot_time, Udeepc);
        figure(3)
        subplot(222);plot(plot_time, P_real.C*Xdeepc,'b');ylabel('y');hold on
        grid on; xlabel('Time [sec]')
   end
end
FailureRatio = 1-length(MAEtable)/Trial.no;
MAE = mean(MAEtable);
figure(2)
subplot(2,1,1), title('Output')
subplot(2,1,2), title('Input')
text = sprintf('DeePC: lam_y = %g, Tini = %d, An = %g, Nd = %d, q = %d, FR = %1.2f, MAE = %1.3f', ...
    deepc.lambda_y,deepc.Tini,P_real.An,deepc.Nd,deepc.q,FailureRatio,MAE);
xlabel(text)
figure(3);
subplot(222);title('DeePC');
if FailureRatio==1
    title('DeePC: FailureRatio = 1')
else
    plot(plot_time, P_real.C*Xmpc,'r--','Linewidth',2);hold off
end
MAEd2pc         % to display
MAEdeepc = MAE  % to display

return



%%%%%%%%%%%%%%%%%%%%
%%% Function Library
%%%%%%%%%%%%%%%%%%%%

function v_drawSIM(X,U,nsim,Ts)

plot_time = (0:nsim-1) * Ts;
figure
subplot(2,1,1)
plot(plot_time, X);
subplot(2,1,2)
plot(plot_time, U);

end    

function OP = f_MPC_setup(OP,plant)
%%%
%%% Cast MPC problem to a QP
%%% reference: https://osqp.org (MPC example)
%%%
%%% minimize (1/2)*x'*P_qp*x + q_qp'*x  
%%% where x = (x(0), ..., x(N-1), u(0), ... , u(N-1))
%%%
%%% This function updates OP.P_qp, OP.q_qp, OP.A_qp, OP.l_qp, OP.u_qp
%%%

N = OP.N;
Q = OP.Q;
R = OP.R;
r = OP.r;
A = plant.A;
B = plant.B;
C = plant.C;
D = plant.D;
nx = plant.nx;
nu = plant.nu;
x0 = plant.x0;


if strcmp(OP.solver,'handful')
    Abig = eye(nx);
    Bbig = zeros(nx, nu*N);
    for i = 2:N
        Abig = [Abig; A*Abig(end-nx+1:end,:)];
        Bbig = [Bbig; A*Bbig(end-nx+1:end,:)];
        Bbig((i-1)*nx+1:i*nx,(i-2)*nu+1:(i-1)*nu) = B;
    end
    OP.S0 = inv(Bbig'*kron(speye(N), C'*Q*C)*Bbig + kron(speye(N), R)) ...
        * Bbig'*kron(speye(N), C'*Q*C)*Abig;
    OP.S1 = inv(Bbig'*kron(speye(N), C'*Q*C)*Bbig + kron(speye(N), R)) ...
        * Bbig'*kron(speye(N), C'*Q)*kron(ones(N,1), r);
else
    OP.P_qp = [kron(speye(N), C'*Q*C), kron(speye(N), C'*Q*D);
        kron(speye(N), D'*Q*C), kron(speye(N), R + D'*Q*D) ];
    OP.q_qp = [repmat(-C'*Q*r, N, 1); repmat(-D'*Q*r, N, 1)];
    % with linear dynamics
    Ax = kron(speye(N), -speye(nx)) + kron(sparse(diag(ones(N-1, 1), -1)), A);
    Bu = kron(sparse(diag(ones(N-1, 1), -1)), B);
    % to have Aeq*x = leq, let leq <= Aeq*x <= leq
    Aeq = [Ax, Bu];
    leq = [-x0; zeros((N-1)*nx, 1)];  % the current state enters here
    ueq = leq;
    % to have lineq .<= Aineq*x .<= uineq
    Aineq = [kron(speye(N), sparse(C)), kron(speye(N), sparse(D));
        sparse(zeros(nu*N,nx*N)), speye(nu*N) ];
    lineq = [repmat(OP.ymin, N, 1); repmat(OP.umin, N, 1)];
    uineq = [repmat(OP.ymax, N, 1); repmat(OP.umax, N, 1)];
    % translate to l_qp .<= A_qp .<= u_qp
    OP.A_qp = [Aeq; Aineq];
    OP.l_qp = [leq; lineq];
    OP.u_qp = [ueq; uineq];
end

if strcmp(OP.solver,'osqp')
    % Create an OSQP object
    OP.prob = osqp;
    
    % Setup workspace (parameters can be adjusted)
    OP.prob.setup(OP.P_qp, OP.q_qp, OP.A_qp, OP.l_qp, OP.u_qp, ...
        'polish', 1, 'eps_rel',1e-4, 'eps_abs',1e-4,...
        'max_iter',200000, 'warm_start', true,  'verbose', false);
end

end

function [u0, err_code] = f_MPC_solve(OP,plant,x_cur)

nx = plant.nx;
nu = plant.nu;
N = OP.N;
l_qp = OP.l_qp;
u_qp = OP.u_qp;

% Update the constraints
l_qp(1:nx) = -x_cur;
u_qp(1:nx) = -x_cur;
    
if strcmp(OP.solver,'osqp')
    OP.prob.update('l', l_qp, 'u', u_qp);
end

% Solve
switch OP.solver
    case 'osqp'
        res = OP.prob.solve();
        if ~strcmp(res.info.status, 'solved')
            disp('OSQP did not solve the problem!')
            err_code = 1;
            return;
        end
        u0 = res.x(N*nx+1:N*nx+nu);

    case 'cvx'
        cvx_begin quiet
            variable x( N*nx + N*nu )
            minimize ( (1/2)*quad_form(x,OP.P_qp) + OP.q_qp'*x )
            OP.A_qp*x >= l_qp;
            OP.A_qp*x <= u_qp;
        cvx_end
        if ~strcmp(cvx_status, 'Solved')
            disp('CVX did not solve the problem!')
            err_code = 1;
            return;
        end
        u0 = x(N*nx+1:N*nx+nu);

    case 'handful'
        x = -OP.S0*x_cur + OP.S1;
        u0 = x(1:nu);
        
    otherwise
        error('Option simulation.solver is incorrect.')
end
err_code = 0;

end

function [CMat,DMat] = f_D2PC_sampling(plant, d2pc)
% MISO implementation
% Averaging by Nd implemented

NBar = d2pc.NBar;
Nd = d2pc.Nd;
T = d2pc.Tsampling;
ny = plant.ny;
nx = plant.nx;
nu = plant.nu;

if NBar < nx
    error("Your guess of plant's order is less than the actual.")
end

% placeholder (that will be averaged) for [CMat,DMat] of output j
CDMat = zeros(ny, NBar + nu*NBar + nu);
for k = 1:Nd
    
    % Generate sampling input (should be PE); (nu by Tsampling)
    U = d2pc.gen_sampling_input(nu,T);
    
    % Get the output of the plant
    x = plant.x0;
    Y = [];
    for i = 1:T
        Y = [Y, plant.C * x + plant.D * U(:,i)];
        x = plant.A * x + plant.B * U(:,i);
    end
    Y = plant.gen_noise(Y, plant.An);   % Y is now corruped with noise

    nrow = NBar+nu*NBar; ncol = T-NBar+1;  
    CDMat_ = zeros(ny, NBar + nu*NBar + nu);
    for j = 1:ny
        % Get Ad and Bd and accumulate for output j
        XX = zeros(nrow,ncol);
        for i=1:ncol
            XX(1:NBar,i) = Y(j, i:i+NBar-1); 
            u_ = U(:, i:i+NBar-1);
            XX(NBar+1:end,i) = u_(:);
        end
        X0 = XX(:,1:end-1);
        X1 = XX(:,2:end);
        XU = [X0; U(:, end-(ncol-1)+1:end)]; % (NBar + nu*NBar + 1) by (T - NBar)

        % Taking pseudo-inverse of XU matrix
        AdBd = X1*pinv(XU);
        CDMat_(j,:) = AdBd(NBar,:);
    end
    CDMat = CDMat + CDMat_; % accumulating
end
CDMat = CDMat/Nd;   % taking average
CMat = CDMat(:,1:end-nu);
DMat = CDMat(:,end-nu+1:end);
end

function OP = f_D2PC_setup(OP,d2pc,CMat,DMat)
%%%
%%% Cast D2PC problem to a QP
%%%
%%% minimize (1/2)*x'*P_qp*x + q_qp'*x  
%%% where decision variable x = (y(0), ..., y(N-1), u(0), ... , u(N-1))
%%%
%%% This function updates OP.P_qp, OP.q_qp, OP.A_qp, OP.l_qp, OP.u_qp
%%%

N = OP.N;
Q = OP.Q;
R = OP.R;
r = OP.r;
NBar = d2pc.NBar;
%CMat = d2pc.CMat;
%DMat = d2pc.DMat;
[ny,nu] = size(DMat);

OP.P_qp = blkdiag( kron(speye(N), Q), kron(speye(N), R) );
OP.q_qp = [repmat(-Q*r, N, 1); zeros(N*nu, 1)];
Af = kron(speye(N),-eye(ny));
for i = 1:min(N-1,NBar)
    Af = Af + kron( sparse(diag( ones(N-i,1), -i )), ...
        diag( CMat(:,NBar-i+1) ) );
end
Bf = kron(speye(N),DMat);
for i = 1:min(N-1,NBar)
    Bf = Bf + kron( sparse(diag( ones(N-i,1), -i )), ...
        CMat(:,NBar+nu*(NBar-i)+1:NBar+nu*(NBar-i+1)) );
end
Ap = sparse(zeros(ny*NBar,ny*NBar));
Bp = sparse(zeros(ny*NBar,nu*NBar));
for i = 1:NBar
    Ap = Ap + kron( sparse(diag( ones(NBar-i+1,1), i-1 )), ...
        diag( CMat(:,i) ) );
    Bp = Bp + kron( sparse(diag( ones(NBar-i+1,1), i-1 )), ...
        CMat(:, NBar+nu*(i-1)+1:NBar+nu*i) );
end
OP.ApBp = sparse(zeros(ny*N, ny*NBar+nu*NBar));
OP.ApBp(1:ny*min(NBar,N),:) = [Ap(1:ny*min(NBar,N),:), Bp(1:ny*min(NBar,N),:)];

% to have Aeq*x = leq, let leq <= Aeq*x <= leq
Aeq = [Af, Bf];
leq = zeros(ny*N,1);  % the current state enters here
ueq = leq;
% to have lineq .<= Aineq*x .<= uineq
Aineq = speye(N*(ny+nu));
lineq = [repmat(OP.ymin, N, 1); repmat(OP.umin, N, 1)];
uineq = [repmat(OP.ymax, N, 1); repmat(OP.umax, N, 1)];
% translate to l_qp .<= A_qp .<= u_qp
OP.A_qp = [Aeq; Aineq];
OP.l_qp = [leq; lineq];
OP.u_qp = [ueq; uineq];

if strcmp(OP.solver,'osqp')
    % Create an OSQP object
    OP.prob = osqp;
    
    % Setup workspace (parameters can be adjusted)
    OP.prob.setup(OP.P_qp, OP.q_qp, OP.A_qp, OP.l_qp, OP.u_qp, ...
        'polish', 1, 'eps_rel',1e-4, 'eps_abs',1e-4,...
        'max_iter',200000, 'warm_start', true,  'verbose', false);
elseif strcmp(OP.solver,'handful')
    % to make optimal solution: x = S0 + S1*leq
    % with no constraints => linear feedback control
    iP = inv(OP.P_qp);
    iAPA = inv(Aeq*iP*Aeq');
    OP.S0 = (iP*Aeq'*iAPA*Aeq*iP - iP)*OP.q_qp;
    OP.S1 = iP*Aeq'*iAPA;
end

end

function [u0, err_code] = f_D2PC_solve(OP,chi_cur)

ApBp = OP.ApBp;
nu = size(OP.R,1);
ny = size(OP.Q,1);
N = OP.N;
l_qp = OP.l_qp;
u_qp = OP.u_qp;

% Update the constraints
l_qp(1:ny*N) = -ApBp*chi_cur;
u_qp(1:ny*N) = -ApBp*chi_cur;
    
if strcmp(OP.solver,'osqp')
    OP.prob.update('l', l_qp, 'u', u_qp);
end

% Solve
switch OP.solver
    case 'osqp'
        res = OP.prob.solve();
        if ~strcmp(res.info.status, 'solved')
            disp('OSQP did not solve the problem!')
            u0 = [];
            err_code = 1;
            return;
        end
        u0 = res.x(ny*N+1:ny*N+nu);

    case 'cvx'
        cvx_begin quiet
            variable x( N*nx + N*nu )
            minimize ( (1/2)*quad_form(x,OP.P_qp) + OP.q_qp'*x )
            OP.A_qp*x >= l_qp;
            OP.A_qp*x <= u_qp;
        cvx_end
        if ~strcmp(cvx_status, 'Solved')
            disp('CVX did not solve the problem!')
            err_code = 1;
            return;
        end
        u0 = x(ny*N+1:ny*N+nu);
        
    case 'handful'
        x = OP.S0 + OP.S1*l_qp(1:ny*N);
        u0 = x(ny*N+1:ny*N+nu);

    otherwise
        error('Option simulation.solver is incorrect.')
end
err_code = 0;

end

function [MAE, Xd2pc, Ud2pc, err_code] = SimulateD2PC(P_real,d2pc,OP,Simulation,Trial)

%%% Offline sampling of real plant (data-driven model id)
[CMat,DMat] = f_D2PC_sampling(P_real, d2pc);
%d2pc.CMat = CMat;
%d2pc.DMat = DMat; 
%d2pc.DMat = zeros(size(DMat));    % Remove this line later!

%%% Setup D2PC
OP = f_D2PC_setup(OP,d2pc,CMat,DMat);

%%% Simulate the closed-loop

% free-flight of the plant to obtain initial chi(0)
[y0, x_cur] = f_FreeFlight(P_real,d2pc.NBar);
chi_cur = [y0; zeros(P_real.nu*d2pc.NBar,1)]; % initial condition for D2PC

Xd2pc = [];  % for state recording
Ud2pc = [];  % for input recording
for i = 1 : Simulation.nsim
    
    % record the state of real plant
    Xd2pc = [Xd2pc, x_cur];
    
    % compute feedback input u0 & record
    fprintf('=== D2PC: trial %d, step %d out of %d === \n', Trial.i, i, Simulation.nsim);
    [u0, err_code] = f_D2PC_solve(OP,chi_cur);
    if err_code == 1, disp('! D2PC failure'); MAE = []; return; end
    Ud2pc = [Ud2pc, u0];
    
    % get the output of the plant and advance the plant to the next x_cur
    y_cur = P_real.C*x_cur + P_real.D*u0;
    x_cur = P_real.A*x_cur + P_real.B*u0;
    
    % make next chi_cur with noise
    chi_cur = [chi_cur(P_real.ny+1:P_real.ny*d2pc.NBar); ...
        P_real.gen_noise(y_cur,P_real.An); ...
        chi_cur(P_real.ny*d2pc.NBar+P_real.nu+1:end); u0];
end
% temporary plot of simulation result
%v_drawSIM(Xd2pc,Ud2pc,Simulation.nsim,P_real.Ts)

MAEsum = 0;
Y_err = P_real.C*(Xd2pc - Trial.Xmpc);
for i=1:Simulation.nsim
    MAEsum = MAEsum + norm( Y_err(:,i), 2);
end
MAE = MAEsum/Simulation.nsim;
err_code = 0;

end

function [Up,Yp,Uf,Yf] = f_DeePC_sampling(plant,deepc,OP)

Nd = deepc.Nd;
q = deepc.q;
T = deepc.Tsampling;
Tini = deepc.Tini;
N = OP.N;
ny = plant.ny;
nx = plant.nx;
nu = plant.nu;

UY = zeros( (nu+ny)*(Tini+N), q*(T-(Tini+N)+1) );
for k = 1:Nd
    UY_ = [];
    for i = 1:q
        UY_ = [UY_, f_GetUY(plant,T,Tini,N,deepc.gen_sampling_input)];
    end
    UY = UY + UY_;
end
UY = UY/Nd; % Now UY is ready
Up = UY(1:nu*Tini, :);
Uf = UY(nu*Tini+1:nu*(Tini+N), :);
Yp = UY(nu*(Tini+N)+1:nu*(Tini+N)+ny*Tini, :);
Yf = UY(nu*(Tini+N)+ny*Tini+1:end, :);

end

function OP = f_DeePC_setup(OP,deepc,Up,Yp,Uf,Yf)
%%%
%%% Cast DeePC problem to a QP
%%%
%%% minimize (1/2)*x'*P_qp*x + q_qp'*x  
%%% where x = (g(1),g(2),...,g(ng))
%%%

N = OP.N;
Q = OP.Q;
R = OP.R;
r = OP.r;
Tini = deepc.Tini;
ny = size(Q,1);
nu = size(R,1);
lambda_g = deepc.lambda_g;
lambda_y = deepc.lambda_y;

ng = size(Up,2);    % size of g

OP.P_qp = blkdiag( Yf'*kron(speye(N), Q)*Yf + Uf'*kron(speye(N), R)*Uf ) ...
    + sparse(lambda_g*eye(ng) + lambda_y*Yp'*Yp);
OP.q_qp0 = -Yf' * kron(speye(N), Q) * repmat(r,N,1);
OP.q_qp1 = lambda_y * Yp'; 
OP.q_qp = OP.q_qp0; % this updates each time as q_qp = q_qp0 - q_qp1*y_ini

if lambda_y == 0    % Regularization is OFF
    Aeq = [Up; Yp];
    leq = zeros( (nu+ny)*Tini, 1 ); % placeholder for [u_ini; y_ini]    
    ueq = leq;
else % Regularization is ON
    Aeq = Up;
    leq = zeros( nu*Tini, 1 ); % placeholder for [u_ini]
    ueq = leq;
end

% input and state constraints
Aineq = [Uf; Yf];
lineq = [repmat(OP.umin, N, 1); repmat(OP.ymin, N, 1)];
uineq = [repmat(OP.umax, N, 1); repmat(OP.ymax, N, 1)];

% final constraints
OP.A_qp = [Aeq; Aineq];
OP.l_qp = [leq; lineq];
OP.u_qp = [ueq; uineq];

OP.Uf1 = Uf(1:nu,:);    % is used when u0 is computed after optimization

if strcmp(OP.solver,'osqp')
    % Create an OSQP object
    OP.prob = osqp;
    
    % Setup workspace (parameters can be adjusted)
    OP.prob.setup(OP.P_qp, OP.q_qp, OP.A_qp, OP.l_qp, OP.u_qp, ...
        'polish', 1, 'eps_rel',1e-4, 'eps_abs',1e-4,...
        'max_iter',200000, 'warm_start', true,  'verbose', false);

elseif strcmp(OP.solver,'handful')
    % to make optimal solution: x = S0 + S1*leq
    % with no constraints => linear feedback control
    rank(OP.P_qp)
    size(OP.P_qp)
    %pause
    iP = inv(OP.P_qp);
    iAPA = inv(Aeq*iP*Aeq');
    OP.S0 = (iP*Aeq'*iAPA*Aeq*iP - iP)*OP.q_qp0;
    OP.S1 = (iP*Aeq'*iAPA*Aeq*iP - iP)*OP.q_qp1;
    OP.S2 = iP*Aeq'*iAPA;
end

end

function [u0, err_code] = f_DeePC_solve(OP,deepc,u_ini,y_ini)

q_qp0 = OP.q_qp0;
q_qp1 = OP.q_qp1;
nu = size(OP.R,1);
ny = size(OP.Q,1);
l_qp = OP.l_qp;
u_qp = OP.u_qp;
N = OP.N;
Tini = deepc.Tini;
lambda_y = deepc.lambda_y;

% Update the constraints
q_qp = q_qp0 - q_qp1*y_ini;
if lambda_y == 0    % Regularization is OFF
    l_qp(1:(nu+ny)*Tini) = [u_ini; y_ini];
    u_qp(1:(nu+ny)*Tini) = [u_ini; y_ini];
    l_handful = [u_ini; y_ini]; % for handful computation
else                % Regularization is ON
    l_qp(1:nu*Tini) = [u_ini];
    u_qp(1:nu*Tini) = [u_ini];
    l_handful = [u_ini]; % for handful computation
end
    
if strcmp(OP.solver,'osqp')
    OP.prob.update('q', q_qp, 'l', l_qp, 'u', u_qp);
end

% Solve
switch OP.solver
    case 'osqp'
        res = OP.prob.solve();
        if ~strcmp(res.info.status, 'solved')
            disp('OSQP did not solve the problem!')
            u0 = [];
            err_code = 1;
            return;
        end
        u0 = OP.Uf1 * res.x;

    case 'cvx'
        cvx_begin quiet
            variable x( N*nx + N*nu )
            minimize ( (1/2)*quad_form(x,OP.P_qp) + OP.q_qp'*x )
            OP.A_qp*x >= l_qp;
            OP.A_qp*x <= u_qp;
        cvx_end
        if ~strcmp(cvx_status, 'Solved')
            disp('CVX did not solve the problem!')
            err_code = 1;
            return;
        end
        u0 = OP.Uf1 * x;
        
    case 'handful'
        x = OP.S0 - OP.S1*y_ini + OP.S2*l_handful;
        u0 = OP.Uf1 * x;

    otherwise
        error('Option simulation.solver is incorrect.')
end
err_code = 0;

end

function [MAE, Xdeepc, Udeepc, err_code] = SimulateDeePC(P_real,deepc,OP,Simulation,Trial)

%%% Offline sampling of real plant (data-driven model id)
[Up,Yp,Uf,Yf] = f_DeePC_sampling(P_real, deepc, OP);

%%% Setup DeePC
OP = f_DeePC_setup(OP,deepc,Up,Yp,Uf,Yf);

% Simulate the closed-loop
u_ini = zeros(P_real.nu*deepc.Tini, 1);
[y_ini, x_cur] = f_FreeFlight(P_real,deepc.Tini);

Xdeepc = [];  % for state recording
Udeepc = [];  % for input recording
for i = 1 : Simulation.nsim
    
    % record the state of real plant
    Xdeepc = [Xdeepc, x_cur];
    
    % compute feedback input u0 & record
    fprintf('=== DeePC: trial %d, step %d out of %d === \n', Trial.i, i, Simulation.nsim);
    [u0, err_code] = f_DeePC_solve(OP,deepc,u_ini,y_ini);
    if err_code == 1, disp('! DeePC failure'); MAE = []; return; end
    Udeepc = [Udeepc, u0];
    
    % get the output of the plant and advance the plant to the next x_cur
    y_cur = P_real.C*x_cur + P_real.D*u0;
    x_cur = P_real.A*x_cur + P_real.B*u0;
    
    % make next u_ini and y_ini with noise
    u_ini = [u_ini(P_real.nu+1:end); u0];
    y_ini = [y_ini(P_real.ny+1:end); P_real.gen_noise(y_cur,P_real.An)];
end
% temporary plot of simulation result
%v_drawSIM(Xdeepc,Udeepc,Simulation.nsim,P_real.Ts)

% old: MAE = norm( P_real.C*(Xdeepc - Trial.Xmpc), 1)/Simulation.nsim;
% new:
MAEsum = 0;
Y_err = P_real.C*(Xdeepc - Trial.Xmpc);
for i=1:Simulation.nsim
    MAEsum = MAEsum + norm( Y_err(:,i), 2);
end
MAE = MAEsum/Simulation.nsim;

err_code = 0;
end

function H = f_MakeHankelMat(X,p)
% Make a Hankel matrix of p row elements from X
% p: integer
% X: m by T
% H: m*p by T-p+1

T = size(X,2);

if T<p
    error('Size error to make Hankel matrix.');
end

H = [];
for i = 1:p
    H = [H; X(:,i:T-p+i)];
end

end

function UY = f_GetUY(plant,T,Tini,N,f_gen_sampling_input)
    
nu = plant.nu;
ny = plant.ny;

% Generate sampling input (should be PE); (nu by Tsampling)
U = f_gen_sampling_input(nu,T);

% Get the output of the plant
x = plant.x0;
Y = [];
for i = 1:T
    Y = [Y, plant.C * x + plant.D * U(:,i)];
    x = plant.A * x + plant.B * U(:,i);
end
Y = plant.gen_noise(Y, plant.An);   % Y is now corruped with noise

U_ = f_MakeHankelMat(U,Tini+N);
Y_ = f_MakeHankelMat(Y,Tini+N);
UY = [U_; Y_];
    
end

function [y0, x_cur] = f_FreeFlight(plant,T)
% free-flight of the plant to obtain initial condition for D2PC/DeePC
% under zero inputs of length T
% y0: the last T outputs (ny*T by 1)
% x_cur: the state at time 0 (nx by 1)
% the output is measured under noise

y0 = [];
x0 = plant.x0;
for i = 1:T
    y0 = [y0; plant.gen_noise(plant.C*x0,plant.An)];
    x_cur = x0;   % for saving the current state when completed
    x0 = plant.A*x0;
end

end
    
function plant = f_DefinePlant(name)
%
% plant.A, plant.B, plant.C, plant.D: DT parameters
% plant.Ts: Sampling time
% plant.An: Noise Level
% plant.nx, plant.nu, plant.ny: sizes of x, u, and y
% plant.x0: Initial condition
% plant.gen_noise:
%
switch name
    case 'InvPen'
        %%% Inverted Pendulum Modeling
        m = 0.23;       % mass of the rod
        M = 0.57;       % mass of the cart
        L = 0.64;       % length of the rod
        g = 9.8;        % gravity
        bp = 0.0024;    % lateral friction coefficient for cart
        bc = 4.3;       % rotational friction coefficient for rod
        % With x1: \theta   x2: \dot theta   x3: \x   x4: \dot x
        % get the plant (A, B, C, D) in the discrete-time 
        A0 = [0 1 0 0; 
            (M+m)/M/(L/2)*g -(M+m)/M/m/(L/2)^2*bp 0 bc/M/(L/2); 
            0 0 0 1; 
            -m/M*g bp/M/(L/2) 0 -bc/M];
        B0 = [0; -1/M/(L/2); 0; 1/M];   
        plant.Ts = 0.1;   % Sampling time for CT to DT conversion of the plant
        [plant.A, plant.B] = c2d(A0,B0,plant.Ts);  % Converting CT to DT system
        plant.C = [0 0 1 0];  % cart position is the output
        plant.D = 0;
        [plant.nx, plant.nu] = size(plant.B);    % nx = size of x, nu = size of u
        plant.ny = size(plant.C,1);         % ny = size of y
        plant.x0 = zeros(size(plant.B,1),1);  % initial condition of the plant
        plant.gen_noise = @(X,An) X + An * 2 * (rand(size(X))-0.5);
        %plant.An = 1e-4;   % Measurement Noise Level
        plant.An = 0;   % Measurement Noise Level

    case 'TwoCart'
        %%% Two Cart Modeling
        k = 2; b = 0; m1 = 1; m2 = 0.1;
        A0 = [0 1 0 0; -k/m1 -b/m1 k/m1 0; 0 0 0 1; k/m2 0 -k/m2 -b/m2];
        B0 = [0; 1/m1; 0; 0];
        plant.Ts = 0.1;   % Sampling time for CT to DT conversion of the plant
        [plant.A, plant.B] = c2d(A0,B0,plant.Ts);  % Converting CT to DT system
        plant.C = [0 0 1 0];  % cart position is the output
        plant.D = 0;
        [plant.nx, plant.nu] = size(plant.B);    % nx = size of x, nu = size of u
        plant.ny = size(plant.C,1);         % ny = size of y
        plant.x0 = zeros(size(plant.B,1),1);  % initial condition of the plant
        plant.gen_noise = @(X,An) X + An * 2 * (rand(size(X))-0.5);
        %plant.An = 1e-1;   % Measurement Noise Level
        plant.An = 1e-2;   % Measurement Noise Level        
        
    case 'FourTank'
        %%% Four Tank system 
        plant.A = [0.921 0 0.041 0; 0 0.918 0 0.033; 0 0 0.924 0; 0 0 0 0.937];
        plant.B = [0.017 0.001; 0.001 0.023; 0 0.061; 0.072 0];
        plant.C = [1 0 0 0; 0 1 0 0];
        plant.D = zeros(2,2);
        [plant.nx, plant.nu] = size(plant.B);    % nx = size of x, nu = size of u
        plant.ny = size(plant.C,1);         % ny = size of y
        plant.x0 = zeros(size(plant.B,1),1);  % initial condition of the plant
        plant.gen_noise = @(X,An) X + An * 2 * (rand(size(X))-0.5);
        plant.Ts = 1;
        plant.An = 1e-1;   % Measurement Noise Level
end

end


