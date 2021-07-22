%%% 
%%% Simulation for D2PC (Data-Driven Predictive Control)
%%%
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
%%% Specify your solver option at OP.solver.
%%% If you don't have OSQP or CVX, choose the solver option = 'handful',
%%% but the control problem is solved without constraints.

clear

% Select your model; 
% (You can define your own at the end of this file in f_DefinePlant)
model = 'TwoCart';
%model = 'FourTank';

switch model
        
    case 'TwoCart'
        P_real = f_DefinePlant('TwoCart');

        OP.r = 1;       % target set-point of the output (ny by 1)
        OP.Q = 200;     % size: (ny by ny)
        OP.R = 1;       % size: (nu by nu)
        OP.N = 20;      % prediction horizon (step) of cost function (from 0 to N-1)
        OP.umax = 2;    % size: (nu by 1) Input constraint
        OP.umin = -2;   % size: (nu by 1) Input constraint
        OP.ymax = Inf;  % size: (ny by 1) Output constraint
        OP.ymin = -Inf; % size: (ny by 1) Output constraint
        OP.solver = 'osqp'; % solver = 'osqp' or 'cvx' or 'handful'
        
        % for Simulation
        Simulation.nsim = 60;   % simulation step

        % Parameters for d2pc
        d2pc.NBar = 20;          % Estimated order of the plant
        d2pc.Nd = 1;            % # of averaging for Ad and Bd
        d2pc.Tsampling = 100;   % length of sampling input (step)
        d2pc.gen_sampling_input = @(nu,T) 2*rand(nu,T)-1;   % generating function of sampling input 

    case 'FourTank'
        P_real = f_DefinePlant('FourTank');

        OP.r = [0.65; 0.77];    % target set-point of the output (ny by 1)
        OP.Q = 3*eye(2);        % size: (ny by ny)
        OP.R = 0.01*eye(2);     % size: (nu by nu)
        OP.N = 30;      % prediction horizon (step) of cost function (from 0 to N-1)
        OP.umax = [Inf; Inf];   % size: (nu by 1) Input constraint
        OP.umin = -OP.umax;     % size: (nu by 1) Input constraint
        OP.ymax = [Inf; Inf];   % size: (ny by 1) Output constraint
        OP.ymin = -OP.ymax;     % size: (ny by 1) Output constraint
        OP.solver = 'osqp'; % solver = 'osqp' or 'cvx' or 'handful'
        
        % for Simulation
        Simulation.nsim = 60;   % simulation step

        % Parameters for d2pc
        d2pc.NBar = 30;          % Estimated order of the plant
        d2pc.Nd = 5;            % # of averaging for Ad and Bd
        d2pc.Tsampling = 400;   % length of sampling input (step)
        d2pc.gen_sampling_input = @(nu,T) 2*rand(nu,T)-1;   % generating function of sampling input 
        

end

OP.P_qp = [];   % for internal use
OP.q_qp = [];   % for internal use
OP.A_qp = [];   % for internal use
OP.l_qp = [];   % for internal use
OP.u_qp = [];   % for internal use
OP.prob = [];   % for use of osqp
OP.ApBp = [];   % for use of d2pc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate: D2PC 
%%% Most of the codes in this block is for displaying.
%%% Key functions are performed in SimulateD2PC function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng('default');     % reset random number generator


%%% Offline sampling of real plant (data-driven model id)
[CMat,DMat] = f_D2PC_sampling(P_real, d2pc);

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
    fprintf('=== D2PC: step %d out of %d === \n', i, Simulation.nsim);
    [u0, err_code] = f_D2PC_solve(OP,chi_cur);
    if err_code == 1, 
        error('! D2PC failure in the solver'); 
    end
    Ud2pc = [Ud2pc, u0];
    
    % get the output of the plant and advance the plant to the next x_cur
    y_cur = P_real.C*x_cur + P_real.D*u0;
    x_cur = P_real.A*x_cur + P_real.B*u0;
    
    % make next chi_cur with noise
    chi_cur = [chi_cur(P_real.ny+1:P_real.ny*d2pc.NBar); ...
        P_real.gen_noise(y_cur,P_real.An); ...
        chi_cur(P_real.ny*d2pc.NBar+P_real.nu+1:end); u0];
end
% plot of simulation result
v_drawSIM(Xd2pc,Ud2pc,Simulation.nsim,P_real.Ts)

plot_time = (0:Simulation.nsim-1) * P_real.Ts;
% Now, plot_time, Xd2pc, and Ud2pc contain all the simulation information.

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
        plant.Ts = 1;   % Sampling time
        plant.An = 1e-1;   % Measurement Noise Level
end

end


