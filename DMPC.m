clear all
close all
clc
import casadi.*

addpath(genpath('D:\adosreis\Documents\MATLAB'))
addpath(genpath('D:\adosreis\Documents\Optimization'))
addpath(genpath('D:\adosreis\Desktop\researchcodes\MPC & RG\'))

%% Define system, constraint sets, MPC
opt.N           = 30;
opt.n_controls  = 2;
opt.n_states    = 3;
opt.dt = 0.1;
opt.continuous_model.integration = 'RK4';
opt.model.function	= @(x,u) [u(1)*cos(x(3));
                              u(1)*sin(x(3));
                              u(2)];
N_robots = 3;

% define parameters of the optimization and their dimensions
maxrows = 10; % limit the lines of the description Ax<=b for the Voronoi cell
opt.parameters.name = {'pos_ref','local_targ','pos_obs','pos_neighbors1','pos_neighbors2','VoroA','Vorob'};
opt.parameters.dim = [opt.n_states 1; opt.n_states-1 1; opt.n_states-1 1; opt.n_states-1 1; opt.n_states-1 1; maxrows 2; maxrows 1];

% define cost functions
opt.costs.stage.parameters = {'pos_ref','local_targ'};

Q1 = 1000*ones(2);
Q2 = blkdiag(100000,10000);
Q3 = 1000;
R = 0.001*eye(2);
rdest = 0;
dobs_safe = 0.1;
dcol_safe = 0.1;
lambda_obs = 100000000;
lambda_col = 0.1;
Pobs = 0.5*eye(2);
Pcol = 1*eye(2);
pos_obs = [-1.4 4 5; 1.8 -2 6]; % obstacles with fixed position

% this cost takes every x_k into account
opt.costs.stage.function = @(x,u,varargin) (x([1:2])-varargin{2})'*Q1*(x([1:2])-varargin{2}) + ...
                                           (x(3)-varargin{1}(3))'*Q3*(x(3)-varargin{1}(3)) + u'*R*u;


% this cost does not take x_k into account, only the other variables
opt.costs.general.parameters = {'pos_ref','local_targ','pos_neighbors1','pos_neighbors2'}; 

% loop for collision with several agents
func_col  = @(varargin) 0;
for i = 1:(N_robots-1)
    func_col  = @(varargin) func_col(varargin{:}) + lambda_col/(1+exp((varargin{2}-varargin{2+i})'*Pcol*(varargin{2}-varargin{2+i})-dcol_safe));
end

% loop for collision with several obstacles
func_obs  = @(x,varargin) 0;
for i = 1:length(pos_obs)
    func_obs  = @(x,varargin) func_obs(x,varargin{:}) + lambda_obs/(1+exp((varargin{2}-pos_obs(:,i))'*Pcol*(varargin{2}-pos_obs(:,i))-dcol_safe));
end

opt.costs.general.function = @(x,varargin) (varargin{2}-varargin{1}(1:2))'*Q2*(varargin{2}-varargin{1}(1:2)) - rdest^2 + ...
                                           func_obs(x,varargin{:}) + ... 
                                           func_col(varargin{:});

% define hard constraints (states and inputs)
opt.constraints.states.upper  =   [50;50;2*pi];
opt.constraints.states.lower  =  -[50;50;2*pi];
opt.constraints.control.upper =   [1;pi/6];
opt.constraints.control.lower =  -[0;pi/6];

% define general constraints parameters
opt.constraints.general.parameters  = {'local_targ','VoroA','Vorob'};

% one for the end point constraint
opt.constraints.general.function{1} = @(x,varargin) (x([1:2],end-1)-varargin{1}(1:2));
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1}     = 'equality';

% another one for the Voronoi appartenance
opt.constraints.general.function{2} = @(x,varargin) ([varargin{2}(1:maxrows,1) varargin{2}(maxrows+1:maxrows*2,1)]*varargin{1} - varargin{3});
opt.constraints.general.elements{2} = 'N';
opt.constraints.general.type{2}     = 'inequality';

% parameters that are constrained
opt.constraints.parameters = {'local_targ'};

% define inputs (vectors and matrices)
opt.input.vector = {'pos_ref','pos_neighbors1','pos_neighbors2','Vorob'}; %
opt.input.matrix = {'VoroA'};

% get solver
opt.solver      = 'ipopt';
[solver,args]   = build_mpc(opt);

%% Simulation
tmax = 500;


% initial conditions
robot(:,1,1)    = [-9; -2; pi/4];
robot(:,1,2)    = [-5; -5; pi/4];
robot(:,1,3)    = [0; -5; pi/4];

for j = 1:N_robots
    x0(:,j) = zeros(length(args.vars{1}),1);
end

voroA = zeros(maxrows*2,1);
vorob = zeros(maxrows,1);

%  simulation loop
for t = 1:tmax
    % reference position/orientation
    ref_pos(:,t,1) = [10;10;pi/4];
    ref_pos(:,t,2) = [10;10;pi/4];
    ref_pos(:,t,3) = [10;10;pi/4];
    
    % obstacle at origin

    for j = 1:N_robots
        
        % inputs to OCP
        others = find([1:N_robots]~=j);
        args.x0 = x0(:,j);
        [voro] = voronoicell(robot(:,t,j),robot(:,t,others));
        voroA = reshape(vertcat(voro.A,zeros(maxrows-length(voro.A),2)),maxrows*2,1);
        vorob = vertcat(voro.b,zeros(maxrows-length(voro.b),1));
        args.p          = [robot(:,t,j); ref_pos(:,t,j); vertcat(robot(1:2,t,others(1)),robot(1:2,t,others(2)));vorob;voroA]; %
        [usol,xpred,rest,ct,fullsol] = solve_mpc(opt,solver,args);
        u(:,t,j) = usol(:,1);
        cost(j,t) = ct;
        ref(:,t,j) = rest;
        robot(:,t+1,j) = robot(:,t,j) + opt.dt*opt.model.function(robot(:,t,j),u(:,t));
        x0(:,j) = fullsol;
        t
    end
end

%% Get plots
close all
figure(1)
xs = linspace(-10,20,50);
ys = linspace(-10,15,50);
hold on
func = @(x,y) 0;
for i = 1:length(pos_obs)
    func = @(x,y) func(x,y) + lambda_obs/(1+exp(([x;y]-pos_obs(:,i))'*Pobs*([x;y]-pos_obs(:,i))-dobs_safe));
end
for i = 1:50
    for j = 1:50
        Z(i,j) = func(xs(i),ys(j));
    end
end

% plot negative Z so graph can be read
mesh(xs,ys,-Z')

for i = 1:length(pos_obs)
plot(pos_obs(1,i),pos_obs(2,i),'xk','LineWidth',2)
end

plot(robot(1,:,1),robot(2,:,1),'-ok')
plot(robot(1,:,2),robot(2,:,2),'-or')
plot(robot(1,:,3),robot(2,:,3),'-ob')
plot(ref(1,:,1),ref(2,:,1),'-xk')
plot(ref(1,:,2),ref(2,:,2),'-xr')
plot(ref(1,:,3),ref(2,:,3),'-xb')

% plot distances
for t = 1:tmax
dist(1,t) = norm(robot(:,t,1)-robot(:,t,2));
dist(2,t) = norm(robot(:,t,1)-robot(:,t,3));
dist(3,t) = norm(robot(:,t,2)-robot(:,t,3));
end
figure
plot(dist(1,:),'r');
hold on
plot(dist(2,:),'g');
plot(dist(3,:),'b');