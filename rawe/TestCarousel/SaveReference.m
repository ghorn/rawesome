clear all
% close all
clc

addpath('Simulation','Matlabfunctions')

%% Most of the initialization is just useless, the only important part is
% the initial condition (just after this comment)

% Initial condition
MPC.Ref.z = -0.1189362777884522 ;         % reference trajectory height (relative to the arm)
MPC.Ref.r = 1.2;            % tether length
MPC.Ref.delta = 0;          % initial carousel angle
MPC.Ref.RPM = 60;           % carousel rotational velocity
MPC.Ref.ddelta = 2*pi*MPC.Ref.RPM/60.;
Sim.r = MPC.Ref.r;          % copy the tether length (just for the ease of use)

%% Dummy part, just to let the functions not complain (they do more than what we need...)
% The only parameter that could be useful to change is the sampling time

% MPC settings
MPC.Tc = 0.5; % horizon in seconds
MPC.Ncvp = 5; % number of cvp
MPC.Ts  = MPC.Tc/MPC.Ncvp; % sampling time
MPC.Nref = MPC.Ncvp+1; % 500;% % number of elements in the reference file
% Multiplying factor for the terminal cost
MPC.Sfac = 1; % factor multiplying the terminal cost


% Set integration parameters
Sim.intoptions.AbsTol = 1e-8;
Sim.intoptions.RelTol = 1e-8;
Sim.intoptions.MaxStep = 1/10;
Sim.intoptions.Events = 'off';

% Wind parameters
Sim.W0 = 0;         % Wind known by the controller
Sim.DeltaW = 0;     % Windgusts or unknown wind

%% Collect the parameters defining the reference in one array and comput the stuff
MPC.Ref.vars(1,1:3)   = [MPC.Ref.z MPC.Ref.r MPC.Ref.ddelta];

% Compute the reference and the terminal cost
[MPC,Sim] = initializeMPC(MPC,Sim);

% Reference trajectory
Xref = MPC.Xref;
Uref = MPC.Uref;
Tref = MPC.Tref;

% Terminal cost matrix
S = MPC.S;