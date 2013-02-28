function [MPC,Sim] = initializeMPC(MPC,Sim)


%% Build up a reference and an initial guess (starting from an equilibrium)

% Read the first reference
MPC.Ref.z = MPC.Ref.vars(1,1);
MPC.Ref.r = MPC.Ref.vars(1,2);
% MPC.Ref.delta = 0;
MPC.Ref.ddelta = MPC.Ref.vars(1,3);

MPC.Ref.delta0 = 0;

% Compute the equilibrium and generate the reference
MPC = generate_reference(MPC);

% Initial guess (for the first OCP)
MPC.Tguess = MPC.Tref(1:MPC.Ncvp+1);
MPC.Xguess = MPC.Xref(1:MPC.Ncvp+1,:);
MPC.Uguess = MPC.Uref(1:MPC.Ncvp,:);


%% BUILD REFERENCE
 
% COMPUTE NEW TERMINAL COST
% Linearize the system
MPC.stateNext = MPC.Xref(end,:);
MPC.uNext = MPC.Uref(end,:);

P.W0 = 0;           % Reference always built with no wind
P.tu = [0 MPC.uNext];
P.delta = MPC.stateNext(1,19);
P.ddelta = MPC.stateNext(1,20);
P.dddelta = MPC.uNext(1); 
P.r = MPC.Ref.r;

[MPC.A,MPC.B,~,~] = jacobianI(MPC.stateNext,P);  % Linearize the system at the equilibrium

% % Need to eliminate the part of the matrices relative to the dynamics of the rotation of the carousel
% MPC.A = MPC.A([1:18,21:22],[1:18,21:22]);
% MPC.B = MPC.B([1:18,21:22],2:3);

% Define the weights for the LQR (they will be used also by the NMPC)
We = 10;
Ww = 0.1;
% Control weights
MPC.R = 1e-4*diag([100, 100, 100, 100]);
% State weights
MPC.Q = 1e-2*diag([25,  25,  25, 25, 25, 25,    ...
                We,  We,  We,  We,  We,  We,  We,  We,  We, ...
                Ww,  Ww,  Ww, 100, 100,  100, 100, 250, 250]);

% % Weights without the carousel rotational dynamics
% MPC.Q1 = MPC.Q([1:18,21:22],[1:18,21:22]);
% MPC.R1 = MPC.R(2:3,2:3);

% Compute the LQR
[K,MPC.S,e] = lqr(MPC.A,MPC.B,MPC.Q,MPC.R);

% % Add to the terminal cost matrix the missing lines / columns
% MPC.S = [ S(1:18,1:18), zeros(18,2),  S(1:18,19:20);
%            zeros(2,18), 1e-4*eye(2),     zeros(2,2);
%          S(19:20,1:18),  zeros(2,2), S(19:20,19:20);];

% Multiply the terminal cost by the desired factor
MPC.S = MPC.Sfac*MPC.S-MPC.Q;
MPC.K = K;
MPC.storeS = MPC.S; % Store all the terminal costs (just for debugging)


%% Save some useful variables in the MPC structure, to use them for plotting and postprocessing
MPC.X0 = MPC.Xguess(1,:);
P.W0 = Sim.W0;
OutPred = Model_integ_ACADO(0,MPC.X0,'output',P);
Out = Model_integ_ACADO(0,MPC.X0,'IMU',P);
ddX = Out.';
const = Model_integ_ACADO(0,MPC.X0,'const',P);

MPC.time = 0;
MPC.state = MPC.X0;
MPC.ddX = ddX;
MPC.controls = [];
MPC.Out = OutPred.';
MPC.constraints = const;
MPC.constraints0 = const;
MPC.CLpred = [0, OutPred(7)];
MPC.KKT = [];

% Store the reference points (for postprocessing)
MPC.refX = MPC.Xref(1:2,:);
MPC.refU = MPC.Uref(1:2,:);
MPC.refT = MPC.Tref(1:2,:);

MPC.tu = [];
MPC.tSlack = [];

% Store the same things for the "real" plant
Sim.time  = 0;
Sim.state = MPC.X0;
Sim.ddX = ddX;
Sim.controls = [];
Sim.Out = OutPred.';
Sim.constraints = const;

