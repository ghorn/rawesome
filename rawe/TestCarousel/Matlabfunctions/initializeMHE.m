function MHE = initializeMHE(MHE,Sim)

%% Read the measurements and build up a reference and initial guess

MHE = warm_start(MHE);

% Define and save the weighting matrices
MHE = MHEweights(MHE,Sim);
 
% generate the measurements
MHE = MHEmeasurements(MHE,Sim);


% Generate an initial guess for the states and controls
% Initial guess (for the first OCP)
MHE.Tguess = MHE.Tref(1:MHE.Ncvp+1);
MHE.Xguess = MHE.Xref(MHE.t_index:MHE.Ncvp+MHE.t_index,:);
MHE.Uguess = MHE.Uref(MHE.t_index:MHE.Ncvp+MHE.t_index-1,:);

% Save some values in the MHE structure for plotting and postprocessing
MHE.time = [];
MHE.state = [];
MHE.control = [];
MHE.constraints = [];
MHE.ddXIMU = [];
MHE.parameter = [];
MHE.KKT = [];

% Store the measurements (for postprocessing)
MHE.measM = MHE.Measurements(:,2:end);
MHE.measT = MHE.Tref-(MHE.Nref-1)*MHE.Ts;


