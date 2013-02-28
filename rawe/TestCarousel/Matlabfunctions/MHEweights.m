function MHE = MHEweights(MHE,Sim)

% Read the covariances from the script Covariance.m:
% MHE.covariance_accelerations
% MHE.covariance_position
% MHE.covariance_rotation
% MHE.covariance_camera
% MHE.covariance_rotationalVelocity
% MHE.covariance_r
% MHE.covariance_dr
% MHE.covariance_delta
% MHE.covariance_ddelta
% MHE.covariance_controlSurfaces
% MHE.covariance_dddelta
% MHE.covariance_controls
Covariance

% Store the covariance in the matrices needed by ACADO: they come from Covariance.m
MHE.Q = diag([1/MHE.covariance_dddelta^2, 1/MHE.covariance_ddr^2 ones(1,2)/MHE.covariance_controls^2]);
MHE.R = diag([ones(1,12)/MHE.covariance_camera^2, ones(1,3)/MHE.covariance_rotationalVelocity^2, ones(1,3)/MHE.covariance_accelerations^2, 1/MHE.covariance_r^2, 1/MHE.covariance_delta^2, ones(1,2)/MHE.covariance_controlSurfaces^2]);
    

factor = 1;
% Put together the weights for the LSQ cost function in the MHE
% MHE.QQ = diag([diag(MHE.R);diag(MHE.Q);diag(MHE.QW)])*factor;
MHE.QQ = diag([diag(MHE.R);diag(MHE.Q)])*factor;

% MHE.QN = diag([diag(MHE.Q);diag(MHE.QW)])*factor;
MHE.QN = diag([diag(MHE.Q)])*factor;
MHE.RN = MHE.R*factor;

%% Save the weights

MHE.Qend = diag([ones(1,3)/MHE.covariance_rotationalVelocity^2, ones(1,3)/MHE.covariance_accelerations^2, 1/MHE.covariance_r^2, 1/MHE.covariance_delta^2, ones(1,2)/MHE.covariance_controlSurfaces^2]);


end

