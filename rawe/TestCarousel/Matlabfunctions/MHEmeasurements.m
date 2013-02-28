function MHE = MHEmeasurements(MHE,Sim)

%%% Generate measurements
MHE.t_index = 1;

% states have a different ordering in the measurements
MHE.state_index2 = [19 21 23:24];


P.W0 = 0;           % Measurements always built with no wind
P.r = Sim.r;

% Create the markers
Markers = [];
for k = MHE.t_index:MHE.Nref-1+MHE.t_index
    
    P.tu = [0 MHE.Uref(k,:)];
    P.delta = MHE.Xref(k,21);
    P.ddelta = MHE.Xref(k,22);
    P.dddelta = MHE.Uref(k,1);
    P.r = Sim.r;
    
    out = Model_integ_ACADO( 0, MHE.Xref(k,:), 'markers', P );
    Markers = [Markers;out(1:12).';];
end

% Generate IMU measurements
P.W0 = 0;
P.r = Sim.r;
MHE.ddXref = [];
for k = 1:length(MHE.Tref)
    P.tu = [0      MHE.Uref(k,:);
            MHE.Ts MHE.Uref(k,:);];
    Out = Model_integ_ACADO(0,MHE.Xref(k,:),'IMU',P);
    MHE.ddXref = [MHE.ddXref; Out.'];
end

measurements = [Markers MHE.ddXref(MHE.t_index:MHE.Nref-1+MHE.t_index,:) MHE.Xref(MHE.t_index:MHE.Nref-1+MHE.t_index,MHE.state_index2) MHE.Uref(MHE.t_index:MHE.Nref-1+MHE.t_index,1:4)];

% Add the noise
isNoise = Sim.Noise.is;
factorNoise = Sim.Noise.factor;

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

% Scale the noise according to the measurement covariance
MHE.Scaling = [ones(1,12)*MHE.covariance_camera ones(1,3)*MHE.covariance_rotationalVelocity ones(1,3)*MHE.covariance_accelerations MHE.covariance_r MHE.covariance_delta ones(1,2)*MHE.covariance_controlSurfaces MHE.covariance_dddelta MHE.covariance_ddr ones(1,2)*MHE.covariance_controls];

Noise = randn(size(measurements))*factorNoise*diag(MHE.Scaling);
measurements = measurements + Noise*isNoise;

MHE.Measurements = [MHE.Tref(1:MHE.Nref,:) measurements];


end

