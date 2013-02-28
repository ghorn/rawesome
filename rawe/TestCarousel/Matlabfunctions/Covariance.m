MHE.covariance_accelerations = 10;
MHE.covariance_position = 1;
MHE.covariance_rotation = 1;
MHE.covariance_camera = 1e3;
MHE.covariance_rotationalVelocity = 1;
MHE.covariance_r = 1e-3;                % Covariance of the tether length
MHE.covariance_dr = 1e-3;               % Covariance of the tether velocity
MHE.covariance_delta = 1e-3;            % Covariance of the carousel angle
MHE.covariance_ddelta = 1e-3;           % Covariance of the carousel rotational velocity
MHE.covariance_controlSurfaces = 1e-6;
MHE.covariance_dddelta = 1e-4;          % Covariance of the carousel rotational acceleration
MHE.covariance_ddr = 1e-4;              % Covariance of the tether acceleration
MHE.covariance_controls = 1e-6;            % Covariance of the carousel control surfaces time derivative

