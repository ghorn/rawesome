function MHE = newMeasurement(MHE,MPC,Sim)

    
    %%%% Shift measurements %%%%
    
    % Take the last n-1 measurements (and shift them one timestep back)
    MHE.Measurements = [MHE.Tref(1:end-1,:) MHE.Measurements(2:end,2:end)]; %MHE.ddXref(MHE.t_index:MHE.Ncvp+MHE.t_index,:) MHE.Xref(MHE.t_index:MHE.Ncvp+MHE.t_index,MHE.state_index) MHE.Uref(MHE.t_index:MHE.Ncvp+MHE.t_index,:)];
    
    stateLast = Sim.state(end,:);
    
    % Compute the accelerations
    P.W0 = Sim.W0;
    P.tu = [0      MPC.U(1,:);
            MHE.Ts MPC.U(1,:);];
    P.r = Sim.r;
    Out = Model_integ_ACADO(0,stateLast,'IMU',P);
    ddXref = Out.';
    
    % Compute the markers
    P.delta = stateLast(19);
    P.ddelta = stateLast(20);
    P.dddelta = MHE.Uref(1,1);
    P.r = Sim.r;
    out = Model_integ_ACADO( 0, stateLast, 'markers', P );
    newMarkers = out(1:12).';
    
    % New measurement (from the simulation)
    newMeasurement = [newMarkers Sim.ddX(end,:) stateLast(:,MHE.state_index2) MPC.U(1,:)]; 

    % Add noise
    isNoise = Sim.Noise.is;
    factorNoise = Sim.Noise.factor;
    
    Noise = randn(size(newMeasurement))*factorNoise*diag(MHE.Scaling);
    newMeasurement = newMeasurement + Noise*isNoise;
    
    
    MHE.Measurements(end,end-3:end) = MPC.U(1,:);
    
    MHE.Measurements = [MHE.Measurements;
                        MHE.Tref(end,:) newMeasurement]; 

                    
    % New initial guess for the states (last trajectory + the state predicted by the MPC)
    MHE.Xguess = [MHE.X(2:end,:);MPC.nextInitialState]; 
        
    % For the controls pick the last MPC control
    uNext = MPC.U(1,:);
    
    uOLD = MHE.U(1,:);
    MHE.Uguess = [MHE.U(2:end,:); uNext;];

    
    MHE.measM = [MHE.measM; MHE.Measurements(end,2:end)];
    MHE.measT = [MHE.measT; MHE.measT(end)+MHE.Ts];

    
