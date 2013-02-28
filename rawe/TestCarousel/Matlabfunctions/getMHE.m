function MHE = getMHE(MHE,Sim)

    MHE.t_index = MHE.t_index + 1;
    

%     MHE.U = MHE.U.';
%     MHE.X = MHE.X.';
    
    MHE.KKT = [MHE.KKT; MHE.stepKKT];

    % Define the "time"
    t = 0:MHE.Ts:MHE.Tc;

    % Save the controls
    tu = [t.' [MHE.U;MHE.U(end,:)]];
    P.tu = tu;
    
    MHE.control = [MHE.control;MHE.U(end,:)];
    MHE.previousControl = MHE.U(end,:);
    
    % Save the state
    stateProcess = MHE.X;
    
    MHE.state = [MHE.state;stateProcess(end,:)];
    
    MHE.previousState = stateProcess;
    

    if length(MHE.time) == 0
        MHE.time = 0;
    else
        MHE.time   = [MHE.time; MHE.time(end) + MHE.Ts];
    end
    
    P.W0 = Sim.W0;
    P.r = Sim.r;
    ConstrProcess = Model_integ_ACADO(MHE.time(end),MHE.state(end,:),'const',P);
    MHE.constraints = [MHE.constraints;
                          ConstrProcess];
    
    ddXProcess = Model_integ_ACADO(MHE.time(end),MHE.state(end,:),'IMU',P);
    MHE.ddXIMU = [MHE.ddXIMU;
                     ddXProcess.';];
    
