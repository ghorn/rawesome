function [MPC,Sim] = getMPC(MPC,Sim)

    MPC.KKT = [MPC.KKT;MPC.stepKKT];

    %MPC.U = MPC.U.';
    %MPC.X = MPC.X.';
    
    % Define the "time"
    t = 0:MPC.Ts:MPC.Tc;
    
    % Save the controls
    MPC.controls = [MPC.controls;
                    MPC.time(end) MPC.U(1,:)];
    
    
    %%%%%%%%%%%%%%%  SIMULATE THE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P.W0 = Sim.W0 + Sim.DeltaW;
    P.tu = [0      MPC.U(1,:);
            MPC.Ts MPC.U(1,:);];
    P.r = Sim.r;
    
    [timeInter,stateInter] = ode45('Model_integ_ACADO',[0 MPC.Ts],Sim.state(end,:),Sim.intoptions,P);
    timeInter = timeInter + Sim.time(end);

    % Get the outputs for all times
    OutInter = [];
    ddX = [];
    const = [];
    for k = 1:length(timeInter)
        OutInter(k,:) = Model_integ_ACADO(timeInter(k),stateInter(k,:),'output',P);
        Out = Model_integ_ACADO(timeInter(k),stateInter(k,:),'IMU',P);     
        ddX(k,:)  = Out.';
        const(k,:) = Model_integ_ACADO(timeInter(k),stateInter(k,:),'const',P);     
    end
    
    Sim.controls = [Sim.controls;
                    Sim.time(end) MPC.U(1,:)];
                    
    Sim.time  = [Sim.time;timeInter];
    Sim.state = [Sim.state;stateInter];
    Sim.ddX   = [Sim.ddX; ddX;];
    Sim.Out   = [Sim.Out;OutInter];
    Sim.constraints   = [Sim.constraints;const];
   
%%%%%%%% GENERATE NEW INITIAL GUESS %%%%%%%%%%%%%
    
    %%%% Shift reference trajectories and inputs %%%%
    
    IG = stateInter(end,[3,21,22]);
    z = IG(1); delta = IG(2); 
    ddelta = IG(3); ddelta = 2*pi*MPC.Ref.RPM/60.;
    
    MPC.Ref.delta0 = delta;
    
    % Set the reference height and tether length
    timeINDEX = find( MPC.Ref.time > MPC.time(end), 1, 'first' ) - 1; % find the index pointing to the right reference
    MPC.Ref.z = MPC.Ref.vars(timeINDEX,1);
    MPC.Ref.r = MPC.Ref.vars(timeINDEX,2);
    MPC.Ref.ddelta = MPC.Ref.vars(timeINDEX,3);
    
    % Compute the new reference
    MPC = generate_reference(MPC);
    
    % Get the states
    stateProcess = MPC.X;
    
    OutInter = [];
    ddX = [];
    const = [];
    for k = 1:length(t)
        OutInter(k,:) = Model_integ_ACADO(t(k),stateProcess(k,:),'output',P);
        Out = Model_integ_ACADO(t(k),stateProcess(k,:),'IMU',P);     
        ddX(k,:)  = Out.';
        const(k,:) = Model_integ_ACADO(t(k),stateProcess(k,:),'const',P);     
    end
    
    MPC.time  = [MPC.time; MPC.time(end) + MPC.Ts];
    MPC.state = [MPC.state;stateProcess(2,:)];
    MPC.ddX   = [MPC.ddX; ddX(2,:);];
    MPC.Out   = [MPC.Out;OutInter(2,:)];
    MPC.constraints   = [MPC.constraints;const(end,:)];
    MPC.constraints0   = [MPC.constraints0;const(1,:)];
    
    
    % New guess: integrate the last interval starting from the last state
    % use the last control input
    uNext = MPC.U(end,:);
    
    P.tu = [0      uNext;
            MPC.Ts uNext];
    
    P.W0 = Sim.W0;
    [timeNext,stateNext] = ode45('Model_integ_ACADO',[0 MPC.Ts],stateProcess(end,:),Sim.intoptions,P);
    
    uNext = MPC.Uref(MPC.Ncvp+1,:);
    
    % New initial guess (last trajectory + the newly integrated part)
    stateProcess = [stateProcess(2:end,:); stateNext(end,:)];

    MPC.nextInitialState = stateProcess(1,:);
    
    MPC.Xguess = stateProcess;
    
    uOLD = MPC.U(1,:);
    MPC.Uguess = [MPC.U(2:end,:);uNext];
    
    
    % COMPUTE NEW TERMINAL COST
    % Linearize the system
    MPC.stateNext = MPC.Xref(MPC.Ncvp+1,:);
    
    stateLQR = MPC.stateNext;
    
    % Linearize and compute the LQR
    P.delta = stateLQR(19);
    P.ddelta = stateLQR(20);
    P.dddelta = uNext(1); 
    [MPC.A,MPC.B,~,~] = jacobianI(stateLQR,P);
    
%     MPC.A = MPC.A([1:18,21:22],[1:18,21:22]);
%     MPC.B = MPC.B([1:18,21:22],2:3);
    
    [K,MPC.S,e] = lqr(MPC.A,MPC.B,MPC.Q,MPC.R);
%     MPC.S = [ S(1:18,1:18), zeros(18,2),  S(1:18,19:20);
%                zeros(2,18), 1e-4*eye(2),     zeros(2,2);
%              S(19:20,1:18),  zeros(2,2), S(19:20,19:20);];
    
    MPC.S = MPC.Sfac*MPC.S-MPC.Q;
    MPC.K = K;
    MPC.storeS = [MPC.storeS;MPC.S;];
    

    % Compute the predicted lift coefficient
    k = 1;
    P.tu = [t(1) uOLD;
            t(2) uOLD;];
    P.W0 = Sim.W0;
    OutProcess = Model_integ_ACADO(t(k+1),stateProcess(k,:),'output',P);     
    CL = OutProcess(7);
    
    MPC.CLpred = [MPC.CLpred;[MPC.time(end), CL]];
    
    % Store the reference for postproc
    MPC.refX = [MPC.refX; MPC.Xref(2,:)];
    MPC.refU = [MPC.refU; MPC.Uref(2,:)];
    MPC.refT = [MPC.refT; MPC.time(end)+MPC.Ts];
    

    
