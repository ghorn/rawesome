function MPC = getEstimatedState(MHE,MPC,Sim)


% Initial condition
if Sim.decoupleMPC_MHE
    MPC.X0 = Sim.state(end,:);
else
    MPC.X0 = MHE.state(end,:);
end
