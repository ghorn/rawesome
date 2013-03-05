the RAWESOME Airborne Wind Energy Simulation, Optimization and Modeling Environment

This is the top level folder. Items of interest:

kite.proto - google protobuf definitions that many programs use
mkprotos.hs - a script which compiles kite.proto for kitevis and rawe
aero/ - some aerodynamic analysis stuff
rawe/ - python libraries for dynamics and optimal control (and kite stuff)
studies/ - python code which uses rawe to do fun stuff!
wtfviz/ - Wind Turbine Flight Visualizer: live 3d simulation/optimization visualization using protobufs and zeromq in haskell
