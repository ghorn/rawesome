the RAWESOME Airborne Wind Energy Simulation, Optimization and Modeling Environment

This is the top level folder. Items of interest:

kite.proto - google protobuf definitions that many programs use
mkprotos.hs - a script which compiles kite.proto for kitevis and rawe
aero/ - some aerodynamic analysis stuff
rawe/ - python libraries for dynamics and optimal control
rawekite/ - kite specific python utilities
examples/ - official examples which should all work
studies/ - mostly broken python code which uses rawe,
           will eventually be deleted or moved to examples/
plot-ho-matic/ - live 2d plotting
wtfviz/ - Wind Turbine Flight Visualizer: live 3d simulation/optimization visualization using protobufs and zeromq in haskell
