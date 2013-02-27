import casadi as C

def unpack(d,name,msg):
    try:
        return d[name]
    except KeyError:
        raise Exception('you forgot to put "'+name+'" in '+msg)
            

class Sim(object):
    def __init__(self, dae, ts):
#        x = C.DMatrix([x0[name] for name in self.xNames()])
    
#        print "creating outputs function"
#        (_,(fNoZ,outputNames)) = self.outputsFun()
#        fNoZ.init()
#    
#        print "creating communicator"
#        communicator = simutils.Communicator(fNoZ,outputNames)

        print "creating integrator"
        self.dae = dae
        self.integrator = C.IdasIntegrator(self.dae.casadiDae())
        self.integrator.setOption("reltol",1e-6)
        self.integrator.setOption("abstol",1e-8)
        self.integrator.setOption("t0",0)
        self.integrator.setOption("tf",ts)
        self.integrator.setOption('name','integrator')
#        self.integrator.setOption("linear_solver_creator",C.CSparse)
#        self.integrator.setOption("linear_solver","user_defined")
        self.integrator.init()
        
        print "creating outputs function"
        (fAll, (f0,outputs0names)) = self.dae.outputsFun()
        self.outputsFunAll = fAll
        self.outputsFun0 = f0
        self.outputs0names = outputs0names

    def step(self, x, u, p):
        xVec = C.veccat([unpack(x,name,"the states") for name in self.dae.xNames()])
        uVec = C.veccat([unpack(u,name,"the controls") for name in self.dae.uNames()])
        pVec = C.veccat([unpack(p,name,"the parameters") for name in self.dae.pNames()])
        self.integrator.setInput(xVec,C.INTEGRATOR_X0)
        self.integrator.setInput(C.veccat([uVec,pVec]),C.INTEGRATOR_P)
        self.integrator.evaluate()
        xNext = C.DMatrix(self.integrator.output())
        ret = {}
        for k,name in enumerate(self.dae.xNames()):
            ret[name] = xNext[k].at(0)
        return ret

    def getOutputs(self, x, u, p):
        xVec = C.veccat([unpack(x,name,"the states") for name in self.dae.xNames()])
        uVec = C.veccat([unpack(u,name,"the controls") for name in self.dae.uNames()])
        pVec = C.veccat([unpack(p,name,"the parameters") for name in self.dae.pNames()])
        self.outputsFun0.setInput(xVec, 0)
        self.outputsFun0.setInput(uVec, 1)
        self.outputsFun0.setInput(pVec, 2)
        self.outputsFun0.evaluate()
        outs = C.DMatrix(self.integrator.output())
        ret = {}
        for k,name in enumerate(self.outputs0names):
            ret[name] = outs[k]
        return ret
