import model
import joy
import casadi as C

#tc0 = 2*389.970797939731

x0 = C.DMatrix( [ 1.154244772411
                , -0.103540608242
                , -0.347959211327
                , 0.124930983341
                , 0.991534857363
                , 0.035367725910
                , 0.316039689643
                , -0.073559821379
                , 0.945889986864
                , 0.940484536806
                , -0.106993361072
                , -0.322554269411
                , 0.000000000000
                , 0.000000000000
                , 0.000000000000
                , 0.137035790811
                , 3.664945343102
                , -1.249768772258
                , 0.000000000000
                , 3.874600000000
                ])

r = 1.2

ts = 0.02

if __name__=='__main__':
    print "creating model"
    (dae, others) = model.model()
    dae.init()
    
    print "creating integrator"
    f = C.IdasIntegrator(dae)
    f.setOption("reltol",1e-6)
    f.setOption("abstol",1e-8)
    f.setOption("t0",0)
    f.setOption("tf",ts)
    f.init()
    
    js = joy.Joy()
    def simFun(x):
        axes = js.getAxes()
        u0 = -axes[0]*0.02
        u1 =  axes[1]*0.05
        tc = 600*(1 - axes[6])

        u = C.DMatrix([tc,u0,u1])
        f.setInput(x,C.INTEGRATOR_X0)
        f.setInput(u,C.INTEGRATOR_P)
        f.evaluate()
        return C.DMatrix(f.output())

    x = x0
    for k in range(0,100):
        x = simFun(x)
        print (k,x)
