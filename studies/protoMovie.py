import time
import zmq
import pickle
import sys

if __name__=='__main__':
#    filename = 'data/crosswind_protos.dat'
#    filename = 'data/crosswind_protos_unconstrainedCl.dat'
#    filename = 'data/crosswind_protos_4_loops_unconstrainedCl.dat'
    filename = 'data/crosswind_protos_4_loops_unconstrainedCl_sparse.dat'
    f=open(filename,'r')
    protos = pickle.load(f)
    f.close()

    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    time.sleep(0.1)
    if len(sys.argv)>1 and sys.argv[1]=='-fst':
        publisher.send_multipart(["multi-carousel", protos[0]])
        sys.exit()

    k = 0
    dt = 1.0/20.0
    tNext = time.time() + dt
    for mcStr in protos:
        publisher.send_multipart(["multi-carousel", mcStr])
        print "wooo ",k
        k += 1
        toSleep = tNext - time.time()
        if toSleep > 0:
            time.sleep(toSleep)
        tNext += dt
