import os
import subprocess_tee
import casadi as C

def getIndices(sym):
    (nr,nc) = sym.shape
    if (nr,nc) == (1,1):
        indices = ['']
    elif 1 in [nr,nc]:
        indices = ['__'+str(j) for j in range(nr*nc)]
    else:
        indices = []
        for jr in range(nr):
            for jc in range(nc):
                indices.append('__'+str(jr)+'_'+str(jc))
    return indices

def daeNames(dae,measurements,measurementsEnd):
    return [('DifferentialStates',dae.xNames()),
            ('AlgebraicVars',dae.zNames()),
            ('Controls',dae.uNames()),
            ('Parameters',dae.pNames()),
            ('Outputs', dae.outputNames()),
            ('Measurements',measurements),
            ('MeasurementsEnd',measurementsEnd)]

def simpleMessage(dae, messagename, fieldnames, qualifier='required'):
    ret = []
    ret.append('message '+messagename+' {')
    k = 1
    for name in fieldnames:
        for idx in getIndices(dae[name]):
            ret.append('  '+qualifier+' double '+name+idx+' = '+str(k)+';')
            k += 1
    ret.append('}\n\n')
    return '\n'.join(ret)

def writeProtoSpec(topname, dae, measurements, measurementsEnd):
    protobufs = 'package '+topname+';\n\n'
    for (msgname,fieldnames) in daeNames(dae, measurements, measurementsEnd):
        protobufs += simpleMessage(dae, msgname, fieldnames)
        #if msgname == 'Measurements':
        #    protobufs += simpleMessage(dae, msgname, fieldnames,qualifier='optional')
        #else:
        #    protobufs += simpleMessage(dae, msgname, fieldnames)

    protobufs += '''\
message Dae {
  required DifferentialStates differentialStates = 1;
  optional AlgebraicVars algebraicVars = 2;
  required Controls controls = 3;
  required Parameters parameters = 4;
  optional Outputs outputs = 5;
}

message Trajectory {
  repeated Dae traj = 1;
  required int32 iteration = 2;
  repeated string messages = 3;
}
'''
    return protobufs

def writePythonGenerator(topname,dae):
    lines = ['import '+topname+'_pb2']
    lines.append('')
    lines.append('def toProto(lookup):')
    lines.append('    dkp = '+topname+'_pb2.'+topname+'State()')
    lines.append('')
    for (field,names) in [('x',dae.xNames()),
                          ('u',dae.uNames()),
                          ('p',dae.pNames())]:
        for name in names:
            lines.append('    dkp.'+field+'.'+name+' = lookup(\''+name+'\')')
        lines.append('')
    for (field,names) in [('z',dae.zNames()),
                          ('outs',dae.outputNames())]:
        for name in names:
            lines.append('    try:')
            lines.append('        dkp.'+field+'.'+name+' = lookup(\''+name+'\')')
            lines.append('    except Exception:')
            lines.append('        pass')
    lines.append('')
    lines.append('    return dkp')
    return ('\n'.join(lines)).strip()+'\n'

def writeStruct(vecname,fieldnames,dae,realType='double'):
    ret = []
    ret.append('\ntypedef struct {')
    k = 0
    for name in fieldnames:
        for idx in getIndices(dae[name]):
            ret.append('  '+realType+' '+name+idx+'; /* '+str(k)+' */')
            k += 1
    ret.append('} '+vecname+';\n')
    return '\n'.join(ret)

def writeStructs(dae, topname, measurements, measurementsEnd):
    structs = []
    structs.append('#ifndef __'+topname+'_STRUCTS_H__')
    structs.append('#define __'+topname+'_STRUCTS_H__')
    for (structName, fieldNames) in daeNames(dae, measurements, measurementsEnd):
        structs.append(writeStruct(structName,fieldNames,dae))
    structs.append('#endif // __'+topname+'_STRUCTS_H__')
    return '\n'.join(structs)

def writeProtoConverter(topname, vecname, fieldnames, dae):
    ret0 = []
    ret1 = []
    prototype0 = 'void from'+vecname+'('+topname+'::'+vecname+' * proto, const '+vecname+' * data)'
    prototype1 = 'void to'+vecname+'('+vecname+' * data, const '+topname+'::'+vecname+' * proto)'
    ret0.append(prototype0+' {')
    ret1.append(prototype1+' {')
    for name in fieldnames:
        for idx in getIndices(dae[name]):
            ret0.append('  proto->set_'+name.lower()+idx+'(data->'+name+idx+');')
            ret1.append('  data->'+name+idx+' = proto->'+name.lower()+idx+'();')
    ret0.append('}\n\n')
    ret1.append('}\n\n')
    return ('\n'.join(ret0), prototype0+';\n','\n'.join(ret1), prototype1+';\n')

def writeProtoConverters(dae, topname, measurements, measurementsEnd):
    protoConverters = '#include "protoConverters.h"\n\n'
    protoConverterHeader = []
    def write(blah):
        protoConverterHeader.append(blah)
    write('''\
#ifndef __PROTO_CONVERTER_H__
#define __PROTO_CONVERTER_H__

#include "%(topname)s_structs.h"
#include "%(topname)s.pb.h"

#ifdef __cplusplus
extern "C" {
#endif

''' % {'topname':topname})
    for (vecname,fieldnames) in daeNames(dae, measurements, measurementsEnd):
        (pc0,pcpt0,pc1,pcpt1) = writeProtoConverter(topname, vecname, fieldnames, dae)
        protoConverters += pc0
        protoConverters += pc1
        write(pcpt0)
        write(pcpt1)
        write('\n')
    write('''
#ifdef __cplusplus
}
#endif

#endif // __PROTO_CONVERTER_H__
''')
    return (protoConverters,''.join(protoConverterHeader))


def writeDimensions(topname, dae, meas, measEnd, mheHorizN, mpcHorizN):
    nMeas    = C.veccat([dae[n] for n in meas]).size()
    nMeasEnd = C.veccat([dae[n] for n in measEnd]).size()
    ret = []
    ret.append('#ifndef __'+topname+'_DIMENSIONS_H__')
    ret.append('#define __'+topname+'_DIMENSIONS_H__')
    ret.append('\n')
    ret.append('#define NUM_DIFFSTATES       '+str(len(dae.xNames())))
    ret.append('#define NUM_ALGVARS          '+str(len(dae.zNames())))
    ret.append('#define NUM_CONTROLS         '+str(len(dae.uNames())))
    ret.append('#define NUM_PARAMETERS       '+str(len(dae.pNames())))
    ret.append('#define NUM_MEASUREMENTS     '+repr(nMeas))
    ret.append('#define NUM_MEASUREMENTS_END '+repr(nMeasEnd))
    ret.append('#define NUM_MHE_HORIZON      '+str(mheHorizN))
    ret.append('#define NUM_MPC_HORIZON      '+str(mpcHorizN))
    ret.append('\n')
    ret.append('#endif // __'+topname+'_DIMENSIONS_H__')
    return '\n'.join(ret)

def writeAll(dae, topname, autogenDir, measurements=[], measurementsEnd=[], mheHorizN=0, mpcHorizN=0):
    # make autogen directory if it doesn't exist
    if not os.path.exists(autogenDir):
        os.makedirs(autogenDir)

    # write files in autogen directory
    protobufs = writeProtoSpec(topname,dae,measurements,measurementsEnd)

    protobufs += '''
message Debug {
  optional double d0 = 1;
  optional double d1 = 2;
  optional double d2 = 3;
  optional double d3 = 4;
  optional double d4 = 5;
  optional double d5 = 6;
  optional double d6 = 7;
  optional double d7 = 8;
  optional double d8 = 9;
  optional double d9 = 10;
}

message Mhe {
  repeated DifferentialStates x = 1;
  repeated Controls u = 2;
  repeated Measurements y = 3;
  required MeasurementsEnd yN = 4;
  required double kkt = 5;
  required double objective = 6;
  required double prepTime = 7;
  required double fbTime = 8;
}

message Mpc {
  repeated DifferentialStates x = 1;
  repeated Controls u = 2;
  required DifferentialStates x0 = 3;
  //repeated Dae referenceTrajectory = 4;
  required double kkt = 5;
  required double objective = 6;
  required double prepTime = 7;
  required double fbTime = 8;
}

message MheMpcHorizons {
  required Mhe mhe = 1;
  required Mpc mpc = 2;
  required DifferentialStates mheXN = 3;
  required Controls mpcU0 = 4;
  optional DifferentialStates simX = 5;
  optional Controls simU = 6;
  repeated string messages = 7;
  required Debug debug = 8;
}
'''

    f = open(os.path.join(autogenDir,topname+'.proto'),'w')
    f.write(protobufs)
    f.close()

    # dimensions
    dims = writeDimensions(topname, dae, measurements, measurementsEnd, mheHorizN, mpcHorizN)
    f = open(os.path.join(autogenDir, topname+'_dimensions.h'),'w')
    f.write(dims)
    f.close()

    # structs
    structs = writeStructs(dae, topname, measurements, measurementsEnd)
    f = open(os.path.join(autogenDir, topname+'_structs.h'),'w')
    f.write(structs)
    f.close()
    
    # to python proto
    toProtobufs = writePythonGenerator(topname,dae)
    f = open(os.path.join(autogenDir, 'to'+topname+'Proto.py'),'w')
    f.write(toProtobufs)
    f.close()
    
    # __init__.py
    module = '''\
import to%(topname)sProto
import %(topname)s_pb2
''' % {'topname':topname}
    f = open(os.path.join(autogenDir, '__init__.py'),'w')
    f.write(module)
    f.close()

    # C++ proto converters
    (protoConverters, protoConverterHeader) = \
        writeProtoConverters(dae, topname, measurements, measurementsEnd)
    f = open(os.path.join(autogenDir, 'protoConverters.cpp'),'w')
    f.write(protoConverters)
    f.close()
    f = open(os.path.join(autogenDir, 'protoConverters.h'),'w')
    f.write(protoConverterHeader)
    f.close()



#    # call protoc to make python/C++ code
#    (ret, msgs) = subprocess_tee.call(
#        ['protoc','--cpp_out=.','--python_out=.',os.path.join(autogenDir,topname+'.proto')])
#    if ret != 0:
#        raise Exception('protoc fail\n'+msgs)
#    for haskellDir in haskellDirs:
#        (ret,msgs) = subprocess_tee.call(
#            ['hprotoc','-I../'+autogenDir,'--haskell_out=src',topname+'.proto'],cwd=haskellDir)
#        if ret != 0:
#            raise Exception('hrotoc fail\n'+msgs)

