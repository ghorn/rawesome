import os
import subprocess

def _makeProtoBuf(protoname, fieldnames, qualifier='required'):
    ret = []
    ret.append('message '+protoname+' {')
    for k,name in enumerate(fieldnames):
        ret.append('  '+qualifier+' double '+name+' = '+str(k+1)+';')
    ret.append('}\n\n')
    return '\n'.join(ret)

def _toProtos(topname,dae):
    protobufs = ''
    for (protoname,fieldnames) in \
        [('DifferentialStates',dae.xNames()),
         ('Controls',dae.uNames()),
         ('Parameters',dae.pNames())]:
        protobufs += _makeProtoBuf(protoname, fieldnames)

    for (protoname,fieldnames) in \
         [('AlgebraicVars',dae.zNames()),
          ('Outputs', dae.outputNames())]:
        protobufs += _makeProtoBuf(protoname, fieldnames, qualifier='optional')

    protobufs += '''\
message %(topname)sState {
  required DifferentialStates x = 1;
  optional AlgebraicVars z = 2;
  required Controls u = 3;
  optional Parameters p = 4;
  optional Outputs outs = 5;
}

message %(topname)sTrajectory {
  repeated %(topname)sState traj = 1;
  required int32 iteration = 2;
  repeated string messages = 3;
}
''' % {'topname':topname}
    return protobufs

def _toToProto(topname,dae):
    lines = ['import '+topname+'_pb2']
    lines.append('')
    lines.append('def toProto(lookup):')
    lines.append('    dkp = '+topname+'_pb2.'+topname+'State()')
    lines.append('')
    for (field,names) in [('x',dae.xNames()),
                          ('u',dae.uNames())]:
        for name in names:
            lines.append('    dkp.'+field+'.'+name+' = lookup(\''+name+'\')')
        lines.append('')
    for (field,names) in [('z',dae.zNames()),
                          ('p',dae.pNames()),
                          ('outs',dae.outputNames())]:
        for name in names:
            lines.append('    try:')
            lines.append('        dkp.'+field+'.'+name+' = lookup(\''+name+'\')')
            lines.append('    except Exception:')
            lines.append('        pass')
    lines.append('')
    lines.append('    return dkp')
    return ('\n'.join(lines)).strip()+'\n'

def writeAll(dae, topname, autogenDir,haskellDirs=[],extraProtos=None,package=None):
    protobufs = _toProtos(topname,dae)
    if package is not None:
        protobufs = 'package '+package+';\n'+protobufs
    if extraProtos is not None:
        protobufs += '\n'+extraProtos
    toProtobufs = _toToProto(topname,dae)
    module = '''\
import to%(topname)sProto
import %(topname)s_pb2
''' % {'topname':topname}

    # make autogen directory if it doesn't exist
    if not os.path.exists(autogenDir):
        os.makedirs(autogenDir)

    # write files in autogen directory
    f = open(os.path.join(autogenDir,topname+'.proto'),'w')
    f.write(protobufs)
    f.close()
    f = open(os.path.join(autogenDir, 'to'+topname+'Proto.py'),'w')
    f.write(toProtobufs)
    f.close()
    f = open(os.path.join(autogenDir, '__init__.py'),'w')
    f.write(module)
    f.close()

    # call protoc to make protos
    subprocess.call(['protoc','--python_out=.','autogen/'+topname+'.proto'])
    for haskellDir in haskellDirs:
        subprocess.call(['hprotoc','-I../'+autogenDir,'--haskell_out=src',topname+'.proto'],cwd=haskellDir)
