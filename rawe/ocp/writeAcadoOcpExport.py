# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import casadi as C

replace0 = {'real':'IntermediateState',
            'work':'_work',
            'init':'',
            'output':'_output'
            }

def writeAcadoAlgorithm(ocp, dae):
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs = [dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot]
    outputs = []
    constraintData = []
    for (lhs,comparison,rhs) in ocp._constraints:
        k = len(outputs)
        outputs.append(rhs-lhs)
        constraintData.append((k,comparison,'ALWAYS'))
    for (lhs,comparison,rhs) in ocp._constraintsEnd:
        k = len(outputs)
        outputs.append(rhs-lhs)
        constraintData.append((k,comparison,'AT_END'))
    for (lhs,comparison,rhs) in ocp._constraintsStart:
        k = len(outputs)
        outputs.append(rhs-lhs)
        constraintData.append((k,comparison,'AT_START'))
    assert len(outputs) == len(constraintData), 'the "impossible" happened'

    # add dae residual and objectives
    residual = list(dae.getResidual())

    assert len(residual) == len(dae.xNames()) + len(dae.zNames()), \
        'ocp export error: residual.size() != self.xSize() + self.zSize() ==> (%d != %d + %d)' % (len(residual),len(dae.xNames()), len(dae.zNames()))

    idx_residual = len(outputs)
    len_residual = len(residual)
    outputs += residual


    minLsq = list(ocp._minLsq)
    minLsqEndTerm = list(ocp._minLsqEndTerm)
    assert len(minLsq) == ocp._minLsq.size()
    assert len(minLsqEndTerm) == ocp._minLsqEndTerm.size()

    idx_minLsq = len(outputs)
    len_minLsq = len(minLsq)
    outputs += minLsq

    idx_minLsqEndTerm = len(outputs)
    len_minLsqEndTerm = len(minLsqEndTerm)
    outputs += minLsqEndTerm

    f = C.SXFunction( inputs, outputs )
    f.init()

    inputNames = [dae.xNames(), dae.zNames(), dae.uNames(), dae.pNames(), ['dot( '+name+' )' for name in dae.xNames()]]

    # error dictionary
    errorMap = {}
    for key,value in C.__dict__.iteritems():
        if key[0:3] == "OP_":
            errorMap[value] = key

    algStrings = []
    initializedWorkVars = set();
    def write(blah):
        algStrings.append(blah)
    def makeUnary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d );' %  replace)
    def makeInfixBinary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(work)s_%(i2)d %(op)s %(work)s_%(i3)d;' %  replace)
    def makePrefixBinary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d, %(work)s_%(i3)d );' %  replace)

    # Loop over the algorithm
    for i in range(f.getAlgorithmSize()):

        # Get the atomic operation
        op = f.getAtomicOperation(i)
        i1 = f.getAtomicOutput(i)

        replace = dict(replace0.items() + {'i1':i1}.items())
        if op != C.OP_OUTPUT:
            if i1 not in initializedWorkVars:
                initializedWorkVars.add(i1)
                replace['init'] = replace['real']+' '
        if(op==C.OP_CONST):
            replace['const'] = repr(f.getAtomicInputReal(i))
            write( '%(init)s%(work)s_%(i1)d = %(const)s;' % replace )
        else:
            i2,i3 = f.getAtomicInput(i)
            replace['i2'] = i2
            replace['i3'] = i3
            if op==C.OP_INPUT:
                #assert i2==0, "oh noes, INPUT IS MULTIDIMENSIONAL!!!"
                replace['input'] = inputNames[i2][i3]
                write( '%(init)s%(work)s_%(i1)d = %(input)s;' % replace)
            elif op==C.OP_OUTPUT:
#                assert i1==0, "oh noes, OUTPUT IS MULTIDIMENSIONAL!!!"
#                write( '%(spaces)sf << 0 == %(work)s_%(i2)d;' % replace )
                rowidx = f.output(i1).sparsity().getRow()[i3]
                colidx = f.output(i1).sparsity().col()[i3]
                assert colidx==0 and rowidx==0 and i3==0, 'non-scalars not supported in ocp constraints, colIdx: '+str(colidx)+', rowidx: '+str(rowidx)+', i3: '+str(i3)

                # different names and indexes for {constraints, dae residual, lsq, lsqEnd}
                if i1 >= idx_residual:
                    if i1 < idx_minLsq:
                        replace['output'] += '_dae_residual'
                        replace['i1'] -= idx_residual
                    else:
                        if i1 < idx_minLsqEndTerm:
                            replace['output'] += '_lsq'
                            replace['i1'] -= idx_minLsq
                        else:
                            replace['output'] += '_lsqEndTerm'
                            replace['i1'] -= idx_minLsqEndTerm

                write( '%(real)s %(output)s_%(i1)d = %(work)s_%(i2)d;' % replace )


            ########## BINARY ########
            elif op==C.OP_ADD:
                makeInfixBinary('+',replace)
            elif op==C.OP_SUB:
                makeInfixBinary('-',replace)
            elif op==C.OP_MUL:
                makeInfixBinary('*',replace)
            elif op==C.OP_SQ:
                replace['i3'] = replace['i2']
                makeInfixBinary('*',replace)
            elif op==C.OP_DIV:
                makeInfixBinary('/',replace)
            elif op==C.OP_ATAN2:
                makePrefixBinary('atan2',replace)
            elif op in [C.OP_POW,C.OP_CONSTPOW]:
                makePrefixBinary('pow',replace)

            ########## UNARY ########
            elif op==C.OP_ACOS:
                makeUnary('acos',replace)
            elif op==C.OP_ATAN:
                makeUnary('atan',replace)
            elif op==C.OP_COS:
                makeUnary('cos',replace)
            elif op==C.OP_EXP:
                makeUnary('exp',replace)
            elif op==C.OP_INV:
                makeUnary('1.0/',replace)
            elif op==C.OP_SIN:
                makeUnary('sin',replace)
            elif op==C.OP_ACOSH:
                makeUnary('acosh',replace)
            elif op==C.OP_COSH:
                makeUnary('cosh',replace)
            elif op==C.OP_FABS:
                makeUnary('fabs',replace)
            elif op==C.OP_SINH:
                makeUnary('sinh',replace)
            elif op==C.OP_ATANH:
                makeUnary('tanh',replace)
            elif op==C.OP_NEG:
                makeUnary('-',replace)
            elif op==C.OP_LOG:
                makeUnary('log',replace)
            elif op==C.OP_SQRT:
                makeUnary('sqrt',replace)
            elif op==C.OP_ASIN:
                makeUnary('asin',replace)
            elif op==C.OP_ASINH:
                makeUnary('asinh',replace)
            elif op==C.OP_TAN:
                makeUnary('tan',replace)
            elif op==C.OP_ERF:
                makeUnary('erf',replace)
            elif op==C.OP_ERFINV:
                makeUnary('erfinv',replace)
            elif op==C.OP_SIGN:
                makeUnary('sign',replace)
            elif op==C.OP_TANH:
                makeUnary('tanh',replace)
            elif op==C.OP_ASSIGN:
                makeUnary('',replace)
            elif op==C.OP_PARAMETER:
                raise KeyError('Oh man, there is a free parameter in your SXFunction')
            else:
                raise KeyError('Unknown operation: '+ errorMap[op])

    return (algStrings, constraintData)

def generateAcadoOcp(ocp, integratorOptions, ocpOptions):
    dae = ocp.dae
    #print "WARNING: RE-ENABLE PARAMETER UNSUPPORTED ASSERTION"
    assert len(dae.pNames()) == 0, 'parameters not supported by acado codegen'

    lines = []
    lines.append('/* comment the following in to enable terminal barf: */')
    
    lines.append("""\
    string path( exportDir );
    string _stdout = path + "/_stdout.txt";
    
    std::ofstream out( _stdout.c_str() );
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); // redirect std::cout to the text file
    
    Logger::instance().setLogLevel( LVL_DEBUG );
""")

    lines.append('/* differential states */')
    for name in dae.xNames():
        lines.append('DifferentialState '+name+';')
    lines.append('')
    lines.append('/* algebraic variables */')
    for name in dae.zNames():
        lines.append('AlgebraicState '+name+';')
    lines.append('')
    lines.append('/* control inputs */')
    for name in dae.uNames():
        lines.append('Control '+name+';')
    lines.append('')
    lines.append('/* parameters */')
    for name in dae.pNames():
        lines.append('Parameter '+name+';')
    lines.append('')

    (alg, constraintData) = writeAcadoAlgorithm(ocp, dae)
    lines.extend( alg )
    # differential equation
    lines.append('')
    lines.append('DifferentialEquation _differentialEquation;')
    for k in range(len(dae.xNames() + dae.zNames())):
        lines.append('_differentialEquation << 0 == _output_dae_residual_%d ;' % k)

    # minLsq
    lines.append('')
    lines.append('Function _lsqAcadoSymbolics;')
    for k in range(len(list(ocp._minLsq))):
        lines.append('_lsqAcadoSymbolics << _output_lsq_%d;' % k)

    # minLsqEndTerm
    lines.append('')
    lines.append('Function _lsqEndTermAcadoSymbolics;')
    for k in range(len(list(ocp._minLsqEndTerm))):
        lines.append('_lsqEndTermAcadoSymbolics << _output_lsqEndTerm_%d;' % k)

    # ocp
    lines.append('''
/* setup OCP */
const int N = %(N)d;
const double Ts = 1.0;
OCP _ocp(0, N * Ts, N);
_ocp.setModel( "model", "rhs", "rhsJacob" );
_ocp.setDimensions( %(nx)d, %(nx)d, %(nz)d, %(nu)d, 0, 0 );
//_ocp.subjectTo( _differentialEquation );
''' % {'nx':len(dae.xNames()), 'nz':len(dae.zNames()), 'nu':len(dae.uNames()),'N':ocp.N})

    lines.append('/* complex constraints */')
    for (k, comparison, when) in constraintData:
        if when == 'ALWAYS':
            whenStr = ''
        elif when == 'AT_END':
            whenStr = 'AT_END, '
        elif when == 'AT_START':
            whenStr = 'AT_START, '
        else:
            raise Exception('the "impossible" happened, unrecognized "when": '+str(when))
        lines.append(
            '_ocp.subjectTo( %(whenStr)s0 %(comparison)s %(output)s_%(k)d );'
            % { 'output':replace0['output'], 'comparison':comparison, 'whenStr':whenStr,
                'k':k })
    # ocp
    lines.append('')
    lines.append('/* simple constraints */')
    # the simple constraints won't work if you do
    # subjectTo( -2 <= x);
    # subjectTo(  x <= 2);
    # you have to do:
    # subjectTo( -2 <= x <= 2 );

    # bnds is a list of [(key, BLAH)]
    # where `BLAH` will be written as ocp.subjectTo( BLAH )
    # and `key` is just a key for sorting, so that we can
    # write the constraints in a nice order
    bnds = []
    for name,lbnd in ocp._lbndmap.items():
        lbnd = str(lbnd)
        if name in ocp._ubndmap:
            ubnd = str(ocp._ubndmap[name])
            bnds.append( ((name,0,1), lbnd + ' <= ' + name + ' <= ' + ubnd ) )
        else:
            bnds.append( ((name,0,1), lbnd + ' <= ' + name ) )

    for name,lbnd in ocp._lbndmapStart.items():
        lbnd = str(lbnd)
        if name in ocp._ubndmapStart:
            ubnd = str(ocp._ubndmapStart[name])
            bnds.append( ((name,1,1), 'AT_START, ' + lbnd + ' <= ' + name + ' <= ' + ubnd ) )
        else:
            bnds.append( ((name,1,1), 'AT_START, ' + lbnd + ' <= ' + name ) )

    for name,lbnd in ocp._lbndmapEnd.items():
        lbnd = str(lbnd)
        if name in ocp._ubndmapEnd:
            ubnd = str(ocp._ubndmapEnd[name])
            bnds.append( ((name,2,1), 'AT_END, ' + lbnd + ' <= ' + name + ' <= ' + ubnd ) )
        else:
            bnds.append( ((name,2,1), 'AT_END, ' + lbnd + ' <= ' + name ) )

    for name,bnd in ocp._ebndmap.items():
        bnds.append( ((name,0,0),  name + ' == ' + str(bnd) ))
    for name,bnd in ocp._ebndmapStart.items():
        bnds.append( ((name,1,0),  'AT_START, ' + name + ' == ' + str(bnd) ))
    for name,bnd in ocp._ebndmapEnd.items():
        bnds.append( ((name,2,0),  'AT_END, ' + name + ' == ' + str(bnd) ))

    bnds.sort()
    for _, bnd in bnds:
        lines.append('_ocp.subjectTo( '+bnd+' );')
    lines.append('')

    # objective
    lines.append('/* set objective */')
    lines.append('BMatrix  _W(%(size)d, %(size)d); _W.setAll( true );' % {'size': ocp._minLsq.size()})
    lines.append('BMatrix _WN(%(size)d, %(size)d); _WN.setAll( true );' % {'size': ocp._minLsqEndTerm.size()})

    lines.append('std::string _lsqExtern = "lsqExtern";')
    lines.append('std::string _lsqEndTermExtern = "lsqEndTermExtern";')
    
    # Use CasADi symbolics
    lines.append('_ocp.minimizeLSQ(        _W,        _lsqExtern);')
    lines.append('_ocp.minimizeLSQEndTerm(_WN, _lsqEndTermExtern);')
    
    # ... or use ACADO symbolics as an alternative
    lines.append('//_ocp.minimizeLSQ(        _W,        _lsqAcadoSymbolics);')
    lines.append('//_ocp.minimizeLSQEndTerm(_WN, _lsqEndTermAcadoSymbolics);')

    lines.append('''
/* setup OCPexport */
OCPexport _ocpe( _ocp );\
''')

    # write integrator options
    lines.append('\n/* integrator options */')
    for name,val in iter(sorted(integratorOptions.getAcadoOpts().items())):
        # multiply NUM_INTEGRATOR_STEPS by number of control intervals
        if name == 'NUM_INTEGRATOR_STEPS':
            val = repr(integratorOptions['NUM_INTEGRATOR_STEPS']*ocp.N)
        lines.append('_ocpe.set( '+name+', '+val+' );')

    lines.append('\n/* ocp options */')
    # write default ocp options if user hasn't set them
    for name,val in iter(sorted(ocpOptions.getAcadoOpts().items())):
        lines.append('_ocpe.set( '+name+', '+val+' );')


    lines.append('''
/* export the code */
_ocpe.exportCode( exportDir );

std::cout.rdbuf(coutbuf); //reset to standard output again

return 0
;''' )

    lines = '\n'.join(['    '+l for l in ('\n'.join(lines)).split('\n')])
    lines = '''\
#include <acado_toolkit.hpp>

extern "C" int exportOcp(const char * exportDir);

using namespace std;
USING_NAMESPACE_ACADO

int exportOcp(const char * exportDir){
''' + lines + '\n}\n'

    return lines
