import casadi as C

replace0 = {'real':'IntermediateState',
            'work':'_work',
            'init':'',
            'output':'_output'
            }

def writeAcadoAlgorithm(ocp, dae):
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs = [dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot]
    outputs = [ocp._minLsq, ocp._minLsqEndTerm]
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
    assert len(outputs)-2 == len(constraintData), 'the "impossible" happened'

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
                if i1 == 0:
                    write( '_obj << %(work)s_%(i2)d;' % replace )
                elif i1 == 1:
                    write( '_objEnd << %(work)s_%(i2)d;' % replace )
                else:
                    rowidx = f.output(i1).sparsity().getRow()[i3]
                    colidx = f.output(i1).sparsity().col()[i3]
                    assert colidx==0 and rowidx==0 and i3==0, 'non-scalars not supported in ocp constraints, colIdx: '+str(colidx)+', rowidx: '+str(rowidx)+', i3: '+str(i3)
    
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

def generateAcadoOcp(ocp, acadoOptions):
    dae = ocp._dae
    #print "WARNING: RE-ENABLE PARAMETER UNSUPPORTED ASSERTION"
    assert len(dae.pNames()) == 0, 'parameters not supported by acado codegen'

    lines = []

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

    # constraints and objective algorithm
    lines.append('/* setup objective function and constraints */')
    lines.append('Function _obj;')
    lines.append('Function _objEnd;')
    lines.append('')
    (alg, constraintData) = writeAcadoAlgorithm(ocp, dae)
    lines.extend( alg )
    lines.append('')
    lines.append('''\
/* setup OCP */
OCP _ocp(0, N * Ts, N);
_ocp.setModel( "model", "rhs", "rhs_jac" );
_ocp.setDimensions( %(nx)d, %(nx)d, %(nz)d, %(nup)d );\
''' % {'nx':len(dae.xNames()), 'nz':len(dae.zNames()), 'nup':len(dae.uNames())+len(dae.pNames())})

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
    bnds = []
    for name,bnd in ocp._lbndmap.items():
        bnds.append( (name,           '', '>=', bnd) )
    for name,bnd in ocp._lbndmapStart.items():
        bnds.append( (name, 'AT_START, ', '>=', bnd) )
    for name,bnd in ocp._lbndmapEnd.items():
        bnds.append( (name,   'AT_END, ', '>=', bnd) )

    for name,bnd in ocp._ubndmap.items():
        bnds.append( (name,           '', '<=', bnd) )
    for name,bnd in ocp._ubndmapStart.items():
        bnds.append( (name, 'AT_START, ', '<=', bnd) )
    for name,bnd in ocp._ubndmapEnd.items():
        bnds.append( (name,   'AT_END, ', '<=', bnd) )

    for name,bnd in ocp._ebndmap.items():
        bnds.append( (name,           '', '==', bnd) )
    for name,bnd in ocp._ebndmapStart.items():
        bnds.append( (name, 'AT_START, ', '==', bnd) )
    for name,bnd in ocp._ebndmapEnd.items():
        bnds.append( (name,   'AT_END, ', '==', bnd) )

    bnds.sort()
    for name, whenStr, comparison, bnd in bnds:
        lines.append('_ocp.subjectTo( %(whenStr)s%(name)s %(comparison)s %(bnd)s );'
                     % {'whenStr':whenStr, 'name':name, 'comparison':comparison, 'bnd':bnd})
    lines.append('')

    # objective
    lines.append('/* set objective */')
    lines.append('ExportVariable  _SS( "SS", %(size)d, %(size)d);' % {'size': ocp._minLsq.size()})
    lines.append('ExportVariable _SSN("SSN", %(size)d, %(size)d);' % {'size': ocp._minLsqEndTerm.size()})
    lines.append('_ocp.minimizeLSQ(_SS, _obj);')
    lines.append('_ocp.minimizeLSQEndTerm(_SSN, _objEnd);')
    
    lines.append('''
/* setup OCPexport */
OCPexport _ocpe( _ocp );\
''')

    for name,val in acadoOptions:
        lines.append('_ocpe.set( '+name+', '+val+' );')

    lines.append('''
/* export the code */
_ocpe.exportCode( exportDir );

return 0;''' )

    lines = '\n'.join(['    '+l for l in ('\n'.join(lines)).split('\n')])
    lines = '''\
#include <acado_toolkit.hpp>
#include <ocp_export.hpp>

extern "C" int exportOcp(int N, double Ts, const char * exportDir);

using namespace std;
USING_NAMESPACE_ACADO

int exportOcp(int N, double Ts, const char * exportDir){
''' + lines + '\n}'
    return lines
