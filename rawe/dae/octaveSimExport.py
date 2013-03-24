import casadi as C

replace0 = {'real':'',
            'work':'work',
            'init':''
            }
def writeAcadoAlgorithm(dae, fm):
#    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
#    inputs = [dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot]
#    outputs = [dae[name] for name in dae.outputNames()]

    inputNames = ['x','u','p']

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
    for i in range(fm.getAlgorithmSize()):
      
        # Get the atomic operation
        op = fm.getAtomicOperation(i)
        i1 = fm.getAtomicOutput(i)

        replace = dict(replace0.items() + {'i1':i1}.items())
        if op != C.OP_OUTPUT:
            if i1 not in initializedWorkVars:
                initializedWorkVars.add(i1)
                replace['init'] = replace['real']+''
        if(op==C.OP_CONST):
            replace['const'] = repr(fm.getAtomicInputReal(i))
            write( '%(init)s%(work)s_%(i1)d = %(const)s;' % replace )
        else:
            i2,i3 = fm.getAtomicInput(i)
            replace['i2'] = i2
            replace['i3'] = i3
            if op==C.OP_INPUT:
                #assert i2==0, "oh noes, INPUT IS MULTIDIMENSIONAL!!!"
                replace['input'] = inputNames[i2]
                replace['i3'] += 1
                write( '%(init)s%(work)s_%(i1)d = %(input)s(%(i3)d);' % replace)
            elif op==C.OP_OUTPUT:
                rowidx = fm.output(i1).sparsity().getRow()[i3]
                colidx = fm.output(i1).sparsity().col()[i3]
                replace['row'] = rowidx+1
                replace['col'] = colidx+1
                if i1==0:
                    # blah
                    write( 'f(%(row)d,%(col)d) = %(work)s_%(i2)d;' % replace )
                elif i1==1:
                    # mass matrix
                    write( 'MM(%(row)d,%(col)d) = %(work)s_%(i2)d;' % replace )
                else:
                    raise Exception('too many outputs....')
            
            ########## BINARY ########
            elif op==C.OP_ADD:
                makeInfixBinary('+',replace)
            elif op==C.OP_SUB:
                makeInfixBinary('-',replace)
            elif op==C.OP_MUL:
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

    return algStrings

def generateOctaveSim(dae, functionName):
#    return C.SXFunction( C.daeIn( x=self.xVec(),
#                                  z=C.veccat([self.zVec(),xdot]),
#                                  p=C.veccat([self.uVec(),self.pVec()])
#                                  ),
#                         C.daeOut( alg=f, ode=xdot) )

    # get the residual fg(xdot,x,z)
    fg = dae.getResidual()

    # take the jacobian w.r.t. xdot,z
    z = dae.zVec()
    jac = C.jacobian(fg,C.veccat([dae.xDotVec(), z]))

    # make sure that it was linear in {xdot,z}, i.e. the jacobian is not a function of {xdot,z}
    testJac = C.SXFunction([dae.xVec(),dae.uVec(),dae.pVec()], [jac])
    testJac.init()
    assert len(testJac.getFree()) == 0, "can't generate octave sim, jacobian a function of {xdot,z}"

    # it was linear, so export the jacobian
    fg_fun = C.SXFunction([dae.xVec(),dae.zVec(),dae.uVec(),dae.pVec(),dae.xDotVec()], [fg])
    fg_fun.init()

    # get the constant term
    [fg_zero] = fg_fun.eval([dae.xVec(),0*dae.zVec(),dae.uVec(),dae.pVec(),0*dae.xDotVec()])
    testFun = C.SXFunction([dae.xVec(),dae.uVec(),dae.pVec()], [jac])
    testFun.init()
    assert len(testFun.getFree()) == 0, "can't generate octave sim, function line linear in {xdot,z}"

    fm = C.SXFunction([dae.xVec(), dae.uVec(), dae.pVec()],[fg_zero, jac])
    fm.init()
    lines = []
    lines.append('function [f,MM] = '+functionName+'_modelAndJacob(x,u,p)')
    lines.append('')
    lines.append('MM = zeros'+str(jac.shape)+';')
    lines.append('f = zeros('+str(fg_zero.size())+',1);')
    lines.append('')
    # dae residual
    lines.extend( writeAcadoAlgorithm(dae, fm) )
    lines.append('')
    lines.append('end\n')

    return ('\n'.join(lines))
