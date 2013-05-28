import os
import hashlib
import shutil
import tempfile
import multiprocessing

rawesomeDataPath = os.path.expanduser("~/.rawesome")

def makeJobs():
    return '-j'+str(multiprocessing.cpu_count())

# Given a recursive dict filename:source, return a unique directory with
# those files written to it. If these exact files were already memoized,
# return the existing directory (possible with other stuff, like objects built my make).
def memoizeFiles(genfiles,prefix=''):
    # make ~/.rawesome if it doesn't exist
    if not os.path.exists(rawesomeDataPath):
        os.makedirs(rawesomeDataPath)

    try:
        def doAssertions(blah):
            assert isinstance(blah, dict), 'blah not a dict'
            for name,src in blah.items():
                assert isinstance(name, str), "name not a string: "+str(name)
                if isinstance(src,dict):
                    doAssertions(src)
                else:
                    assert isinstance(src, str), "src not a string: "+str(src)
        doAssertions(genfiles)
    except Exception as e:
        raise Exception("memoizeFiles needs a (possibly recursive) dict of filename:source input\n"+\
                        str(e))

    # hash the files
    exportpath = _getExportPath(genfiles,prefix)

    writeDifferentFiles(exportpath, genfiles)

    return exportpath

# Try to write all the files.
# If a file already exists, only overwrite if they are different.
def writeDifferentFiles(path,gfs):
    assert isinstance(gfs,dict), 'not a dict'
    assert isinstance(path,str), 'not a string'

    if not os.path.exists(path):
        os.makedirs(path)

    for filename,filestring in gfs.items():
        newpath = os.path.join(path,filename)
        if isinstance(filestring,dict):
            # recursive directory
            writeDifferentFiles(newpath,filestring)
        else:
            # normal file
            def writeFile():
                f = open(newpath,'w')
                f.write(filestring)
                f.close()
            try:
                with open(newpath, 'r') as f:
                    existingfile = f.read()
                if existingfile != filestring:
                    writeFile()
            except:
                writeFile()

def _getExportPath(genfiles,prefix):
    flattened = flattenFileDict(genfiles)
    flattened.sort()
    return os.path.join(rawesomeDataPath,
                        prefix+hashlib.md5(''.join([''.join(x) for x in flattened])).hexdigest())

def flattenFileDict(fd,prefix=''):
        assert isinstance(fd,dict)
        ret = []
        for name,src in fd.items():
            if isinstance(src,str):
                ret.append( (os.path.join(prefix,name),src) )
            else:
                ret += flattenFileDict(src, prefix=os.path.join(prefix,name))
        return ret

def directoryToDict(top):
    ret = {}
    for name in os.listdir(top):
        filename = os.path.join(top, name)
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                ret[name] = f.read()
        elif os.path.isdir(filename) and not os.path.islink(filename):
            ret[name] = directoryToDict(filename)
        else:
            raise Exception(filename+" is not a file or directory (aka WTF IS GOING ON)")
    return ret


# Input a command which takes a directory as argument.
# Create a temp directory at tmpdirPath and call command(tmpdirPath).
# Return all files created in tmpdirPath as a recursive dictionary and remove tmpdirPath.
def withTempdir(command):
    shutil.rmtree # make sure this exists before making the tmp dir
    # make temporary directory
    tmppath = tempfile.mkdtemp()
    try:
        # run the user's command
        command(tmppath)

        # grab all files generated in the temp directory
        allfiles = directoryToDict(tmppath)

    finally:
        shutil.rmtree(tmppath)

    return allfiles

def writeCCode(f, name):
    def callme(tmpdir):
        f.generateCode( os.path.join(tmpdir,'generatedCode.c') )
    namespace = "ALRIGHT_GUYS_LETS_EXPORT_"+name
    codestring = withTempdir(callme)['generatedCode.c']
    codestring = '\n'.join(['  '+c for c in codestring.split('\n')])
    codestring = "namespace "+namespace+"{\n"+\
                 codestring +\
                 "} // namespace "+ namespace +"\n\n"

    real = 'double'
    args0 = ['const '+real+'* x'+str(k) for k in range(f.getNumInputs())]
    args0 += [real+'* r'+str(k) for k in range(f.getNumOutputs())]
    args0 = '(' + ', '.join(args0) + ')'
    args1 = ['x'+str(k) for k in range(f.getNumInputs())]
    args1 += ['r'+str(k) for k in range(f.getNumOutputs())]
    args1 = '(' + ', '.join(args1) + ')'

    proto0 = 'void '+name+args0
    fun0 = proto0+'{\n' + \
           '  '+namespace+'::evaluate'+args1+';\n}\n'

    proto1 = 'int '+name+'Wrap(const '+real+'** x, '+real+'** r)'
    fun1 = proto1+'{\n' + \
           '  return '+namespace+'::evaluateWrap(x, r);\n}\n'
    codestring += fun0 + fun1
    header = '''\
#ifndef __%(name)s_H__
#define __%(name)s_H__

#ifdef __cplusplus
extern "C" {
#endif

%(protos)s

#ifdef __cplusplus
}
#endif

#endif // __%(name)s_H__
''' % {'name':name, 'protos':proto0 + ';\n' + proto1 + ';'}
    return (codestring, header)
