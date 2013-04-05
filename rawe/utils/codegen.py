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
def memoizeFiles(genfiles):
    # make ~/.rawesome if it doesn't exist
    if not os.path.exists(rawesomeDataPath):
        os.makedirs(rawesomeDataPath)

    try:
        def doAssertions(blah):
            assert isinstance(blah, dict), 'blah not a dict'
            for name,src in blah.items():
                assert isinstance(name, str), "name not a string"
                if isinstance(src,dict):
                    doAssertions(src)
                else:
                    assert isinstance(src, str), "src not a string"
        doAssertions(genfiles)
    except Exception:
        raise Exception("memoizeFiles needs a (possibly recursive) dict of filename:source input")

    # hash the files
    exportpath = getExportPath(genfiles)

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

def getExportPath(genfiles):
    flattened = flattenFileDict(genfiles)
    flattened.sort()
    return os.path.join(rawesomeDataPath,
                        hashlib.md5(''.join([''.join(x) for x in flattened])).hexdigest())

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
