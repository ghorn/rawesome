import os
import hashlib

rawesomeDataPath = os.path.expanduser("~/.rawesome")

# Given a list of (filename, string) source files, return a unique directory
# with those files written to it. If these exact files were already memoized,
# return the existing directory (possible with other stuff, like objects built my make).
def memoizeFiles(genfiles):
    # make ~/.rawesome if it doesn't exist
    if not os.path.exists(rawesomeDataPath):
        os.makedirs(rawesomeDataPath)

    try:
        assert isinstance(genfiles, list)
        for (name,src) in genfiles:
            assert isinstance(name, str)
            assert isinstance(src, str)
    except Exception:
        raise Exception("memoizeFiles needs a list of (filename, source) tuples as input")

    
    # hash the files
    genfiles.sort()
    exportpath = os.path.join(rawesomeDataPath,
                              hashlib.md5(''.join([''.join(x) for x in genfiles])).hexdigest())

    def writeFiles():
        for (filename, filestring) in genfiles:
            f = open(os.path.join(exportpath,filename),'w')
            f.write(filestring)
            f.close()

    # if no directory named by this hash exists, create it and compile the project there
    if not os.path.exists(exportpath):
        os.makedirs(exportpath)
        writeFiles()

    # if the directory already exists, check if the contents match
    else:
        unmatched = []
        for (filename, filestring) in genfiles:
            f = open(os.path.join(exportpath,filename), 'r')
            contents = f.read()
            f.close()
            if contents != filestring:
                unmatched += filename
        # if the files don't match, remove all files and write new files
        if len(unmatched) > 0:
            # remove any existing files
            for f in os.listdir(exportpath):
                os.remove(os.path.join(exportpath,f))
            # write new files
            writeFiles()

    return exportpath
