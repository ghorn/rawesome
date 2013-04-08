import subprocess
import sys
import select

def call(args, cwd='.'):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)

    msgs = []

    while True:
        reads = [p.stdout.fileno(), p.stderr.fileno()]
        ret = select.select(reads, [], [])
    
        for fd in ret[0]:
            if fd == p.stdout.fileno():
                read = p.stdout.readline()
                sys.stdout.write(read)
                msgs.append(read)
            if fd == p.stderr.fileno():
                read = p.stderr.readline()
                sys.stderr.write(read)
                msgs.append(read)
    
        ret = p.poll()
        if ret != None:
            return (ret, (''.join(msgs)).strip())
