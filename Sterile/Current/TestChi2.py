import subprocess, sys


core = int(sys.argv[1])

#f = open('./tmp/core%d.dat'%core, 'r')
#lines = f.readlines()
#lastline = lines[-1]

#ncore = lastline.split()[1]
#start = int(ncore) + 1

for i in range(0,300):
	subprocess.call("./test %d %d"%(core, i), shell=True)
