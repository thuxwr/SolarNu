import subprocess, sys

core = int(sys.argv[1])
for i in range(0,300):
	subprocess.call("./test %d %d"%(core, i), shell=True)
