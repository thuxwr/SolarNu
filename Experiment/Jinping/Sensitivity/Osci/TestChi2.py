import subprocess, sys

core = int(sys.argv[1])
for i in range(0,300):
	subprocess.call("./full %d %d"%(core, i), shell=True)
