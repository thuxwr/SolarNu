import subprocess, sys

core = int(sys.argv[1])
for i in range(0,300):
	#subprocess.call("./TestSNO %d %d"%(core, i), shell=True)
	subprocess.call("sleep 10", shell=True)
