import subprocess, sys

theta = int(sys.argv[1])
for i in range(200):
	subprocess.call("./DrawProb %d %d"%(theta, i), shell=True)
