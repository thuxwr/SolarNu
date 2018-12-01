import subprocess, sys

core = int(sys.argv[1])
initeng = 0.05
for i in range(0,20):
	eng = initeng + core * 0.2 + i * 0.05 * 0.2
	subprocess.call("./TestProb %d %lf"%(core, eng), shell=True)
	
