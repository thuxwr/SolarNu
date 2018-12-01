import subprocess

#for i in range(0,25):
#	subprocess.call("bsub -q short -n 12 python run.py %d"%i, shell=True)
for i in range(0,300):
	subprocess.call("bsub -n 1 python TestChi2.py %d"%i, shell=True)
