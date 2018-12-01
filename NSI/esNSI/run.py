import subprocess
for i in range(0,300):
	subprocess.call("bsub python TestChi2.py %d"%i, shell=True)
