import subprocess
for i in range(200):
	subprocess.call("bsub python DrawProb.py %d"%i, shell=True)
