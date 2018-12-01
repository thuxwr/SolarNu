import subprocess

for i in range(30):
	subprocess.call("bsub -n 12 -q hpc_linux python run.py %d"%i, shell=True)
