import subprocess

for i in range(1,360):
	subprocess.call("cp flux0.dat flux%d.dat"%i, shell=True)
	#subprocess.call("mkdir %d"%i, shell=True)
