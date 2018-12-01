#!/usr/bin/python
import subprocess
for i in range(0,100):
	subprocess.call("bsub python TestSurv.py %d" %i, shell=True)
