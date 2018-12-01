import subprocess, sys

for i in range(1000):
	subprocess.call("../../Simulation/generateMC", shell=True)
	subprocess.call("./FitFlux >> result/flux.log", shell=True)
