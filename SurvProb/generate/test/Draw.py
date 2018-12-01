from ROOT import *
import os, sys

energy = []
SurvProb = []
#canvas = TCanvas()
graph = TGraph(2000)
for core in range(0,100):
	file = open("./tmp/core%d.dat"%core, "r")
	j = 0
	while True:
		line = file.readline()
		if not line:
			break
		energy_tmp, SurvProb_tmp = [float(i) for i in line.split()]
		k = 20 * core + j
		graph.SetPoint(k,energy_tmp,SurvProb_tmp)
		j = j+1
	file.close()
	os.remove("./tmp/core%d.dat"%core)

#graph.Draw()
graph.SaveAs("graph.root")




