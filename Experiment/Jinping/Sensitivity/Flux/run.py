import subprocess, multiprocessing, sys, random, time

job = int(sys.argv[1])

def run(num):
	for j in range(200):
		k = random.random()
		time.sleep(k)
		subprocess.call("../../Simulation/generateMC /work/wangzhe9_work/neutrino/Experiment/Jinping/Sensitivity/Flux/result/tmp/%d"%num, shell=True)
		subprocess.call("./FitFlux %d >> /work/wangzhe9_work/neutrino/Experiment/Jinping/Sensitivity/Flux/result/tmp/flux%d.dat"%(num,num), shell=True)
		#subprocess.call("./FitFlux %d"%num, shell=True)

p = multiprocessing.Pool(12)
for i in range(job*12, 12+job*12):
	p.apply_async(run, args=(i,))
p.close()
p.join()

