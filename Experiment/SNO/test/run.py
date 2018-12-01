import subprocess, multiprocessing, sys

job = int(sys.argv[1])

def run(num):
	subprocess.call("python TestChi2.py %d"%num, shell=True)

p = multiprocessing.Pool(12)
for i in range(job*12,12+job*12):
	p.apply_async(run, args=(i,))
p.close()
p.join()

print('Finished')
