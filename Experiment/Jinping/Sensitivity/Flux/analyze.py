import math

f = open('./result/Flux.dat', 'r')

nflux = []
nfluxerr = []
nrelaerr = []


nline = 0
nMC = -1

for line in f.readlines():
	nline += 1
	m = nline%19
	if m == 1:
		nMC += 1
		nflux.append([])
		nfluxerr.append([])
		nrelaerr.append([])
	if m==1 or m==2 or m==5 or m==6 or m==11 or m==16 or m==18:
		continue
	else:
		# Order: pp, pep, Be7__862, B8, N13, O15, Kr85, Bi210, C11, C14, Tl208, ExtTl208 

		
		nflux[nMC].append(line.split()[1])
		nfluxerr[nMC].append(line.split()[3])
		nrelaerr[nMC].append(line.split()[4])

comp = 5
flux = 0
fluxerr = 0

errnum = 0
for i in range(nMC+1):
	R = float(nfluxerr[i][comp])/float(nflux[i][comp])
	if R > 8:
		errnum += 1
		continue
	flux += float(nflux[i][comp])
	fluxerr += float(nfluxerr[i][comp])
	if R > 4:
		print(i*19)

print(errnum)
flux /= (nMC+1-errnum)
fluxerr /= (nMC+1-errnum)
print(flux)
print(fluxerr)

relaerr = fluxerr / flux

print(relaerr)


