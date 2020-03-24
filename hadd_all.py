import os,sys
from time import sleep
import gc
import fnmatch
import subprocess

batch = True  # True

outdir = "results"
if not os.path.isdir(outdir):
	os.makedirs(outdir)

if batch:
	ch = "dimuon_Legacy2018_BB"
	funcs = [ "default", "MRBW", "MBWZg"]  # defaultlow
	masses = [200, 250, 300, 350, 400, 500,600,700,800,900,1000]  # 600
	nToys = 200
	nsigfacs = [0,3]
	script = "script_BiasTestSigBkg_batch.sh"
	for f in funcs:
		for mass in masses:
			for nsigfac in nsigfacs:
				print "\n"
				inputdir = 'results_%s'%f
				inputtag_ = "nsig_%s_fix1_nsigfac%d_%s_M%.0f_nToys%d_nEv*--Job0.root" % (ch,nsigfac,f,mass,nToys)
				for ff in os.listdir(inputdir):
					if fnmatch.fnmatch(ff, inputtag_):
						inputtag_ = str(ff)
						break
				inputtag = inputtag_.replace("Job0","Job*")
				outtag   = inputtag_.replace("--Job0","")

				cmd = 'hadd -f %s/%s %s/%s' % (outdir, outtag, inputdir, inputtag)

				print(cmd)
				subprocess.call(cmd,shell=True)
				sleep(0.5)
				sys.stdout.flush()
				gc.collect()
