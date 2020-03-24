import os,sys
from time import sleep
import gc
import subprocess

batch = True  # True  False

template = '''
Executable = script_BiasTestSigBkg_batch.sh
Arguments  = dimuon_Legacy2018_BB %s %d %d %d %d
Log        = logs/%s.log
Output     = logs/%s.out
Error      = logs/%s.error
+JobFlavour = "tomorrow"
Queue
'''




if not os.path.isdir("logs"):
	os.makedirs("logs")

if batch:
	ch = "dimuon_Legacy2018_BB"
	funcs = [ "default", "MRBW", "MBWZg"]    #  "defaultlow", 
	#funcs = [ "MRBW"]    #  "defaultlow", 
	masses = [200, 250, 300, 350, 400, 500,600,700,800,900,1000]  # 600
	#masses = [1000]  # 600
	nToys = 200
	nJobs = int(2000/nToys)
	nsigfacs = [0,3]
	script = "script_BiasTestSigBkg_batch.sh"
	for f in funcs:
		for nsigfac in nsigfacs:
			for mass in masses:
				for jobId in range(nJobs):
					cname = "%s_nsigfac%d_%s_M%.0f_nToys%d--Job%d" %(ch,nsigfac,f,mass,nToys,jobId)
				        text_file = open("condor.sub", "w")
        				text_file.write(template % (f,mass,nToys,nsigfac,jobId,cname,cname,cname) )
        				#text_file.write(template % (f,mass,nToys,nsigfac,jobId) )
        				text_file.close()

					
					print("condor_submit condor.sub")
					cmd = "condor_submit condor.sub"
					subprocess.call(cmd,shell=True)
					sleep(0.5)
					sys.stdout.flush()
					gc.collect()


else:
	ranges = [
		("dimuon_Legacy2018_BB", "default", 400, 2, 0, 0),
		("dimuon_Legacy2018_BB", "MRBW",    200, 2, 3, 0),
		("dimuon_Legacy2018_BB", "MBWZg",   400, 2, 3, 0)
	]
	for ch, f, mass, nToys, nsigfac, jobId in ranges:
		cname = "%s_nsigfac%d_%s_M%.0f_nToys%d--Job%d" %(ch,nsigfac,f,mass,nToys,jobId)

		cmd = "python runBiasTestSigBkg.py --channel %s -f %s --mass %d --nToys %d --nSigFac %d --jobId %d >&logs/log-%s.log&" % (ch, f, mass, nToys, nsigfac, jobId, cname)
		print cname
		os.system(cmd)

		sys.stdout.flush()
		sleep(0.05)
