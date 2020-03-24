import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from numpy import exp
nBkg = -1
from muonResolution2017 import getResolution as getRes
dataFile = "input/event_list_2017_beee_clean_sort.txt"
def addBkgUncertPrior(ws,label,channel,uncert):

        beta_bkg = RooRealVar('beta_%s_%s'%(label,channel),'beta_%s_%s'%(label,channel),0,-5,5)
        getattr(ws,'import')(beta_bkg,ROOT.RooCmdArg())
        uncert = 1. + uncert
        bkg_kappa = RooRealVar('%s_%s_kappa'%(label,channel),'%s_%s_kappa'%(label,channel),uncert)
        bkg_kappa.setConstant()
        getattr(ws,'import')(bkg_kappa,ROOT.RooCmdArg())
        ws.factory("PowFunc::%s_%s_nuis(%s_%s_kappa, beta_%s_%s)"%(label,channel,label,channel,label,channel))
        ws.factory("prod::%s_%s_forUse(%s_%s, %s_%s_nuis)"%(label,channel,label,channel,label,channel))


def provideSignalScaling(mass,spin2=False):
	nz   =  23790                      #From Chris
	nsig_scale = 3501.8726591760296       # prescale/eff_z (561/0.1602) -->derives the lumi 
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result


def signalEff(mass,spin2=False):
	if spin2:
		eff_a = 0.211629
		eff_b = 0.124469
		eff_c = 0.885894
		eff_d = -2605.605215
		eff_e = 611.338233	
		return	eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )

	else:

##### default
		if mass <= 450:
			a =  13.39
			b =  6.696
			c = -4.855e+06
			d = -7.431e+06
			e = -108.8
			f = -1.138
			return a - b * exp( -( (mass - c) / d) ) + e * mass**f
		else:
			eff_a     =  0.3148
			eff_b     =  0.04447
			eff_c     =  1.42
			eff_d     = -5108.
			eff_e     =  713.5
			return	eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )



def recoUncert(mass):
	if mass <= 450:
		a =  13.39
		b =  6.696
		c = -4.855e+06
		d = -7.431e+06
		e = -108.8
		f = -1.138
		eff_default = a - b * exp( -( (mass - c) / d) ) + e * mass**f
	else:
		eff_a     =  0.3148
		eff_b     =  0.04447
		eff_c     =  1.42
		eff_d     = -5108.
		eff_e     =  713.5
		eff_default = eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )

	if mass <= 450:
		a =  1.33901e+01
		b =  6.69687e+00
		c = -4.85589e+06
		d = -7.43036e+06
		e = -1.14263e+02
		f = -1.15028e+00
		eff_syst= a - b * exp( -( (mass - c) / d) ) + e * mass**f
	else:
		eff_a     =  3.07958e-01
		eff_b     =  4.63280e-02
		eff_c     =  1.35632e+00
		eff_d     = -5.00475e+03
		eff_e     =  7.38088e+02
		eff_syst =  eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )



	uncertReco = 1. - eff_default/eff_syst


	uncertUp = 0.
	uncertDown = uncertReco
	
	return [1.-uncertDown,1.+uncertUp]

	

def signalEffUncert(mass):
	uncertHLT = 0.01
	uncertID  = 0.01

	uncertUp = (0.01**2 + 0.01**2)**0.5
	uncertDown = (0.01**2 + 0.01**2)**0.5
	
	return [1.-uncertDown,1.+uncertUp]

def massScaleUncert(mass):

	scale = 1.00263 -1.04029e-05*mass + 8.9214e-09*mass**2 -3.4176e-12*mass**3 +  6.07934e-16*mass**4  -3.73738e-20*mass**5

	return 1.-scale



def provideUncertainties(mass):

	result = {}

	result["reco"] = recoUncert(mass)
	result["sigEff"] = signalEffUncert(mass)
	result["massScale"] = massScaleUncert(mass)
	result ["bkgUncert"] = 1.4	
	result ["res"] = 0.085
	result["bkgParams"] = {"bkg_a":0.00149864671535,"bkg_b":0.00653883209809,"bkg_c":0.0463578422191,"bkg_d":0.045733835952,"bkg_e":0.00127869845332, "bkg_a2":0.00495212596353,"bkg_b2":0.0110368347341,"bkg_c2":0.0360014117393,"bkg_e2":0.00657290957759}
	return result

def provideCorrelations():
	result = {}
	''' Combine correlates uncertainties that have the same name. So wa hve to adjust the names to achieve what we want. 
		1) put the full channel name. That will make it uncorrelated with all other channels
		2) keep the channel name but remove the last bit: will correlate between the two subcategories within a year
		3) Just keep the dimuon or dielectron name, so we correlate between the years
		4) To correlate some specific combination of uncertainties, come up with a name and add it to all releavent channel configs
	'''
	#result['sigEff'] = 'dimuon' 
	#result['massScale'] = 'dimuon' 
	#result['bkgUncert'] = 'dimuon_2017_BE' 
	#result['res'] = 'dimuon' 
	#result['bkgParams'] = 'dimuon_2017_BE' 
	result['sigEff'] = 'dimuon' 
	result['massScale'] = 'dimuon' 
	result['bkgUncert'] = 'dimuon_2017_BE' 
	result['res'] = 'dimuon_2017_BE' 
	result['reco'] = 'dimuon_BE' 
	result['bkgParams'] = 'dimuon_2017_BE' 

	return result

def provideUncertaintiesCI(mass):

	result = {}

	result["trig"] = 1.007
	result["zPeak"] = 1.05
	result["xSecOther"] = 1.07
	result["jets"] = 1.5
	result["lumi"] = 1.025
	result["massScale"] = 0.0 ## dummy value
	result["stats"] = 0.0 ## dummy value
	result["res"] = 0.0 ## dummy value
	result["pdf"] = 0.0 ## dummy value
	result["ID"] = 0.0 ## dummy value
	result["PU"] = 0.0 ## dummy value
	return result



def getResolution(mass):
	result = {}
	params = getRes(mass)
	result['alphaL'] = params['alphaL']['BE']
	result['alphaR'] = params['alphaR']['BE']
	result['nL'] = params['nL']['BE']
	result['nR'] = params['nR']['BE']
	result['res'] = params['sigma']['BE']
	result['scale'] = params['scale']['BE']
	return result


def loadBackgroundShape(ws,useShapeUncert=False):


	bkg_a = RooRealVar('bkg_a_dimuon_2017_BE','bkg_a_dimuon_2017_BE', 20.2294579295)
	bkg_b = RooRealVar('bkg_b_dimuon_2017_BE','bkg_b_dimuon_2017_BE',-0.00169012794522)
	bkg_c = RooRealVar('bkg_c_dimuon_2017_BE','bkg_c_dimuon_2017_BE',8.06798491306e-08)
	bkg_d = RooRealVar('bkg_d_dimuon_2017_BE','bkg_d_dimuon_2017_BE',-1.03307501705e-11)
	bkg_e = RooRealVar('bkg_e_dimuon_2017_BE','bkg_e_dimuon_2017_BE',-3.85977306121)
	bkg_a2 = RooRealVar('bkg_a2_dimuon_2017_BE','bkg_a2_dimuon_2017_BE', 7.45339497885)
	bkg_b2 = RooRealVar('bkg_b2_dimuon_2017_BE','bkg_b2_dimuon_2017_BE',-0.0119993639555)
	bkg_c2 = RooRealVar('bkg_c2_dimuon_2017_BE','bkg_c2_dimuon_2017_BE', 5.65225224556e-06)
	bkg_e2 = RooRealVar('bkg_e2_dimuon_2017_BE','bkg_e2_dimuon_2017_BE',-1.18673996788)

	bkg_a.setConstant()
	bkg_b.setConstant()
	bkg_c.setConstant()
	bkg_d.setConstant()
	bkg_e.setConstant()
	bkg_a2.setConstant()
	bkg_b2.setConstant()
	bkg_c2.setConstant()
	bkg_e2.setConstant()
	getattr(ws,'import')(bkg_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_b,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_c,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_d,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_e,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_a2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_b2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_c2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_e2,ROOT.RooCmdArg())
	
	# background systematics
	bkg_syst_a = RooRealVar('bkg_syst_a','bkg_syst_a',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',0.000)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.00016666666666)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

	if useShapeUncert: 
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dimuon_2017_BE",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dimuon_2017_BE",0.1 )

	# background shape
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2017_BE(mass_dimuon_2017_BE, bkg_a_dimuon_2017_BE_forUse, bkg_b_dimuon_2017_BE_forUse, bkg_c_dimuon_2017_BE_forUse,bkg_d_dimuon_2017_BE_forUse,bkg_e_dimuon_2017_BE_forUse, bkg_a2_dimuon_2017_BE_forUse, bkg_b2_dimuon_2017_BE_forUse, bkg_c2_dimuon_2017_BE_forUse,bkg_e2_dimuon_2017_BE_forUse,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2017_BE_forUse, bkg_b_dimuon_2017_BE_forUse, bkg_c_dimuon_2017_BE_forUse,bkg_d_dimuon_2017_BE_forUse,bkg_e_dimuon_2017_BE_forUse,bkg_a2_dimuon_2017_BE_forUse, bkg_b2_dimuon_2017_BE_forUse, bkg_c2_dimuon_2017_BE_forUse,bkg_e2_dimuon_2017_BE_forUse,bkg_syst_a,bkg_syst_b)")	
	else:
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2017_BE(mass_dimuon_2017_BE, bkg_a_dimuon_2017_BE, bkg_b_dimuon_2017_BE, bkg_c_dimuon_2017_BE,bkg_d_dimuon_2017_BE,bkg_e_dimuon_2017_BE, bkg_a2_dimuon_2017_BE, bkg_b2_dimuon_2017_BE, bkg_c2_dimuon_2017_BE,bkg_e2_dimuon_2017_BE,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2017_BE, bkg_b_dimuon_2017_BE, bkg_c_dimuon_2017_BE,bkg_d_dimuon_2017_BE,bkg_e_dimuon_2017_BE,bkg_a2_dimuon_2017_BE, bkg_b2_dimuon_2017_BE, bkg_c2_dimuon_2017_BE,bkg_e2_dimuon_2017_BE,bkg_syst_a,bkg_syst_b)")	

	return ws
