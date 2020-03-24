import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from muonResolution2017 import getResolution as getRes
nBkg = -1

dataFile = "input/event_list_2017_bb_clean_sort.txt"
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
	nz   =  17116                      #From Chris
	nsig_scale = 4865.568083261058       # prescale/eff_z (561/0.1153) -->derives the lumi 
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result
	

def signalEff(mass,spin2=False):

	# spin 2 to be updated!
	if spin2:
		eff_a = 1.020382
		eff_b = -1166.881533
		eff_c = 1468.989496
		eff_d = 0.000044
	 	return eff_a + eff_b / (mass + eff_c) - mass*eff_d
	else:
#### default	
		if mass <= 600:
			a = 2.13
			b = 0.1313
			c = 110.9
			d = 20.31
			e = -2.387
			f = -0.03616
			from math import exp
			return a - b * exp( -(mass - c) / d ) + e * mass**f
		else:

			eff_a     =   4.931
			eff_b     =  -5.55e+04
			eff_c     =  1.157e+04
			eff_d     =  0.0002108

			return	eff_a + eff_b / (mass + eff_c) - mass*eff_d
	
def massScaleUncert(mass):

	scale = 1.00037 -1.70128e-06*mass + 3.09916e-09*mass**2 -1.32999e-12*mass**3 +  2.99434e-16*mass**4 + -2.28065e-20*mass**5

	return 1.-scale	
def recoUncert(mass):
	uncertReco = 0.01

	uncertDown = 1. - uncertReco
	uncertUp =   1. 
	
	return [uncertDown,uncertUp]

def signalEffUncert(mass):
	uncertID = 0.01
	uncertHLT = 0.01

	uncertDown = 1. - (uncertID**2 + uncertHLT**2 )**0.5
	uncertUp =   1. + (uncertID**2 + uncertHLT**2 )**0.5
	
	return [uncertDown,uncertUp]



def provideUncertainties(mass):

	result = {}

	result["reco"] = recoUncert(mass)
	result["sigEff"] = signalEffUncert(mass)
	result["massScale"] = massScaleUncert(mass)
	result["bkgUncert"] = 1.4
	result["res"] = 0.085
	result["bkgParams"] = {"bkg_a":0.00137974806878,"bkg_b":0.0115227291489,"bkg_c":0.0479441531366,"bkg_d":0.0658594210773,"bkg_e":0.00107617884738, "bkg_a2":0.00518651113976,"bkg_b2":0.0115872553246,"bkg_c2":0.0307170160895,"bkg_e2":0.0068019195307}
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
	#result['bkgUncert'] = 'dimuon_2017_BB' 
	#result['res'] = 'dimuon' 
	#result['bkgParams'] = 'dimuon_2017_BB' 
	result['sigEff'] = 'dimuon' 
	result['massScale'] = 'dimuon' 
	result['bkgUncert'] = 'dimuon_2017_BB' 
	result['res'] = 'dimuon_BB' 
	result['reco'] = 'dimuon_2017_BB' 
	result['bkgParams'] = 'dimuon_2017_BB' 


	return result

def provideUncertaintiesCI(mass):

	result = {}

	result["trig"] = 1.003
	result["zPeak"] = 1.05
	result["xSecOther"] = 1.07
	result["jets"] = 1.5
	result["lumi"] = 1.025
	result["stats"] = 0.0 ##dummy values
	result["massScale"] = 0.0 ##dummy values
	result["res"] = 0.0 ## dummy values
	result["pdf"] = 0.0 ## dummy values
	result["ID"] = 0.0 ## dummy values
	result["PU"] = 0.0 ## dummy values

	return result



def getResolution(mass):
	result = {}
	params = getRes(mass)
	result['alphaL'] = params['alphaL']['BB']
	result['alphaR'] = params['alphaR']['BB']
	result['nL'] = params['nL']['BB']
	result['nR'] = params['nR']['BB']
	result['res'] = params['sigma']['BB']
	result['scale'] = params['scale']['BB']

	return result




def loadBackgroundShape(ws,useShapeUncert=False):

	bkg_a = RooRealVar('bkg_a_dimuon_2017_BB','bkg_a_dimuon_2017_BB', 20.08967958)
	bkg_b = RooRealVar('bkg_b_dimuon_2017_BB','bkg_b_dimuon_2017_BB',-0.000726469552191)
	bkg_c = RooRealVar('bkg_c_dimuon_2017_BB','bkg_c_dimuon_2017_BB',-5.30536993745e-08)
	bkg_d = RooRealVar('bkg_d_dimuon_2017_BB','bkg_d_dimuon_2017_BB',-4.87301932e-12)
	bkg_e = RooRealVar('bkg_e_dimuon_2017_BB','bkg_e_dimuon_2017_BB',-4.04060321483)
	bkg_a2 = RooRealVar('bkg_a2_dimuon_2017_BB','bkg_a2_dimuon_2017_BB', 8.74986890797)
	bkg_b2 = RooRealVar('bkg_b2_dimuon_2017_BB','bkg_b2_dimuon_2017_BB',-0.0146386650943)
	bkg_c2 = RooRealVar('bkg_c2_dimuon_2017_BB','bkg_c2_dimuon_2017_BB', 8.66470108875e-06)
	bkg_e2 = RooRealVar('bkg_e2_dimuon_2017_BB','bkg_e2_dimuon_2017_BB',-1.43336705959)
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
	bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',0.0)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.00016666666666)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

	# background shape
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dimuon_2017_BB",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dimuon_2017_BB",0.1 )

		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2017_BB(mass_dimuon_2017_BB, bkg_a_dimuon_2017_BB_forUse, bkg_b_dimuon_2017_BB_forUse, bkg_c_dimuon_2017_BB_forUse,bkg_d_dimuon_2017_BB_forUse,bkg_e_dimuon_2017_BB_forUse,bkg_a2_dimuon_2017_BB_forUse, bkg_b2_dimuon_2017_BB_forUse, bkg_c2_dimuon_2017_BB_forUse,bkg_e2_dimuon_2017_BB_forUse,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2017_BB_forUse, bkg_b_dimuon_2017_BB_forUse, bkg_c_dimuon_2017_BB_forUse,bkg_d_dimuon_2017_BB_forUse,bkg_e_dimuon_2017_BB_forUse, bkg_a2_dimuon_2017_BB_forUse, bkg_b2_dimuon_2017_BB_forUse, bkg_c2_dimuon_2017_BB_forUse,bkg_e2_dimuon_2017_BB,bkg_syst_a,bkg_syst_b)")		

	else:
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2017_BB(mass_dimuon_2017_BB, bkg_a_dimuon_2017_BB, bkg_b_dimuon_2017_BB, bkg_c_dimuon_2017_BB,bkg_d_dimuon_2017_BB,bkg_e_dimuon_2017_BB,bkg_a2_dimuon_2017_BB, bkg_b2_dimuon_2017_BB, bkg_c2_dimuon_2017_BB,bkg_e2_dimuon_2017_BB,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2017_BB, bkg_b_dimuon_2017_BB, bkg_c_dimuon_2017_BB,bkg_d_dimuon_2017_BB,bkg_e_dimuon_2017_BB, bkg_a2_dimuon_2017_BB, bkg_b2_dimuon_2017_BB, bkg_c2_dimuon_2017_BB,bkg_e2_dimuon_2017_BB,bkg_syst_a,bkg_syst_b)")		
	return ws
