import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from muonResolution2018 import getResolution as getRes
nBkg = -1
dataFile = "input/event_list_2018_bb_clean_sort.txt"

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
	nz   =  27379                      #from Chris
	nsig_scale = 4317.789291882556       # prescale/eff_z (500/0.1158) -->derives the lumi 
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result
	

def signalEff(mass,spin2=False):


	if spin2:
		eff_a = 1.020382
		eff_b = -1166.881533
		eff_c = 1468.989496
		eff_d = 0.000044
	 	return eff_a + eff_b / (mass + eff_c) - mass*eff_d
	else:
#### default	
		if mass <= 600:
			a = 2.14
			b = 0.1286
			c = 110.6
			d = 22.44
			e = -2.366
			f = -0.03382
			from math import exp
			return a - b * exp( -(mass - c) / d ) + e * mass**f
		else:

			eff_a     =   5.18
			eff_b     =  -5.845e+04
			eff_c     =  1.157e+04
			eff_d     =  0.0002255

			return	eff_a + eff_b / (mass + eff_c) - mass*eff_d
		
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
def massScaleUncert(mass):

	scale = 0.999032 + 3.36979e-06*mass -3.4122e-09*mass**2 + 1.62541e-12*mass**3  - 3.12864e-16*mass**4 + 2.18417e-20*mass**5

	return 1.-scale	




def provideUncertainties(mass):

	result = {}

	result["reco"] = recoUncert(mass)
	result["sigEff"] = signalEffUncert(mass)
	result["massScale"] = massScaleUncert(mass)
	result["bkgUncert"] = 1.4
	result["res"] = 0.085

	result["bkgParams"] = {"bkg_a":0.00476859850283119658, "bkg_b":0.00722495482027206773, "bkg_c":0.03683614775776260925, "bkg_d":0.00000000000000000000, "bkg_e":0.00221513405352908232, "bkg_b2":0.04667954324047498410, "bkg_c2":0.02414196780603639417, "bkg_d2":0.09715793303505616818, "bkg_thr":0.00756266886512532016}

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
	#result['bkgUncert'] = 'dimuon_Legacy2018ScaleDown_BB' 
	#result['res'] = 'dimuon' 
	#result['bkgParams'] = 'dimuon_Legacy2018ScaleDown_BB' 
	result['sigEff'] = 'dimuon' 
	result['massScale'] = 'dimuon' 
	result['bkgUncert'] = 'dimuon_Legacy2018ScaleDown_BB' 
	result['res'] = 'dimuon_Legacy2018ScaleDown_BB' 
	result['reco'] = 'dimuon_BB' 
	result['bkgParams'] = 'dimuon_Legacy2018ScaleDown_BB' 


	return result


def loadBackgroundShape(ws,useShapeUncert=False):

	bkg_a = RooRealVar('bkg_a_dimuon_Legacy2018ScaleDown_BB','bkg_a_dimuon_Legacy2018ScaleDown_BB',5.549345124)
	bkg_b = RooRealVar('bkg_b_dimuon_Legacy2018ScaleDown_BB','bkg_b_dimuon_Legacy2018ScaleDown_BB',-0.006861620005)
	bkg_c = RooRealVar('bkg_c_dimuon_Legacy2018ScaleDown_BB','bkg_c_dimuon_Legacy2018ScaleDown_BB',8.148765728e-06)
	bkg_d = RooRealVar('bkg_d_dimuon_Legacy2018ScaleDown_BB','bkg_d_dimuon_Legacy2018ScaleDown_BB',0.0)
	bkg_e = RooRealVar('bkg_e_dimuon_Legacy2018ScaleDown_BB','bkg_e_dimuon_Legacy2018ScaleDown_BB',-1.765679588)

	bkg_b2 = RooRealVar('bkg_b2_dimuon_Legacy2018ScaleDown_BB','bkg_b2_dimuon_Legacy2018ScaleDown_BB',-0.0003585771134)
	bkg_c2 = RooRealVar('bkg_c2_dimuon_Legacy2018ScaleDown_BB','bkg_c2_dimuon_Legacy2018ScaleDown_BB',-1.687626869e-07)
	bkg_d2 = RooRealVar('bkg_d2_dimuon_Legacy2018ScaleDown_BB','bkg_d2_dimuon_Legacy2018ScaleDown_BB',6.549990739e-12)
	bkg_thr = RooRealVar('bkg_thr_dimuon_Legacy2018ScaleDown_BB','bkg_thr_dimuon_Legacy2018ScaleDown_BB',390.745276)

	bkg_a.setConstant()
	bkg_b.setConstant()
	bkg_c.setConstant()
	bkg_d.setConstant()
	bkg_e.setConstant()
	getattr(ws,'import')(bkg_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_b,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_c,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_d,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_e,ROOT.RooCmdArg())

	bkg_b2.setConstant()
	bkg_c2.setConstant()
	bkg_d2.setConstant()
	bkg_thr.setConstant()
	getattr(ws,'import')(bkg_b2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_c2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_d2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_thr,ROOT.RooCmdArg())

	# background systematics
	bkg_syst_a = RooRealVar('bkg_syst_a_dimuon_Legacy2018ScaleDown_BB','bkg_syst_a_dimuon_Legacy2018ScaleDown_BB',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dimuon_Legacy2018ScaleDown_BB','bkg_syst_b_dimuon_Legacy2018ScaleDown_BB',0.0)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

	# background shape
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dimuon_Legacy2018ScaleDown_BB",bkgParamsUncert[uncert] )
		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_dimuon_Legacy2018ScaleDown_BB(mass_dimuon_Legacy2018ScaleDown_BB, bkg_a_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_b_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_c_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_d_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_e_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_b2_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_c2_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_d2_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_thr_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_syst_a_dimuon_Legacy2018ScaleDown_BB,bkg_syst_b_dimuon_Legacy2018ScaleDown_BB)")
		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_b_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_c_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_d_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_e_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_b2_dimuon_Legacy2018ScaleDown_BB_forUse, bkg_c2_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_d2_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_thr_dimuon_Legacy2018ScaleDown_BB_forUse,bkg_syst_a_dimuon_Legacy2018ScaleDown_BB,bkg_syst_b_dimuon_Legacy2018ScaleDown_BB)")
	else:
		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_dimuon_Legacy2018ScaleDown_BB(mass_dimuon_Legacy2018ScaleDown_BB, bkg_a_dimuon_Legacy2018ScaleDown_BB, bkg_b_dimuon_Legacy2018ScaleDown_BB, bkg_c_dimuon_Legacy2018ScaleDown_BB,bkg_d_dimuon_Legacy2018ScaleDown_BB,bkg_e_dimuon_Legacy2018ScaleDown_BB, bkg_b2_dimuon_Legacy2018ScaleDown_BB, bkg_c2_dimuon_Legacy2018ScaleDown_BB,bkg_d2_dimuon_Legacy2018ScaleDown_BB,bkg_thr_dimuon_Legacy2018ScaleDown_BB,bkg_syst_a_dimuon_Legacy2018ScaleDown_BB,bkg_syst_b_dimuon_Legacy2018ScaleDown_BB)")
		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_Legacy2018ScaleDown_BB, bkg_b_dimuon_Legacy2018ScaleDown_BB, bkg_c_dimuon_Legacy2018ScaleDown_BB,bkg_d_dimuon_Legacy2018ScaleDown_BB,bkg_e_dimuon_Legacy2018ScaleDown_BB, bkg_b2_dimuon_Legacy2018ScaleDown_BB, bkg_c2_dimuon_Legacy2018ScaleDown_BB,bkg_d2_dimuon_Legacy2018ScaleDown_BB,bkg_thr_dimuon_Legacy2018ScaleDown_BB,bkg_syst_a_dimuon_Legacy2018ScaleDown_BB,bkg_syst_b_dimuon_Legacy2018ScaleDown_BB)")
	return ws
