import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from resolution_cfg_2018 import DCB_para
nBkg = -1
#nBkg = -1

dataFile = "input/eventList_ele_2018_BE.txt"
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
	nz   = 3401386 
	nsig_scale = 1./0.031377
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result

def signalEff(mass,spin2=False):
	eff_a     = 0.01461
	eff_b     = 479.6
	eff_c     = 635.3
	eff_d     = -1.071e+05
	eff_e	  = 8.345e+04
	eff_f 	  = 1.302e+07
	eff_g 	  = 2.337e+07
	
	if spin2:
		eff_a = 0.06212
		eff_b = -7.192
		eff_c = 56.72
		eff_d = -43.69
		eff_e = 822.9
		eff_f = 3.579e08
		eff_g = 3.048e09
		
	return (eff_a+eff_b/(mass+eff_c)+eff_d/(mass*mass+eff_e))+eff_f/(mass**3+eff_g)

		

def provideUncertainties(mass):

	result = {}

	result["sigEff"] = [1.08] # must be list in case the uncertainty is asymmetric
	result["massScale"] = 0.01
	result ["bkgUncert"] = 1.4
	result ["res"] = 0.0
	result ["reco"] = [0.0]

	result["bkgParams"] = {"bkg_a":0.0309509874258958595,"bkg_b":0.0103550392308542384,"bkg_c":0.0,"bkg_d":0.0,"bkg_e":0.0033432250380655902,"bkg_a2":0.0115578563069132900,"bkg_b2":0.0057746173631073837,"bkg_c2":0.0281542841879044679,"bkg_d2":0.0491110796660534921,"bkg_e2":0.0020474522359512676}

	# result["bkgParams"] = {"bkg_a":0.04905186032576128,"bkg_b":0.010813178009224173,"bkg_c":0.0,"bkg_d":0.0,"bkg_e":0.0036873576105158007,"bkg_a2":0.01304294839773059,"bkg_b2":0.0077626144239507314,"bkg_c2":0.03245544113278542,"bkg_d2":0.05457592336699882,"bkg_e2":0.0022675053701084673}
	#result["bkgParams"] = {"bkg_a":0.23843520903507034,"bkg_b":0.042557205520631004,"bkg_c":0.4772108713152145,"bkg_d":0.20039400197000984,"bkg_e":0.02121430756489364,"bkg_a2":0.21935669378841924,"bkg_b2":0.05489687902965316,"bkg_c2":0.11981929926738805,"bkg_d2":0.15078950936605048,"bkg_e2":0.023322536781016343}
	return result

def provideCorrelations():
	result = {}
	''' Combine correlates uncertainties that have the same name. So wa hve to adjust the names to achieve what we want. 
		1) put the full channel name. That will make it uncorrelated with all other channels
		2) keep the channel name but remove the last bit: will correlate between the two subcategories within a year
		3) Just keep the dimuon or dielectron name, so we correlate between the years
		4) To correlate some specific combination of uncertainties, come up with a name and add it to all releavent channel configs
	'''
	#result['sigEff'] = 'dielectron' 
	#result['massScale'] = 'dielectron' 
	#result['bkgUncert'] = 'dielectron_2018_EBEE' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_2018_EBEE' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_2018_EBEE' 
	result['bkgUncert'] = 'dielectron_2018_EBEE' 
	result['res'] = 'dielectron_2018_EBEE' 
	result['reco'] = 'dielectron_2018_EBEE' 
	result['bkgParams'] = 'dielectron_2018_EBEE' 

	return result


def getResolution(mass):
	CBObject = DCB_para("dcb")
	CBObject.get_value(mass,False)

	result = {}

	result["res"] = CBObject.sigma
	result["scale"] = CBObject.mean
	result["nR"] = CBObject.PowerR
	result["nL"] = CBObject.PowerL
	result["alphaL"] = CBObject.CutL
	result["alphaR"] = CBObject.CutR
	if result["nR"] < 0:
		result["nR"] = 0.
	return result

def loadBackgroundShape(ws,useShapeUncert=False):
	# /afs/cern.ch/user/w/wenxing/public/ZprimeRun2/bkg_fit/bkg_fit_newFR.txt
	bkg_a = RooRealVar('bkg_a_dielectron_2018_EBEE','bkg_a_dielectron_2018_EBEE',1.58898e+00)
	bkg_b = RooRealVar('bkg_b_dielectron_2018_EBEE','bkg_b_dielectron_2018_EBEE',-3.61960e-03)
	bkg_c = RooRealVar('bkg_c_dielectron_2018_EBEE','bkg_c_dielectron_2018_EBEE',0)
	bkg_d = RooRealVar('bkg_d_dielectron_2018_EBEE','bkg_d_dielectron_2018_EBEE',0)
	bkg_e = RooRealVar('bkg_e_dielectron_2018_EBEE','bkg_e_dielectron_2018_EBEE',-2.80437e+00)

	bkg_a2 = RooRealVar('bkg_a2_dielectron_2018_EBEE','bkg_a2_dielectron_2018_EBEE',4.28942e+00)
	bkg_b2 = RooRealVar('bkg_b2_dielectron_2018_EBEE','bkg_b2_dielectron_2018_EBEE',-2.44553e-03)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_2018_EBEE','bkg_c2_dielectron_2018_EBEE',1.64722e-07)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_2018_EBEE','bkg_d2_dielectron_2018_EBEE',-1.41340e-11)
	bkg_e2 = RooRealVar('bkg_e2_dielectron_2018_EBEE','bkg_e2_dielectron_2018_EBEE',-3.35137e+00)

	# bkg_a = RooRealVar('bkg_a_dielectron_2018_EBEE','bkg_a_dielectron_2018_EBEE',1.10142e+00)
	# bkg_b = RooRealVar('bkg_b_dielectron_2018_EBEE','bkg_b_dielectron_2018_EBEE',-3.69464e-03)
	# bkg_c = RooRealVar('bkg_c_dielectron_2018_EBEE','bkg_c_dielectron_2018_EBEE',0)
	# bkg_d = RooRealVar('bkg_d_dielectron_2018_EBEE','bkg_d_dielectron_2018_EBEE',0)
	# bkg_e = RooRealVar('bkg_e_dielectron_2018_EBEE','bkg_e_dielectron_2018_EBEE',-2.78286e+00)
	# bkg_a2 = RooRealVar('bkg_a2_dielectron_2018_EBEE','bkg_a2_dielectron_2018_EBEE',4.14381e+00)
	# bkg_b2 = RooRealVar('bkg_b2_dielectron_2018_EBEE','bkg_b2_dielectron_2018_EBEE',-2.56393e-03)
	# bkg_c2 = RooRealVar('bkg_c2_dielectron_2018_EBEE','bkg_c2_dielectron_2018_EBEE',2.42376e-07)
	# bkg_d2 = RooRealVar('bkg_d2_dielectron_2018_EBEE','bkg_d2_dielectron_2018_EBEE',-2.67039e-11)
	# bkg_e2 = RooRealVar('bkg_e2_dielectron_2018_EBEE','bkg_e2_dielectron_2018_EBEE',-3.37982e+00)


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
	
        
	bkg_a2.setConstant()
	bkg_b2.setConstant()
	bkg_c2.setConstant()
	bkg_d2.setConstant()
	bkg_e2.setConstant()
	getattr(ws,'import')(bkg_a2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_b2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_c2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_d2,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_e2,ROOT.RooCmdArg())
	
	# background systematics
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_2018_EBEE','bkg_syst_a_dielectron_2018_EBEE',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_2018_EBEE','bkg_syst_b_dielectron_2018_EBEE',0.000)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.000017)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_2018_EBEE",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dielectron_2018_EBEE",0.1 )

	# background shape
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2018_EBEE(mass_dielectron_2018_EBEE, bkg_a_dielectron_2018_EBEE_forUse, bkg_b_dielectron_2018_EBEE_forUse, bkg_c_dielectron_2018_EBEE_forUse,bkg_d_dielectron_2018_EBEE_forUse,bkg_e_dielectron_2018_EBEE_forUse,bkg_a2_dielectron_2018_EBEE_forUse, bkg_b2_dielectron_2018_EBEE_forUse, bkg_c2_dielectron_2018_EBEE_forUse,bkg_d2_dielectron_2018_EBEE_forUse,bkg_e2_dielectron_2018_EBEE_forUse,bkg_syst_a_dielectron_2018_EBEE,bkg_syst_b_dielectron_2018_EBEE)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2018_EBEE_forUse, bkg_b_dielectron_2018_EBEE_forUse, bkg_c_dielectron_2018_EBEE_forUse,bkg_d_dielectron_2018_EBEE_forUse,bkg_e_dielectron_2018_EBEE_forUse,bkg_a2_dielectron_2018_EBEE_forUse, bkg_b2_dielectron_2018_EBEE_forUse, bkg_c2_dielectron_2018_EBEE_forUse,bkg_d2_dielectron_2018_EBEE_forUse,bkg_e2_dielectron_2018_EBEE_forUse,bkg_syst_a_dielectron_2018_EBEE,bkg_syst_b_dielectron_2018_EBEE)")		
	else:	
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2018_EBEE(mass_dielectron_2018_EBEE, bkg_a_dielectron_2018_EBEE, bkg_b_dielectron_2018_EBEE, bkg_c_dielectron_2018_EBEE,bkg_d_dielectron_2018_EBEE,bkg_e_dielectron_2018_EBEE,bkg_a2_dielectron_2018_EBEE, bkg_b2_dielectron_2018_EBEE, bkg_c2_dielectron_2018_EBEE,bkg_d2_dielectron_2018_EBEE,bkg_e2_dielectron_2018_EBEE,bkg_syst_a_dielectron_2018_EBEE,bkg_syst_b_dielectron_2018_EBEE)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2018_EBEE, bkg_b_dielectron_2018_EBEE, bkg_c_dielectron_2018_EBEE,bkg_d_dielectron_2018_EBEE,bkg_e_dielectron_2018_EBEE,bkg_a2_dielectron_2018_EBEE, bkg_b2_dielectron_2018_EBEE, bkg_c2_dielectron_2018_EBEE,bkg_d2_dielectron_2018_EBEE,bkg_e2_dielectron_2018_EBEE,bkg_syst_a_dielectron_2018_EBEE,bkg_syst_b_dielectron_2018_EBEE)")		
	return ws
