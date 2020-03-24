import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from  resolution_cfg_2016 import DCB_para
nBkg = -1
#nBkg = -1

dataFile = "input/eventList_ele_2016_BE.txt"
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
	nz   = 1967746 
	nsig_scale = 1./0.030449 
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result

def signalEff(mass,spin2=False):
	eff_a     = -0.03377
	eff_b     = 735.1
	eff_c     = 1318.
	eff_d     = -8.889e04
	eff_e	  = 7.572e04
	eff_f 	  = 1.424e07
	eff_g 	  = 2.342e07
	
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

	result["bkgParams"] = {"bkg_a":0.0190105863031046554,"bkg_b":0.0093443269569989402,"bkg_c":0.0,"bkg_d":0.0,"bkg_e":0.0033740617608506215,"bkg_a2":0.0100566883019371639,"bkg_b2":0.0052304916121128026,"bkg_c2":0.0190423236910167958,"bkg_d2":0.0304093554971315329,"bkg_e2":0.0020731630263761043}

	# result["bkgParams"] = {"bkg_a":0.040838258275886684,"bkg_b":0.0096497267045262,"bkg_c":0.0,"bkg_d":0.0,"bkg_e":0.0036235797545109355,"bkg_a2":0.010891764745699611,"bkg_b2":0.00722767661043595,"bkg_c2":0.02658968410157474,"bkg_d2":0.041935937439232705,"bkg_e2":0.00212398172820708}
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
	#result['bkgUncert'] = 'dielectron_2016_EBEE' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_2016_EBEE' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_2016_EBEE' 
	result['bkgUncert'] = 'dielectron_2016_EBEE' 
	result['reco'] = 'dielectron_2016_EBEE' 
	result['res'] = 'dielectron_2016_EBEE' 
	result['bkgParams'] = 'dielectron_2016_EBEE' 

	return result

def getResolution(mass):
	result = {}
	dcb = DCB_para('dcb')
	dcb.get_value(mass,False) ##True for BB, False for BE

	result['nL'] = dcb.PowerL
	result['nR'] = dcb.PowerR
	result['alphaL'] = dcb.CutL
	result['alphaR'] = dcb.CutR
	result['res'] = dcb.sigma
	result['scale'] = dcb.mean

	return result
def loadBackgroundShape(ws,useShapeUncert=False):
	# /afs/cern.ch/user/w/wenxing/public/ZprimeRun2/bkg_fit/bkg_fit_newFR.txt
	bkg_a = RooRealVar('bkg_a_dielectron_2016_EBEE','bkg_a_dielectron_2016_EBEE',2.52685e+00)
	bkg_b = RooRealVar('bkg_b_dielectron_2016_EBEE','bkg_b_dielectron_2016_EBEE',-4.08176e-03)
	bkg_c = RooRealVar('bkg_c_dielectron_2016_EBEE','bkg_c_dielectron_2016_EBEE',0)
	bkg_d = RooRealVar('bkg_d_dielectron_2016_EBEE','bkg_d_dielectron_2016_EBEE',0)
	bkg_e = RooRealVar('bkg_e_dielectron_2016_EBEE','bkg_e_dielectron_2016_EBEE',-2.73118e+00)

	bkg_a2 = RooRealVar('bkg_a2_dielectron_2016_EBEE','bkg_a2_dielectron_2016_EBEE',4.79464e+00)
	bkg_b2 = RooRealVar('bkg_b2_dielectron_2016_EBEE','bkg_b2_dielectron_2016_EBEE',-2.70986e-03)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_2016_EBEE','bkg_c2_dielectron_2016_EBEE',2.51703e-07)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_2016_EBEE','bkg_d2_dielectron_2016_EBEE',-2.38106e-11)
	bkg_e2 = RooRealVar('bkg_e2_dielectron_2016_EBEE','bkg_e2_dielectron_2016_EBEE',-3.22868e+00)

	# bkg_a = RooRealVar('bkg_a_dielectron_2016_EBEE','bkg_a_dielectron_2016_EBEE',1.26059e+00)
	# bkg_b = RooRealVar('bkg_b_dielectron_2016_EBEE','bkg_b_dielectron_2016_EBEE',-4.13106e-03)
	# bkg_c = RooRealVar('bkg_c_dielectron_2016_EBEE','bkg_c_dielectron_2016_EBEE',0)
	# bkg_d = RooRealVar('bkg_d_dielectron_2016_EBEE','bkg_d_dielectron_2016_EBEE',0)
	# bkg_e = RooRealVar('bkg_e_dielectron_2016_EBEE','bkg_e_dielectron_2016_EBEE',-2.71784e+00)
	# bkg_a2 = RooRealVar('bkg_a2_dielectron_2016_EBEE','bkg_a2_dielectron_2016_EBEE',4.72748e+00)
	# bkg_b2 = RooRealVar('bkg_b2_dielectron_2016_EBEE','bkg_b2_dielectron_2016_EBEE',-2.59526e-03)
	# bkg_c2 = RooRealVar('bkg_c2_dielectron_2016_EBEE','bkg_c2_dielectron_2016_EBEE',2.84047e-07)
	# bkg_d2 = RooRealVar('bkg_d2_dielectron_2016_EBEE','bkg_d2_dielectron_2016_EBEE',-3.34267e-11)
	# bkg_e2 = RooRealVar('bkg_e2_dielectron_2016_EBEE','bkg_e2_dielectron_2016_EBEE',-3.41510e+00)


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
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_2016_EBEE','bkg_syst_a_dielectron_2016_EBEE',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_2016_EBEE','bkg_syst_b_dielectron_2016_EBEE',0.000)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.000017)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_2016_EBEE",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dielectron_2016_EBEE",0.1 )

	# background shape
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2016_EBEE(mass_dielectron_2016_EBEE, bkg_a_dielectron_2016_EBEE_forUse, bkg_b_dielectron_2016_EBEE_forUse, bkg_c_dielectron_2016_EBEE_forUse,bkg_d_dielectron_2016_EBEE_forUse,bkg_e_dielectron_2016_EBEE_forUse,bkg_a2_dielectron_2016_EBEE_forUse, bkg_b2_dielectron_2016_EBEE_forUse, bkg_c2_dielectron_2016_EBEE_forUse,bkg_d2_dielectron_2016_EBEE_forUse,bkg_e2_dielectron_2016_EBEE_forUse,bkg_syst_a_dielectron_2016_EBEE,bkg_syst_b_dielectron_2016_EBEE)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2016_EBEE_forUse, bkg_b_dielectron_2016_EBEE_forUse, bkg_c_dielectron_2016_EBEE_forUse,bkg_d_dielectron_2016_EBEE_forUse,bkg_e_dielectron_2016_EBEE_forUse,bkg_a2_dielectron_2016_EBEE_forUse, bkg_b2_dielectron_2016_EBEE_forUse, bkg_c2_dielectron_2016_EBEE_forUse,bkg_d2_dielectron_2016_EBEE_forUse,bkg_e2_dielectron_2016_EBEE_forUse,bkg_syst_a_dielectron_2016_EBEE,bkg_syst_b_dielectron_2016_EBEE)")		
	else:	
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2016_EBEE(mass_dielectron_2016_EBEE, bkg_a_dielectron_2016_EBEE, bkg_b_dielectron_2016_EBEE, bkg_c_dielectron_2016_EBEE,bkg_d_dielectron_2016_EBEE,bkg_e_dielectron_2016_EBEE,bkg_a2_dielectron_2016_EBEE, bkg_b2_dielectron_2016_EBEE, bkg_c2_dielectron_2016_EBEE,bkg_d2_dielectron_2016_EBEE,bkg_e2_dielectron_2016_EBEE,bkg_syst_a_dielectron_2016_EBEE,bkg_syst_b_dielectron_2016_EBEE)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2016_EBEE, bkg_b_dielectron_2016_EBEE, bkg_c_dielectron_2016_EBEE,bkg_d_dielectron_2016_EBEE,bkg_e_dielectron_2016_EBEE,bkg_a2_dielectron_2016_EBEE, bkg_b2_dielectron_2016_EBEE, bkg_c2_dielectron_2016_EBEE,bkg_d2_dielectron_2016_EBEE,bkg_e2_dielectron_2016_EBEE,bkg_syst_a_dielectron_2016_EBEE,bkg_syst_b_dielectron_2016_EBEE)")		
	return ws
