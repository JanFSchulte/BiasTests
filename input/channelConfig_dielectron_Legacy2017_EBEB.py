import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from resolution_cfg_2017 import DCB_para
nBkg = -1
dataFile = "input/eventList_ele_2017_BB.txt"

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
	nz   =  6239976 
	nsig_scale = 1./0.079562 # signal efficiency at z peak       
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result

def signalEff(mass,spin2=False):

	eff_a     = 0.5849
	eff_b     = -404.2
	eff_c     = 277.7
	eff_d     = 5.64e+04
	eff_e	  = 9.12e+04
	if spin2:
		eff_a = 0.5043
		eff_b = 1173.
		eff_c = 1.1e04
		eff_d = -1.716e05
		eff_e = 6.212e05
	
	return (eff_a+eff_b/(mass+eff_c)+eff_d/(mass*mass+eff_e))
	


	

def provideUncertainties(mass):

	result = {}

	result["sigEff"] = [1.06]
	result["massScale"] = 0.02
	result ["bkgUncert"] = 1.4
	result ["res"] = 0.0
	result ["reco"] = [0.0]

	result["bkgParams"] = {"bkg_a":0.00088373214832100731, "bkg_b":0.03817062445709217683, "bkg_c":0.07015279460424182767, "bkg_d":0.00000000000000000000, "bkg_e":0.00068184055610270096, "bkg_b2":0.02355327980473226335, "bkg_c2":0.01757757864181181198, "bkg_d2":0.05867476463668089975, "bkg_thr":0.01885201867963065200}

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
	#result['bkgUncert'] = 'dielectron_Legacy2017_EBEB' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_Legacy2017_EBEB' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_Legacy2017_EBEB' 
	result['bkgUncert'] = 'dielectron_Legacy2017_EBEB' 
	result['res'] = 'dielectron_Legacy2017_EBEB' 
	result['reco'] = 'dielectron_Legacy2017_EBEB' 
	result['bkgParams'] = 'dielectron_Legacy2017_EBEB' 


	return result

def getResolution(mass):


	CBObject = DCB_para("dcb")
	CBObject.get_value(mass,True)

	result = {}

	result["res"] = CBObject.sigma
	result["scale"] = CBObject.mean
	result["nR"] = CBObject.PowerR
	result["nL"] = CBObject.PowerL
	result["alphaL"] = CBObject.CutL
	result["alphaR"] = CBObject.CutR

	return result



def loadBackgroundShape(ws,useShapeUncert):

	bkg_a = RooRealVar('bkg_a_dielectron_Legacy2017_EBEB','bkg_a_dielectron_Legacy2017_EBEB',23.37464519)
	bkg_b = RooRealVar('bkg_b_dielectron_Legacy2017_EBEB','bkg_b_dielectron_Legacy2017_EBEB',0.001025434202)
	bkg_c = RooRealVar('bkg_c_dielectron_Legacy2017_EBEB','bkg_c_dielectron_Legacy2017_EBEB',-7.666360313e-06)
	bkg_d = RooRealVar('bkg_d_dielectron_Legacy2017_EBEB','bkg_d_dielectron_Legacy2017_EBEB',0.0)
	bkg_e = RooRealVar('bkg_e_dielectron_Legacy2017_EBEB','bkg_e_dielectron_Legacy2017_EBEB',-4.772471552)

	bkg_b2 = RooRealVar('bkg_b2_dielectron_Legacy2017_EBEB','bkg_b2_dielectron_Legacy2017_EBEB',-0.0005049133135)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_Legacy2017_EBEB','bkg_c2_dielectron_Legacy2017_EBEB',-1.532511265e-07)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_Legacy2017_EBEB','bkg_d2_dielectron_Legacy2017_EBEB',7.407462351e-12)
	bkg_thr = RooRealVar('bkg_thr_dielectron_Legacy2017_EBEB','bkg_thr_dielectron_Legacy2017_EBEB',356.4571337)

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
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_Legacy2017_EBEB','bkg_syst_a_dielectron_Legacy2017_EBEB',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_Legacy2017_EBEB','bkg_syst_b_dielectron_Legacy2017_EBEB',0.0)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

	# background shape
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_Legacy2017_EBEB",bkgParamsUncert[uncert] )
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_dielectron_Legacy2017_EBEB(mass_dielectron_Legacy2017_EBEB, bkg_a_dielectron_Legacy2017_EBEB_forUse, bkg_b_dielectron_Legacy2017_EBEB_forUse, bkg_c_dielectron_Legacy2017_EBEB_forUse,bkg_d_dielectron_Legacy2017_EBEB_forUse,bkg_e_dielectron_Legacy2017_EBEB_forUse, bkg_b2_dielectron_Legacy2017_EBEB_forUse, bkg_c2_dielectron_Legacy2017_EBEB_forUse,bkg_d2_dielectron_Legacy2017_EBEB_forUse,bkg_thr_dielectron_Legacy2017_EBEB_forUse,bkg_syst_a_dielectron_Legacy2017_EBEB,bkg_syst_b_dielectron_Legacy2017_EBEB)")
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_Legacy2017_EBEB_forUse, bkg_b_dielectron_Legacy2017_EBEB_forUse, bkg_c_dielectron_Legacy2017_EBEB_forUse,bkg_d_dielectron_Legacy2017_EBEB_forUse,bkg_e_dielectron_Legacy2017_EBEB_forUse, bkg_b2_dielectron_Legacy2017_EBEB_forUse, bkg_c2_dielectron_Legacy2017_EBEB_forUse,bkg_d2_dielectron_Legacy2017_EBEB_forUse,bkg_thr_dielectron_Legacy2017_EBEB_forUse,bkg_syst_a_dielectron_Legacy2017_EBEB,bkg_syst_b_dielectron_Legacy2017_EBEB)")
	else:
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_dielectron_Legacy2017_EBEB(mass_dielectron_Legacy2017_EBEB, bkg_a_dielectron_Legacy2017_EBEB, bkg_b_dielectron_Legacy2017_EBEB, bkg_c_dielectron_Legacy2017_EBEB,bkg_d_dielectron_Legacy2017_EBEB,bkg_e_dielectron_Legacy2017_EBEB, bkg_b2_dielectron_Legacy2017_EBEB, bkg_c2_dielectron_Legacy2017_EBEB,bkg_d2_dielectron_Legacy2017_EBEB,bkg_thr_dielectron_Legacy2017_EBEB,bkg_syst_a_dielectron_Legacy2017_EBEB,bkg_syst_b_dielectron_Legacy2017_EBEB)")
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_Legacy2017_EBEB, bkg_b_dielectron_Legacy2017_EBEB, bkg_c_dielectron_Legacy2017_EBEB,bkg_d_dielectron_Legacy2017_EBEB,bkg_e_dielectron_Legacy2017_EBEB, bkg_b2_dielectron_Legacy2017_EBEB, bkg_c2_dielectron_Legacy2017_EBEB,bkg_d2_dielectron_Legacy2017_EBEB,bkg_thr_dielectron_Legacy2017_EBEB,bkg_syst_a_dielectron_Legacy2017_EBEB,bkg_syst_b_dielectron_Legacy2017_EBEB)")
	return ws
