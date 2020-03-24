import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from resolution_cfg_2018 import DCB_para
nBkg = -1
dataFile = "input/eventList_ele_2018_BB.txt"

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
	nz   = 8601338 
	nsig_scale = 1./0.079022 # signal efficiency at z peak       
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result

def signalEff(mass,spin2=False):

	eff_a     = 0.5817
	eff_b     = -424.
	eff_c     = 392.2
	eff_d     = 4.547e+04
	eff_e	  = 1.082e+05
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

	result["bkgParams"] = {"bkg_a":0.00088467794485996127, "bkg_b":0.08973132650359907925, "bkg_c":0.08133815743836363132, "bkg_d":0.00000000000000000000, "bkg_e":0.00069138348695835727, "bkg_b2":0.02426859384368137626, "bkg_c2":0.01836827046591601315, "bkg_d2":0.06453900162196625490, "bkg_thr":0.02379977851205319905}

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
	#result['bkgUncert'] = 'dielectron_Legacy2018_EBEB' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_Legacy2018_EBEB' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_Legacy2018_EBEB' 
	result['bkgUncert'] = 'dielectron_Legacy2018_EBEB' 
	result['res'] = 'dielectron_Legacy2018_EBEB' 
	result['reco'] = 'dielectron_Legacy2018_EBEB' 
	result['bkgParams'] = 'dielectron_Legacy2018_EBEB' 

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

	bkg_a = RooRealVar('bkg_a_dielectron_Legacy2018_EBEB','bkg_a_dielectron_Legacy2018_EBEB',22.29027411)
	bkg_b = RooRealVar('bkg_b_dielectron_Legacy2018_EBEB','bkg_b_dielectron_Legacy2018_EBEB',0.0003904954765)
	bkg_c = RooRealVar('bkg_c_dielectron_Legacy2018_EBEB','bkg_c_dielectron_Legacy2018_EBEB',-6.359548657e-06)
	bkg_d = RooRealVar('bkg_d_dielectron_Legacy2018_EBEB','bkg_d_dielectron_Legacy2018_EBEB',0.0)
	bkg_e = RooRealVar('bkg_e_dielectron_Legacy2018_EBEB','bkg_e_dielectron_Legacy2018_EBEB',-4.582621254)

	bkg_b2 = RooRealVar('bkg_b2_dielectron_Legacy2018_EBEB','bkg_b2_dielectron_Legacy2018_EBEB',-0.0004969952869)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_Legacy2018_EBEB','bkg_c2_dielectron_Legacy2018_EBEB',-1.495160196e-07)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_Legacy2018_EBEB','bkg_d2_dielectron_Legacy2018_EBEB',6.865397463e-12)
	bkg_thr = RooRealVar('bkg_thr_dielectron_Legacy2018_EBEB','bkg_thr_dielectron_Legacy2018_EBEB',375.9552662)

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
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_Legacy2018_EBEB','bkg_syst_a_dielectron_Legacy2018_EBEB',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_Legacy2018_EBEB','bkg_syst_b_dielectron_Legacy2018_EBEB',0.0)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

	# background shape
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_Legacy2018_EBEB",bkgParamsUncert[uncert] )
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_dielectron_Legacy2018_EBEB(mass_dielectron_Legacy2018_EBEB, bkg_a_dielectron_Legacy2018_EBEB_forUse, bkg_b_dielectron_Legacy2018_EBEB_forUse, bkg_c_dielectron_Legacy2018_EBEB_forUse,bkg_d_dielectron_Legacy2018_EBEB_forUse,bkg_e_dielectron_Legacy2018_EBEB_forUse, bkg_b2_dielectron_Legacy2018_EBEB_forUse, bkg_c2_dielectron_Legacy2018_EBEB_forUse,bkg_d2_dielectron_Legacy2018_EBEB_forUse,bkg_thr_dielectron_Legacy2018_EBEB_forUse,bkg_syst_a_dielectron_Legacy2018_EBEB,bkg_syst_b_dielectron_Legacy2018_EBEB)")
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_Legacy2018_EBEB_forUse, bkg_b_dielectron_Legacy2018_EBEB_forUse, bkg_c_dielectron_Legacy2018_EBEB_forUse,bkg_d_dielectron_Legacy2018_EBEB_forUse,bkg_e_dielectron_Legacy2018_EBEB_forUse, bkg_b2_dielectron_Legacy2018_EBEB_forUse, bkg_c2_dielectron_Legacy2018_EBEB_forUse,bkg_d2_dielectron_Legacy2018_EBEB_forUse,bkg_thr_dielectron_Legacy2018_EBEB_forUse,bkg_syst_a_dielectron_Legacy2018_EBEB,bkg_syst_b_dielectron_Legacy2018_EBEB)")
	else:
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_dielectron_Legacy2018_EBEB(mass_dielectron_Legacy2018_EBEB, bkg_a_dielectron_Legacy2018_EBEB, bkg_b_dielectron_Legacy2018_EBEB, bkg_c_dielectron_Legacy2018_EBEB,bkg_d_dielectron_Legacy2018_EBEB,bkg_e_dielectron_Legacy2018_EBEB, bkg_b2_dielectron_Legacy2018_EBEB, bkg_c2_dielectron_Legacy2018_EBEB,bkg_d2_dielectron_Legacy2018_EBEB,bkg_thr_dielectron_Legacy2018_EBEB,bkg_syst_a_dielectron_Legacy2018_EBEB,bkg_syst_b_dielectron_Legacy2018_EBEB)")
		ws.factory("ZPrimeEleBkgPdf5::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_Legacy2018_EBEB, bkg_b_dielectron_Legacy2018_EBEB, bkg_c_dielectron_Legacy2018_EBEB,bkg_d_dielectron_Legacy2018_EBEB,bkg_e_dielectron_Legacy2018_EBEB, bkg_b2_dielectron_Legacy2018_EBEB, bkg_c2_dielectron_Legacy2018_EBEB,bkg_d2_dielectron_Legacy2018_EBEB,bkg_thr_dielectron_Legacy2018_EBEB,bkg_syst_a_dielectron_Legacy2018_EBEB,bkg_syst_b_dielectron_Legacy2018_EBEB)")
	return ws
