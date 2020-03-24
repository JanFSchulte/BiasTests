import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from resolution_cfg_2018 import DCB_para
#nBkg = -1
nBkg = -1 
#### to be updated!
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

 	result["bkgParams"] = {"bkg_a":0.0100436247057868674,"bkg_b":0.0160198977749077889,"bkg_c":0.0,"bkg_d":0.4285896314803248064,"bkg_e":0.0015902481687985130,"bkg_a2":0.0095227728113440203,"bkg_b2":0.0050040561594527320,"bkg_c2":0.0,"bkg_d2":0.0288032225740758180,"bkg_e2":0.0010481925785027532}  

 	# result["bkgParams"] = {"bkg_a2":0.010627691816876696,"bkg_b2":0.017969606992918632,"bkg_c2":0.,"bkg_d2":0.23222798587022248,"bkg_e2":0.0018184836580724154,"bkg_a":0.010537979706711924,"bkg_b":0.007702240238193461,"bkg_c":0.0,"bkg_d":0.06350974182012484,"bkg_e":0.0012679963411284364}	
 	#result["bkgParams"] = {"bkg_a":0.24391455519910865,"bkg_b":0.03963999446145325,"bkg_c":44.56738363916309,"bkg_d":0.35495154369249793,"bkg_e":0.001701404156849987,"bkg_a2":0.1575451905377131,"bkg_b2":0.04733064954639549,"bkg_c2":0.93272729080101,"bkg_d2":0.6290363984210399,"bkg_e2":0.0024446111862020054}	
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
	#result['bkgUncert'] = 'dielectron_2018_EBEB' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_2018_EBEB' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_2018_EBEB' 
	result['bkgUncert'] = 'dielectron_2018_EBEB' 
	result['res'] = 'dielectron_2018_EBEB' 
	result['reco'] = 'dielectron_2018_EBEB' 
	result['bkgParams'] = 'dielectron_2018_EBEB' 

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
	# /afs/cern.ch/user/w/wenxing/public/ZprimeRun2/bkg_fit/bkg_fit_newFR.txt
	bkg_a = RooRealVar('bkg_a_dielectron_2018_EBEB','bkg_a_dielectron_2018_EBEB',2.95704e+00)
	bkg_b = RooRealVar('bkg_b_dielectron_2018_EBEB','bkg_b_dielectron_2018_EBEB',-1.82734e-03)
	bkg_c = RooRealVar('bkg_c_dielectron_2018_EBEB','bkg_c_dielectron_2018_EBEB',0.)
	bkg_d = RooRealVar('bkg_d_dielectron_2018_EBEB','bkg_d_dielectron_2018_EBEB',-1.20075e-11)
	bkg_e = RooRealVar('bkg_e_dielectron_2018_EBEB','bkg_e_dielectron_2018_EBEB',-3.65880e+00)

	bkg_a2 = RooRealVar('bkg_a2_dielectron_2018_EBEB','bkg_a2_dielectron_2018_EBEB',3.11424e+00)
	bkg_b2 = RooRealVar('bkg_b2_dielectron_2018_EBEB','bkg_b2_dielectron_2018_EBEB',-1.11682e-03)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_2018_EBEB','bkg_b2_dielectron_2018_EBEB',0)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_2018_EBEB','bkg_c2_dielectron_2018_EBEB',-6.47681e-12)
	bkg_e2 = RooRealVar('bkg_e2_dielectron_2018_EBEB','bkg_d2_dielectron_2018_EBEB',-3.75369e+00)

	# bkg_a = RooRealVar('bkg_a_dielectron_2018_EBEB','bkg_a_dielectron_2018_EBEB',3.25988e+00)
	# bkg_b = RooRealVar('bkg_b_dielectron_2018_EBEB','bkg_b_dielectron_2018_EBEB',-1.74119e-03)
	# bkg_c = RooRealVar('bkg_c_dielectron_2018_EBEB','bkg_c_dielectron_2018_EBEB',0.)
	# bkg_d = RooRealVar('bkg_d_dielectron_2018_EBEB','bkg_d_dielectron_2018_EBEB',-2.19961e-11)
	# bkg_e = RooRealVar('bkg_e_dielectron_2018_EBEB','bkg_e_dielectron_2018_EBEB',-3.68347e+00)

	# bkg_a2 = RooRealVar('bkg_a2_dielectron_2018_EBEB','bkg_a2_dielectron_2018_EBEB',3.27596e+00)
	# bkg_b2 = RooRealVar('bkg_b2_dielectron_2018_EBEB','bkg_b2_dielectron_2018_EBEB',-1.11506e-03)
	# bkg_c2 = RooRealVar('bkg_c2_dielectron_2018_EBEB','bkg_b2_dielectron_2018_EBEB',0)
	# bkg_d2 = RooRealVar('bkg_d2_dielectron_2018_EBEB','bkg_c2_dielectron_2018_EBEB',-6.69998e-12)
	# bkg_e2 = RooRealVar('bkg_e2_dielectron_2018_EBEB','bkg_d2_dielectron_2018_EBEB',-3.74979e+00)
	

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
	
	# background systematics
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_2018_EBEB','bkg_syst_a_dielectron_2018_EBEB',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_2018_EBEB','bkg_syst_b_dielectron_2018_EBEB',0.000)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.000017)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_2018_EBEB",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dielectron_2018_EBEB",0.1 )
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2018_EBEB(mass_dielectron_2018_EBEB, bkg_a_dielectron_2018_EBEB_forUse, bkg_b_dielectron_2018_EBEB_forUse, bkg_c_dielectron_2018_EBEB_forUse,bkg_d_dielectron_2018_EBEB_forUse,bkg_e_dielectron_2018_EBEB_forUse,bkg_a2_dielectron_2018_EBEB_forUse, bkg_b2_dielectron_2018_EBEB_forUse, bkg_c2_dielectron_2018_EBEB_forUse,bkg_d2_dielectron_2018_EBEB_forUse,bkg_e2_dielectron_2018_EBEB_forUse,bkg_syst_a_dielectron_2018_EBEB,bkg_syst_b_dielectron_2018_EBEB)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2018_EBEB_forUse, bkg_b_dielectron_2018_EBEB_forUse, bkg_c_dielectron_2018_EBEB_forUse,bkg_d_dielectron_2018_EBEB_forUse,bkg_e_dielectron_2018_EBEB_forUse,bkg_a2_dielectron_2018_EBEB_forUse, bkg_b2_dielectron_2018_EBEB_forUse, bkg_c2_dielectron_2018_EBEB_forUse,bkg_d2_dielectron_2018_EBEB_forUse,bkg_e2_dielectron_2018_EBEB_forUse,bkg_syst_a_dielectron_2018_EBEB,bkg_syst_b_dielectron_2018_EBEB)")		
	else:
	# background shape
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2018_EBEB(mass_dielectron_2018_EBEB, bkg_a_dielectron_2018_EBEB, bkg_b_dielectron_2018_EBEB, bkg_c_dielectron_2018_EBEB,bkg_d_dielectron_2018_EBEB,bkg_e_dielectron_2018_EBEB,bkg_a2_dielectron_2018_EBEB, bkg_b2_dielectron_2018_EBEB, bkg_c2_dielectron_2018_EBEB,bkg_d2_dielectron_2018_EBEB,bkg_e2_dielectron_2018_EBEB,bkg_syst_a_dielectron_2018_EBEB,bkg_syst_b_dielectron_2018_EBEB)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2018_EBEB, bkg_b_dielectron_2018_EBEB, bkg_c_dielectron_2018_EBEB,bkg_d_dielectron_2018_EBEB,bkg_e_dielectron_2018_EBEB,bkg_a2_dielectron_2018_EBEB, bkg_b2_dielectron_2018_EBEB, bkg_c2_dielectron_2018_EBEB,bkg_d2_dielectron_2018_EBEB,bkg_e2_dielectron_2018_EBEB,bkg_syst_a_dielectron_2018_EBEB,bkg_syst_b_dielectron_2018_EBEB)")		
	return ws
