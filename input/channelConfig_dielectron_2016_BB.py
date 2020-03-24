import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from math import sqrt
from  resolution_cfg_2016 import DCB_para
#nBkg = -1
nBkg = -1 
#### to be updated!
dataFile = "input/eventList_ele_2016_BB.txt"

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
	nz   =  5594513 
	nsig_scale = 1./0.084993 # 1./0.895 signal efficiency at z peak       
	eff = signalEff(mass,spin2)
	result = (nsig_scale*nz*eff)

	return result

def signalEff(mass,spin2=False):

	eff_a     = 0.6386
	eff_b     = -497.7
	eff_c     = 348.7
	eff_d     = 6.957e04
	eff_e	  = 1.154e05
	if spin2:
		eff_a = 0.5043
		eff_b = 1173.
		eff_c = 1.1e04
		eff_d = -1.716e05
		eff_e = 6.212e05
	
	return (eff_a+eff_b/(mass+eff_c)+eff_d/(mass*mass+eff_e))
	
def provideUncertaintiesCI(mass):

	result = {}

	result["trig"] = 1.06
	result["zPeak"] = 1.01
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



	

def provideUncertainties(mass):

	result = {}

	result["sigEff"] = [1.06]
	result["massScale"] = 0.02
	result ["bkgUncert"] = 1.4
	result ["res"] = 0.0
	result ["reco"] = [0.0]

 	result["bkgParams"] = {"bkg_a":0.0084073664461949364,"bkg_b":0.0203735214064399046,"bkg_c":0.0,"bkg_d":0.7017137655046845612,"bkg_e":0.0014447387515206115,"bkg_a2":0.0077223413294618785,"bkg_b2":0.0051041020304907732,"bkg_c2":0.0,"bkg_d2":0.0198210948929096360,"bkg_e2":0.0009489989831686973}  

 	# result["bkgParams"] = {"bkg_a":0.009731895898564085,"bkg_b":0.017884574354597824,"bkg_c":0.,"bkg_d":0.37565702379706756,"bkg_e":0.001678522931594343,"bkg_a2":0.009091887838148778,"bkg_b2":0.007738957950179582,"bkg_c2":0.0,"bkg_d2":0.04290977168556734,"bkg_e2":0.0011365419472114415}	
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
	#result['bkgUncert'] = 'dielectron_2016_EBEB' 
	#result['res'] = 'dielectron' 
	#result['bkgParams'] = 'dielectron_2016_EBEB' 
	result['sigEff'] = 'dielectron' 
	result['massScale'] = 'dielectron_2016_EBEB' 
	result['bkgUncert'] = 'dielectron_2016_EBEB' 
	result['res'] = 'dielectron_2016_EBEB' 
	result['reco'] = 'dielectron_2016_EBEB' 
	result['bkgParams'] = 'dielectron_2016_EBEB' 

	return result

def getResolution(mass):

	result = {}
	dcb = DCB_para('dcb')
	dcb.get_value(mass,True) ##True for BB, False for BE

	result['nL'] = dcb.PowerL
	result['nR'] = dcb.PowerR
	result['alphaL'] = dcb.CutL
	result['alphaR'] = dcb.CutR
	result['res'] = dcb.sigma
	result['scale'] = dcb.mean

	return result



def loadBackgroundShape(ws,useShapeUncert):
	# /afs/cern.ch/user/w/wenxing/public/ZprimeRun2/bkg_fit/bkg_fit_newFR.txt
	bkg_a = RooRealVar('bkg_a_dielectron_2016_EBEB','bkg_a_dielectron_2016_EBEB',3.24607e+00)
	bkg_b = RooRealVar('bkg_b_dielectron_2016_EBEB','bkg_b_dielectron_2016_EBEB',-1.75251e-03)
	bkg_c = RooRealVar('bkg_c_dielectron_2016_EBEB','bkg_c_dielectron_2016_EBEB',0)
	bkg_d = RooRealVar('bkg_d_dielectron_2016_EBEB','bkg_d_dielectron_2016_EBEB',-1.31815e-11)
	bkg_e = RooRealVar('bkg_e_dielectron_2016_EBEB','bkg_e_dielectron_2016_EBEB',-3.70739e+00)

	bkg_a2 = RooRealVar('bkg_a2_dielectron_2016_EBEB','bkg_a2_dielectron_2016_EBEB',3.55422e+00)
	bkg_b2 = RooRealVar('bkg_b2_dielectron_2016_EBEB','bkg_b2_dielectron_2016_EBEB',-1.04753e-03)
	bkg_c2 = RooRealVar('bkg_c2_dielectron_2016_EBEB','bkg_c2_dielectron_2016_EBEB',0.)
	bkg_d2 = RooRealVar('bkg_d2_dielectron_2016_EBEB','bkg_d2_dielectron_2016_EBEB',-8.82479e-12)
	bkg_e2 = RooRealVar('bkg_e2_dielectron_2016_EBEB','bkg_e2_dielectron_2016_EBEB',-3.82561e+00)

	# bkg_a = RooRealVar('bkg_a_dielectron_2016_EBEB','bkg_a_dielectron_2016_EBEB',3.28641e+00)
	# bkg_b = RooRealVar('bkg_b_dielectron_2016_EBEB','bkg_b_dielectron_2016_EBEB',-1.68732e-03)
	# bkg_c = RooRealVar('bkg_c_dielectron_2016_EBEB','bkg_c_dielectron_2016_EBEB',0)
	# bkg_d = RooRealVar('bkg_d_dielectron_2016_EBEB','bkg_d_dielectron_2016_EBEB',-2.19817e-11)
	# bkg_e = RooRealVar('bkg_e_dielectron_2016_EBEB','bkg_e_dielectron_2016_EBEB',-3.72630)
	# bkg_a2 = RooRealVar('bkg_a2_dielectron_2016_EBEB','bkg_a2_dielectron_2016_EBEB',3.52027e+00)
	# bkg_b2 = RooRealVar('bkg_b2_dielectron_2016_EBEB','bkg_b2_dielectron_2016_EBEB',-1.03853e-03)
	# bkg_c2 = RooRealVar('bkg_c2_dielectron_2016_EBEB','bkg_c2_dielectron_2016_EBEB',0.)
	# bkg_d2 = RooRealVar('bkg_d2_dielectron_2016_EBEB','bkg_d2_dielectron_2016_EBEB',-9.29420e-12)
	# bkg_e2 = RooRealVar('bkg_e2_dielectron_2016_EBEB','bkg_e2_dielectron_2016_EBEB',-3.82886e+00)
	

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
	bkg_syst_a = RooRealVar('bkg_syst_a_dielectron_2016_EBEB','bkg_syst_a_dielectron_2016_EBEB',1.0)
	bkg_syst_b = RooRealVar('bkg_syst_b_dielectron_2016_EBEB','bkg_syst_b_dielectron_2016_EBEB',0.000)
	#bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',-0.000017)
	bkg_syst_a.setConstant()
	bkg_syst_b.setConstant()
	getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
	getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())
	if useShapeUncert:
		bkgParamsUncert = provideUncertainties(1000)["bkgParams"]
		for uncert in bkgParamsUncert:
			addBkgUncertPrior(ws,uncert,"dielectron_2016_EBEB",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dielectron_2016_EBEB",0.1 )
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2016_EBEB(mass_dielectron_2016_EBEB, bkg_a_dielectron_2016_EBEB_forUse, bkg_b_dielectron_2016_EBEB_forUse, bkg_c_dielectron_2016_EBEB_forUse,bkg_d_dielectron_2016_EBEB_forUse,bkg_e_dielectron_2016_EBEB_forUse,bkg_a2_dielectron_2016_EBEB_forUse, bkg_b2_dielectron_2016_EBEB_forUse, bkg_c2_dielectron_2016_EBEB_forUse,bkg_d2_dielectron_2016_EBEB_forUse,bkg_e2_dielectron_2016_EBEB_forUse,bkg_syst_a_dielectron_2016_EBEB,bkg_syst_b_dielectron_2016_EBEB)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2016_EBEB_forUse, bkg_b_dielectron_2016_EBEB_forUse, bkg_c_dielectron_2016_EBEB_forUse,bkg_d_dielectron_2016_EBEB_forUse,bkg_e_dielectron_2016_EBEB_forUse,bkg_a2_dielectron_2016_EBEB_forUse, bkg_b2_dielectron_2016_EBEB_forUse, bkg_c2_dielectron_2016_EBEB_forUse,bkg_d2_dielectron_2016_EBEB_forUse,bkg_e2_dielectron_2016_EBEB_forUse,bkg_syst_a_dielectron_2016_EBEB,bkg_syst_b_dielectron_2016_EBEB)")		
	else:
	# background shape
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_dielectron_2016_EBEB(mass_dielectron_2016_EBEB, bkg_a_dielectron_2016_EBEB, bkg_b_dielectron_2016_EBEB, bkg_c_dielectron_2016_EBEB,bkg_d_dielectron_2016_EBEB,bkg_e_dielectron_2016_EBEB,bkg_a2_dielectron_2016_EBEB, bkg_b2_dielectron_2016_EBEB, bkg_c2_dielectron_2016_EBEB,bkg_d2_dielectron_2016_EBEB,bkg_e2_dielectron_2016_EBEB,bkg_syst_a_dielectron_2016_EBEB,bkg_syst_b_dielectron_2016_EBEB)")		
		ws.factory("ZPrimeEleBkgPdf3::bkgpdf_fullRange(massFullRange, bkg_a_dielectron_2016_EBEB, bkg_b_dielectron_2016_EBEB, bkg_c_dielectron_2016_EBEB,bkg_d_dielectron_2016_EBEB,bkg_e_dielectron_2016_EBEB,bkg_a2_dielectron_2016_EBEB, bkg_b2_dielectron_2016_EBEB, bkg_c2_dielectron_2016_EBEB,bkg_d2_dielectron_2016_EBEB,bkg_e2_dielectron_2016_EBEB,bkg_syst_a_dielectron_2016_EBEB,bkg_syst_b_dielectron_2016_EBEB)")		
	return ws
