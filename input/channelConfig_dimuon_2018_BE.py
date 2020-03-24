import ROOT,sys
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1 
from ROOT import *
from numpy import exp
nBkg = -1
from muonResolution2018 import getResolution as getRes
dataFile = "input/event_list_2018_beee_clean_sort.txt"
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
	nz   =  38053                      #From Chris
	nsig_scale = 3080.7147258163895       # prescale/eff_z (500/0.1623) -->derives the lumi 
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
			a =  13.4
			b =  6.693
			c = -4.852e+06
			d = -7.437e+06
			e = -81.43
			f = -1.068
			return a - b * exp( -( (mass - c) / d) ) + e * mass**f
		else:
			eff_a     =  0.3154
			eff_b     =  0.04561
			eff_c     =  1.362
			eff_d     = -4927.
			eff_e     =  727.5
			return	eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )

#### flat after 2 TeV
#		if mass <= 450:
#			a =  13.37
#			b = 6.707
#			c = -4.869e+06
#			d = -7.405e+06
#			e = -1476.
#			f = -1.702
#			return a - b * exp( -( (mass - c) / d) ) + e * mass**f
#		else:
#			eff_a     =  0.2174
#			eff_b     =  0.08822
#			eff_c     =  0.7599
#			eff_d     = -4281.
#			eff_e     =  1209.
#			return	eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )
##### linear
#		if mass <= 450:
#			a =  13.38
#			b = 6.707
#			c = -4.869e+06
#			d = -7.403e+06
#			e = -1472.
#			f = -1.701
#			return a - b * exp( -( (mass - c) / d) ) + e * mass**f
#		else:
#			eff_a     =  0.1924
#			eff_b     =  0.09908
#			eff_c     =  0.6725
#			eff_d     =  -4112.
#			eff_e     =  1356.
#			return	eff_a + eff_b * mass**eff_c * exp(- ((mass - eff_d ) / eff_e) )


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

	uncertUp = (uncertHLT**2 + uncertID**2)**0.5
	uncertDown = (uncertHLT**2 + uncertID**2)**0.5
	
	return [1.-uncertDown,1.+uncertUp]

def massScaleUncert(mass):

	scale = 1.00051 -2.21167e-06*mass + 2.21347e-09*mass**2 -7.72178e-13*mass**3 +  1.28101e-16*mass**4  -8.32675e-21*mass**5

	return 1.-scale	



def provideUncertainties(mass):

	result = {}

	result["reco"] = recoUncert(mass)
	result["sigEff"] = signalEffUncert(mass)
	result["massScale"] = massScaleUncert(mass)
	result ["bkgUncert"] = 1.4	
	result ["res"] = 0.085
	result["bkgParams"] = {"bkg_a":0.00148297943775,"bkg_b":0.00598482435414,"bkg_c":0.0359032853638,"bkg_d":0.0379917270419,"bkg_e":0.00126447031957, "bkg_a2":0.00514470341258,"bkg_b2":0.010327332151,"bkg_c2":0.0318691307695,"bkg_e2":0.00703598659419}
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
	#result['bkgUncert'] = 'dimuon_2018_BE' 
	#result['res'] = 'dimuon' 
	#result['bkgParams'] = 'dimuon_2018_BE' 
	result['sigEff'] = 'dimuon' 
	result['massScale'] = 'dimuon' 
	result['bkgUncert'] = 'dimuon_2018_BE' 
	result['res'] = 'dimuon_BE' 
	result['reco'] = 'dimuon_2018_BE' 
	result['bkgParams'] = 'dimuon_2018_BE' 

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
	result['nL'] = params['nL']['BE']  # was BB
	result['nR'] = params['nR']['BE']  # was BB
	result['res'] = params['sigma']['BE']
	result['scale'] = params['scale']['BE']
	return result


def loadBackgroundShape(ws,useShapeUncert=False):


	bkg_a = RooRealVar('bkg_a_dimuon_2018_BE','bkg_a_dimuon_2018_BE', 19.9085672993)
	bkg_b = RooRealVar('bkg_b_dimuon_2018_BE','bkg_b_dimuon_2018_BE',-0.00178754348461)
	bkg_c = RooRealVar('bkg_c_dimuon_2018_BE','bkg_c_dimuon_2018_BE', 1.00128531134e-07)
	bkg_d = RooRealVar('bkg_d_dimuon_2018_BE','bkg_d_dimuon_2018_BE',-1.18860905314e-11)
	bkg_e = RooRealVar('bkg_e_dimuon_2018_BE','bkg_e_dimuon_2018_BE',-3.79943539967)
	bkg_a2 = RooRealVar('bkg_a2_dimuon_2018_BE','bkg_a2_dimuon_2018_BE', 7.14445384949)
	bkg_b2 = RooRealVar('bkg_b2_dimuon_2018_BE','bkg_b2_dimuon_2018_BE',-0.0127392039732)
	bkg_c2 = RooRealVar('bkg_c2_dimuon_2018_BE','bkg_c2_dimuon_2018_BE', 6.33591389821e-06)
	bkg_e2 = RooRealVar('bkg_e2_dimuon_2018_BE','bkg_e2_dimuon_2018_BE',-1.10332904037)

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
			addBkgUncertPrior(ws,uncert,"dimuon_2018_BE",bkgParamsUncert[uncert] )
			#addBkgUncertPrior(ws,uncert,"dimuon_2018_BE",0.10 )

	# background shape
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2018_BE(mass_dimuon_2018_BE, bkg_a_dimuon_2018_BE_forUse, bkg_b_dimuon_2018_BE_forUse, bkg_c_dimuon_2018_BE_forUse,bkg_d_dimuon_2018_BE_forUse,bkg_e_dimuon_2018_BE_forUse, bkg_a2_dimuon_2018_BE_forUse, bkg_b2_dimuon_2018_BE_forUse, bkg_c2_dimuon_2018_BE_forUse,bkg_e2_dimuon_2018_BE_forUse,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2018_BE_forUse, bkg_b_dimuon_2018_BE_forUse, bkg_c_dimuon_2018_BE_forUse,bkg_d_dimuon_2018_BE_forUse,bkg_e_dimuon_2018_BE_forUse,bkg_a2_dimuon_2018_BE_forUse, bkg_b2_dimuon_2018_BE_forUse, bkg_c2_dimuon_2018_BE_forUse,bkg_e2_dimuon_2018_BE_forUse,bkg_syst_a,bkg_syst_b)")	
	else:
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_dimuon_2018_BE(mass_dimuon_2018_BE, bkg_a_dimuon_2018_BE, bkg_b_dimuon_2018_BE, bkg_c_dimuon_2018_BE,bkg_d_dimuon_2018_BE,bkg_e_dimuon_2018_BE, bkg_a2_dimuon_2018_BE, bkg_b2_dimuon_2018_BE, bkg_c2_dimuon_2018_BE,bkg_e2_dimuon_2018_BE,bkg_syst_a,bkg_syst_b)")		
		ws.factory("ZPrimeMuonBkgPdf2::bkgpdf_fullRange(massFullRange, bkg_a_dimuon_2018_BE, bkg_b_dimuon_2018_BE, bkg_c_dimuon_2018_BE,bkg_d_dimuon_2018_BE,bkg_e_dimuon_2018_BE,bkg_a2_dimuon_2018_BE, bkg_b2_dimuon_2018_BE, bkg_c2_dimuon_2018_BE,bkg_e2_dimuon_2018_BE,bkg_syst_a,bkg_syst_b)")	
	
	return ws
