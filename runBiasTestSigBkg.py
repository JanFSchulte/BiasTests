import os, sys
import gc
from array import array
import numpy as np
import math
import glob
from ROOT import gROOT
   



# massRanges = {"default":[140,300,5],  "MRBW":[140,300,5],  "MBWZg":[140,300,5],  "Exp":[0,10,4],  "Gaussian":[81.2,101.2,3]}
# drawRanges = {"default":[150,300,5],  "MRBW":[150,300,5],  "MBWZg":[150,300,5],  "Exp":[0,10,4],  "Gaussian":[81.2,101.2,3]}
# nBins =      {"default":100,          "MRBW":100,          "MBWZg":100,          "Exp":100     ,  "Gaussian":60}

def makeRatioGraph(ws,f1,f2,xMin,xMax):
	graph = ROOT.TGraph()
	ratioList = {}
	valueList = {}
	for i in range(0,100):
		x = xMin + i*(float(xMax)-xMin)/100
		ws.var("mass_%s"%channel).setVal(x)
		sett = ROOT.RooArgSet(ws.var("mass_%s"%channel))
		if f2.getVal(sett) > 0:
			graph.SetPoint(i,x,f1.getVal(sett)/f2.getVal(sett))
			ratioList[x] = f1.getVal(sett)/f2.getVal(sett)
			valueList[x] = f2.getVal(sett)
		else:
			graph.SetPoint(i,x,0)
			ratioList[x] = 0
			valueList[x] = 0
	return graph, ratioList, valueList

def trueValueList(ws,f1,xMin,xMax):
	valueList = {}
	for i in range(0,100):
		x = xMin + i*(float(xMax)-xMin)/100
		ws.var("mass_%s"%channel).setVal(x)
		sett = ROOT.RooArgSet(ws.var("mass_%s"%channel))
		valueList[x] = f1.getVal(sett)

	return valueList

def setIntegrator(ws,name):
	config = RooNumIntConfig(ws.pdf(name).getIntegratorConfig())
	config.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")
	config.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","61Points")
	config.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",1000)
	config.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")
	config.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","61Points")
	config.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",1000)
	ws.pdf(name).setIntegratorConfig(config)

def addGaussianConstraint(ws,name,value,err):
	ws.factory("%sMeasured[%f]" % (name,value))
	ws.factory("%sMeasuredErr[%f]" % (name,err))
	ws.factory("Gaussian::constraint_%s(%s,%sMeasured,%sMeasuredErr)"%(name,name,name,name))

def getTestFunc(ws,func,channel,i,fix):

	if func == "default":
		bkg_a_forFit = RooRealVar('bkg_a_forFit%d_%s'%(i,channel),'bkg_a_forFit%d_%s'%(i,channel),5.820163099,5.820163099 - 5.820163099*0.00462037843829846921,5.820163099 + 5.820163099*0.00462037843829846921)
		bkg_b_forFit = RooRealVar('bkg_b_forFit%d_%s'%(i,channel),'bkg_b_forFit%d_%s'%(i,channel),-0.006025756781,-0.006025756781 - 0.006025756781*0.00767161505186546560,-0.006025756781 + 0.006025756781*0.00767161505186546560)
		bkg_c_forFit = RooRealVar('bkg_c_forFit%d_%s'%(i,channel),'bkg_c_forFit%d_%s'%(i,channel),7.627325877e-06,7.627325877e-06 - 7.627325877e-06*0.02965324582263157913,7.627325877e-06 + 7.627325877e-06*0.02965324582263157913)
		bkg_d_forFit = RooRealVar('bkg_d_forFit%d_%s'%(i,channel),'bkg_d_forFit%d_%s'%(i,channel),0.0)
		bkg_e_forFit = RooRealVar('bkg_e_forFit%d_%s'%(i,channel),'bkg_e_forFit%d_%s'%(i,channel),-1.852928326,-1.852928326 - 1.852928326*0.00212914959129401320, -1.852928326 + 1.852928326*0.00212914959129401320)

		bkg_b2_forFit = RooRealVar('bkg_b2_forFit%d_%s'%(i,channel),'bkg_b2_forFit%d_%s'%(i,channel),-0.0004744949691, -0.0004744949691 - 0.0004744949691*0.03663001163735626203,-0.0004744949691 + 0.0004744949691*0.03663001163735626203)
		bkg_c2_forFit = RooRealVar('bkg_c2_forFit%d_%s'%(i,channel),'bkg_c2_forFit%d_%s'%(i,channel),-1.439508891e-07,-1.439508891e-07 - 1.439508891e-07*0.02916058499703979776,-1.439508891e-07 + 1.439508891e-07*0.02916058499703979776)
		bkg_d2_forFit = RooRealVar('bkg_d2_forFit%d_%s'%(i,channel),'bkg_d2_forFit%d_%s'%(i,channel),4.469376105e-12,4.469376105e-12 - 4.469376105e-12*0.14697581742228424395,4.469376105e-12 + 4.469376105e-12*0.14697581742228424395)
		bkg_thr_forFit = RooRealVar('bkg_thr_forFit%d_%s'%(i,channel),'bkg_thr_forFit%d_%s'%(i,channel),430.6363697,430.6363697 - 430.6363697*0.00799801495493612097,430.6363697 + 430.6363697*0.00799801495493612097)




		if True:
			bkg_a_forFit.setConstant()
			bkg_b_forFit.setConstant()
			bkg_c_forFit.setConstant()
			bkg_d_forFit.setConstant()
			bkg_e_forFit.setConstant()
			bkg_b2_forFit.setConstant()
			bkg_c2_forFit.setConstant()
			bkg_d2_forFit.setConstant()
			bkg_thr_forFit.setConstant()

		getattr(ws,'import')(bkg_a_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_e_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b2_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c2_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d2_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_thr_forFit,ROOT.RooCmdArg())
		addGaussianConstraint(ws,'bkg_a_forFit%d_%s'%(i,channel),5.820163099,5.820163099*0.00462037843829846921)
		addGaussianConstraint(ws,'bkg_b_forFit%d_%s'%(i,channel),-0.006025756781,0.006025756781*0.00767161505186546560)
		addGaussianConstraint(ws,'bkg_c_forFit%d_%s'%(i,channel),7.627325877e-06,7.627325877e-06*0.02965324582263157913)
		addGaussianConstraint(ws,'bkg_e_forFit%d_%s'%(i,channel),-1.852928326,1.852928326*0.00212914959129401320)
		addGaussianConstraint(ws,'bkg_b2_forFit%d_%s'%(i,channel),-0.0004744949691,0.0004744949691*0.03663001163735626203)
		addGaussianConstraint(ws,'bkg_c2_forFit%d_%s'%(i,channel),-1.439508891e-07,1.439508891e-07*0.02916058499703979776)
		addGaussianConstraint(ws,'bkg_d2_forFit%d_%s'%(i,channel),4.469376105e-12,4.469376105e-12*0.14697581742228424395)
		addGaussianConstraint(ws,'bkg_thr_forFit%d_%s'%(i,channel),430.6363697,430.6363697*0.00799801495493612097)



		bkg_syst_a_forFit = RooRealVar('bkg_syst_a_forFit%d_%s'%(i,channel),'bkg_syst_a_forFit%d_%s'%(i,channel),1.0)
		bkg_syst_b_forFit = RooRealVar('bkg_syst_b_forFit%d_%s'%(i,channel),'bkg_syst_b_forFit%d_%s'%(i,channel),0.0)
		bkg_syst_a_forFit.setConstant()
		bkg_syst_b_forFit.setConstant()
		getattr(ws,'import')(bkg_syst_a_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_syst_b_forFit,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_%s_forFit%d(mass_%s, bkg_a_forFit%d_%s, bkg_b_forFit%d_%s, bkg_c_forFit%d_%s, bkg_d_forFit%d_%s, bkg_e_forFit%d_%s, bkg_b2_forFit%d_%s, bkg_c2_forFit%d_%s, bkg_d2_forFit%d_%s, bkg_thr_forFit%d_%s, bkg_syst_a_forFit%d_%s, bkg_syst_b_forFit%d_%s)" % (channel,i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel, i,channel))


	elif func == "defaultlow":
		bkg_a_forFit = RooRealVar('bkg_a_forFit%d_%s'%(i,channel),'bkg_a_forFit%d_%s'%(i,channel),6,3,10)
		bkg_b_forFit = RooRealVar('bkg_b_forFit%d_%s'%(i,channel),'bkg_b_forFit%d_%s'%(i,channel),-0.006025756781,-1,0)
		bkg_c_forFit = RooRealVar('bkg_c_forFit%d_%s'%(i,channel),'bkg_c_forFit%d_%s'%(i,channel),7.627325877e-06,0,1e-02)
		bkg_d_forFit = RooRealVar('bkg_d_forFit%d_%s'%(i,channel),'bkg_d_forFit%d_%s'%(i,channel),0.0)
		bkg_e_forFit = RooRealVar('bkg_e_forFit%d_%s'%(i,channel),'bkg_e_forFit%d_%s'%(i,channel),-1.852928326,-20,0)
		bkg_thr_forFit = RooRealVar('bkg_thr_forFit%d_%s'%(i,channel),'bkg_thr_forFit%d_%s'%(i,channel),430.6363697)

		if fix:
			bkg_a_forFit.setConstant()
			bkg_b_forFit.setConstant()
			bkg_c_forFit.setConstant()
			bkg_e_forFit.setConstant()

		bkg_d_forFit.setConstant()
		bkg_thr_forFit.setConstant()
		getattr(ws,'import')(bkg_a_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_e_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_thr_forFit,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdf4Low::bkgpdf_%s_forFit%d(mass_%s, bkg_a_forFit%d_%s, bkg_b_forFit%d_%s, bkg_c_forFit%d_%s,bkg_d_forFit%d_%s,bkg_e_forFit%d_%s,bkg_thr_forFit%d_%s)"%(channel,i,channel,i,channel,i,channel,i,channel,i,channel,i,channel,i,channel))

	elif func == "MBWZg":
		bkg_a_forFit = RooRealVar('bkg_a_forFit%d_%s'%(i,channel),'bkg_a_forFit%d_%s'%(i,channel),1)
		bkg_b_forFit = RooRealVar('bkg_b_forFit%d_%s'%(i,channel),'bkg_b_forFit%d_%s'%(i,channel),-0.001864061214,-1e-01,1)
		bkg_c_forFit = RooRealVar('bkg_c_forFit%d_%s'%(i,channel),'bkg_c_forFit%d_%s'%(i,channel),-1.480749517e-05,-1,1)
		bkg_d_forFit = RooRealVar('bkg_d_forFit%d_%s'%(i,channel),'bkg_d_forFit%d_%s'%(i,channel),3.484180767,-10,10)

		if fix:
			bkg_a_forFit.setConstant()
			bkg_b_forFit.setConstant()
			bkg_c_forFit.setConstant()
			bkg_d_forFit.setConstant()

		getattr(ws,'import')(bkg_a_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d_forFit,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdfMBWZg::bkgpdf_%s_forFit%d(mass_%s, bkg_a_forFit%d_%s, bkg_b_forFit%d_%s, bkg_c_forFit%d_%s, bkg_d_forFit%d_%s)"%(channel,i,channel,i,channel,i,channel,i,channel,i,channel))

	elif func == "MRBW":
		bkg_a_forFit = RooRealVar('bkg_a_forFit%d_%s'%(i,channel),'bkg_a_forFit%d_%s'%(i,channel),1.14936683,0,10)
		#bkg_a_forFit = RooRealVar('bkg_a_forFit%d_%s'%(i,channel),'bkg_a_forFit%d_%s'%(i,channel),1)
		bkg_b_forFit = RooRealVar('bkg_b_forFit%d_%s'%(i,channel),'bkg_b_forFit%d_%s'%(i,channel),1.039565181,-10,10)
		bkg_c_forFit = RooRealVar('bkg_c_forFit%d_%s'%(i,channel),'bkg_c_forFit%d_%s'%(i,channel),-1.993409218,-10,10)
		bkg_d_forFit = RooRealVar('bkg_d_forFit%d_%s'%(i,channel),'bkg_d_forFit%d_%s'%(i,channel),0.113374161,-5,5)

		if fix:
			bkg_a_forFit.setConstant()
			bkg_b_forFit.setConstant()
			bkg_c_forFit.setConstant()
			bkg_d_forFit.setConstant()

		getattr(ws,'import')(bkg_a_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c_forFit,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d_forFit,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdfMRBW::bkgpdf_%s_forFit%d(mass_%s, bkg_a_forFit%d_%s, bkg_b_forFit%d_%s, bkg_c_forFit%d_%s, bkg_d_forFit%d_%s)"%(channel,i,channel,i,channel,i,channel,i,channel,i,channel))

	elif func == "Gaussian":
		ws.factory("RooGaussian::bkgpdf_%s_forFit%d(mass_%s, mean_forFit%d[90,80,100],sigma_forFit%d[2.5,1,4])"%(channel,i,channel,i,i))

	elif func == "Exp":
		ws.factory("RooExponential::bkgpdf_%s_forFit%d(mass_%s, const_forFit%d[-0.01,-0.1,0.0])"%(channel,i,channel,i))

def getTruthFunc(ws,func,channel):

	if func == "Gaussian":
		ws.factory("RooGaussian::bkgpdf_%s(mass_%s,mean[90],sigma[2.5])"%(channel,channel))

	elif func == "Exp":
		ws.factory("RooExponential::bkgpdf_%s(mass_%s,const[-0.01])"%(channel,channel))

	elif func == "default":
		bkg_a = RooRealVar('bkg_a_%s'%channel,'bkg_a_%s'%channel,5.820163099)
		bkg_b = RooRealVar('bkg_b_%s'%channel,'bkg_b_%s'%channel,-0.006025756781)
		bkg_c = RooRealVar('bkg_c_%s'%channel,'bkg_c_%s'%channel,7.627325877e-06)
		bkg_d = RooRealVar('bkg_d_%s'%channel,'bkg_d_%s'%channel,0.0)
		bkg_e = RooRealVar('bkg_e_%s'%channel,'bkg_e_%s'%channel,-1.852928326)

		bkg_b2 = RooRealVar('bkg_b2_%s'%channel,'bkg_b2_%s'%channel,-0.0004744949691)
		bkg_c2 = RooRealVar('bkg_c2_%s'%channel,'bkg_c2_%s'%channel,-1.439508891e-07)
		bkg_d2 = RooRealVar('bkg_d2_%s'%channel,'bkg_d2_%s'%channel,4.469376105e-12)
		bkg_thr = RooRealVar('bkg_thr_%s'%channel,'bkg_thr_%s'%channel,430.6363697)

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

		bkg_syst_a = RooRealVar('bkg_syst_a','bkg_syst_a',1.0)
		bkg_syst_b = RooRealVar('bkg_syst_b','bkg_syst_b',0.0)
		bkg_syst_a.setConstant()
		bkg_syst_b.setConstant()
		getattr(ws,'import')(bkg_syst_a,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_syst_b,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdf4::bkgpdf_%s(mass_%s, bkg_a_%s, bkg_b_%s, bkg_c_%s,bkg_d_%s,bkg_e_%s,bkg_b2_%s, bkg_c2_%s,bkg_d2_%s, bkg_thr_%s ,bkg_syst_a,bkg_syst_b)"%(channel,channel,channel,channel,channel,channel,channel,channel,channel,channel,channel))

	elif func == "MBWZg":
		bkg_a = RooRealVar('bkg_a_%s'%(channel),'bkg_a_%s'%(channel),1)
		bkg_b = RooRealVar('bkg_b_%s'%(channel),'bkg_b_%s'%(channel),-0.001864061214,-1e-01,1)
		bkg_c = RooRealVar('bkg_c_%s'%(channel),'bkg_c_%s'%(channel),-1.480749517e-05,-1,1)
		bkg_d = RooRealVar('bkg_d_%s'%(channel),'bkg_d_%s'%(channel),3.484180767,-10,10)

		bkg_a.setConstant()
		bkg_b.setConstant()
		bkg_c.setConstant()
		bkg_d.setConstant()

		getattr(ws,'import')(bkg_a,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdfMBWZg::bkgpdf_%s(mass_%s, bkg_a_%s, bkg_b_%s, bkg_c_%s, bkg_d_%s)"%(channel,channel,channel,channel,channel,channel))

	elif func == "MRBW":
		bkg_a = RooRealVar('bkg_a_%s'%(channel),'bkg_a_%s'%(channel),1.14936683,0,10)
		#bkg_a = RooRealVar('bkg_a_%s'%(channel),'bkg_a_%s'%(channel),1)
		bkg_b = RooRealVar('bkg_b_%s'%(channel),'bkg_b_%s'%(channel),1.039565181,-10,10)
		bkg_c = RooRealVar('bkg_c_%s'%(channel),'bkg_c_%s'%(channel),-1.993409218,-10,10)
		bkg_d = RooRealVar('bkg_d_%s'%(channel),'bkg_d_%s'%(channel),0.113374161,-5,5)

		bkg_a.setConstant()
		bkg_b.setConstant()
		bkg_c.setConstant()
		bkg_d.setConstant()

		getattr(ws,'import')(bkg_a,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_b,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_c,ROOT.RooCmdArg())
		getattr(ws,'import')(bkg_d,ROOT.RooCmdArg())

		ws.factory("ZPrimeMuonBkgPdfMRBW::bkgpdf_%s(mass_%s, bkg_a_%s, bkg_b_%s, bkg_c_%s, bkg_d_%s)"%(channel,channel,channel,channel,channel,channel))

	# else:  # func == "defaultlow":
	# 	bkg_a = RooRealVar('bkg_a_%s'%channel,'bkg_a_%s'%channel,5.820163099)
	# 	bkg_b = RooRealVar('bkg_b_%s'%channel,'bkg_b_%s'%channel,-0.006025756781)
	# 	bkg_c = RooRealVar('bkg_c_%s'%channel,'bkg_c_%s'%channel,7.627325877e-06)
	# 	bkg_d = RooRealVar('bkg_d_%s'%channel,'bkg_d_%s'%channel,0.0)
	# 	bkg_e = RooRealVar('bkg_e_%s'%channel,'bkg_e_%s'%channel,-1.852928326)
	# 	bkg_thr = RooRealVar('bkg_thr_%s'%channel,'bkg_thr_%s'%channel,430.6363697)

	# 	bkg_a.setConstant()
	# 	bkg_b.setConstant()
	# 	bkg_c.setConstant()
	# 	bkg_d.setConstant()
	# 	bkg_e.setConstant()
	# 	bkg_thr.setConstant()
	# 	getattr(ws,'import')(bkg_a,ROOT.RooCmdArg())
	# 	getattr(ws,'import')(bkg_b,ROOT.RooCmdArg())
	# 	getattr(ws,'import')(bkg_c,ROOT.RooCmdArg())
	# 	getattr(ws,'import')(bkg_d,ROOT.RooCmdArg())
	# 	getattr(ws,'import')(bkg_e,ROOT.RooCmdArg())
	# 	getattr(ws,'import')(bkg_thr,ROOT.RooCmdArg())

	# 	ws.factory("ZPrimeMuonBkgPdf4Low::bkgpdf_%s(mass_%s, bkg_a_%s, bkg_b_%s, bkg_c_%s,bkg_d_%s,bkg_e_%s,bkg_thr_%s)"%(channel,channel,channel,channel,channel,channel,channel,channel))

	setIntegrator(ws,"bkgpdf_%s"%channel)

def getSigFunc(ws,massVal,width,channel,config,i):

	params = config.getResolution(massVal)
	effWidth = width + params['res']

	peakName = "%d" % i
	peak = RooRealVar("peak%s"%peakName,"peak%s"%peakName,massVal, massVal, massVal)
	peak.setConstant()
	getattr(ws,'import')(peak,ROOT.RooCmdArg())

	beta_peak = RooRealVar('beta_peak%s'%peakName,'beta_peak%s'%peakName,0,-5,5)
	getattr(ws,'import')(beta_peak,ROOT.RooCmdArg())
	scaleUncert = 1. + config.provideUncertainties(massVal)["massScale"]
	peak_kappa = RooRealVar('peak%s_kappa'%peakName,'peak%s_kappa'%peakName,scaleUncert)
	peak_kappa.setConstant()
	getattr(ws,'import')(peak_kappa,ROOT.RooCmdArg())
	ws.factory("PowFunc::peak_nuis%s(peak%s_kappa, beta_peak%s)"%(peakName,peakName,peakName))
	ws.factory("prod::peak_scaled%s(peak%s, peak_nuis%s)"%(peakName,peakName,peakName))

	res = RooRealVar("res%s"%peakName,"res%s"%peakName, massVal*params['res'])
	res.setConstant()
	getattr(ws,'import')(res,ROOT.RooCmdArg())

	mean = RooRealVar("mean%s"%peakName,"mean%s"%peakName, params['scale'])
	mean.setConstant()
	getattr(ws,'import')(mean,ROOT.RooCmdArg())

	alphaL = RooRealVar("alphaL_%s_%s"%(channel,peakName),"alphaL_%s_%s"%(channel,peakName),params['alphaL'])
	alphaL.setConstant()
	getattr(ws,'import')(alphaL,ROOT.RooCmdArg())

	nL = RooRealVar("nL_%s_%s"%(channel,peakName),"nL_%s_%s"%(channel,peakName),params['nL'])
	nL.setConstant()
	getattr(ws,'import')(nL,ROOT.RooCmdArg())

	alphaR = RooRealVar("alphaR_%s_%s"%(channel,peakName),"alphaR_%s_%s"%(channel,peakName),params['alphaR'])
	alphaR.setConstant()
	getattr(ws,'import')(alphaR,ROOT.RooCmdArg())
	if params['nR'] < 0:
		params['nR'] = 0
	nR = RooRealVar("nR_%s_%s"%(channel,peakName),"nR_%s_%s"%(channel,peakName),params['nR'])
	nR.setConstant()
	getattr(ws,'import')(nR,ROOT.RooCmdArg())

	beta_res = RooRealVar('beta_res%s'%peakName,'beta_res%s'%peakName,0,-5,5)
	getattr(ws,'import')(beta_res,ROOT.RooCmdArg())
	resUncert = 1. + config.provideUncertainties(massVal)["res"]
	res_kappa = RooRealVar('res%s_kappa'%peakName,'res%s_kappa'%peakName,resUncert)
	res_kappa.setConstant()
	getattr(ws,'import')(res_kappa,ROOT.RooCmdArg())
	ws.factory("PowFunc::res_nuis%s(res%s_kappa, beta_res%s)"%(peakName,peakName,peakName))
	ws.factory("prod::res_scaled%s(res%s, res_nuis%s)"%(peakName,peakName,peakName))

	addGaussianConstraint(ws,'beta_res%s'%peakName,0,1)
	addGaussianConstraint(ws,'beta_peak%s'%peakName,0,1)

	ws.var('beta_res%s'%peakName).setConstant()
	ws.var('beta_peak%s'%peakName).setConstant()

	ws.factory("BreitWigner::bw_%s_%s(mass_%s, peak_scaled%s, %.3f)"%(channel,peakName,channel,peakName,massVal*width))
	ws.factory("RooDCBShape::cb_%s_%s(mass_%s, mean_cb[0.0], res_scaled%s, alphaL_%s_%s, alphaR_%s_%s, nL_%s_%s, nR_%s_%s)"%(channel,peakName,channel,peakName,channel,peakName,channel,peakName,channel,peakName,channel,peakName))
	bw = ws.pdf("bw_%s_%s"%(channel,peakName))
	cb = ws.pdf("cb_%s_%s"%(channel,peakName))
	ws.var("mass_%s"%channel).setBins(20000,"cache")
	ws.var("mass_%s"%channel).setMin("cache",0)
	ws.var("mass_%s"%channel).setMax("cache",12500); ## need to be adjusted to be higher than limit setting

	sigpdf = ROOT.RooFFTConvPdf("sig_pdf_%s_%s"%(channel,peakName),"sig_pdf_%s_%s"%(channel,peakName),ws.var("mass_%s"%channel),bw,cb)
	getattr(ws,'import')(sigpdf,ROOT.RooCmdArg())
	setIntegrator(ws,"sig_pdf_%s_%s"%(channel,peakName))

def getSigBkgModel(ws,channel,nEvents,i):
	# nsig = RooRealVar('nsig%d'%i,'nsig%d'%i,0,0,nEvents)
	# getattr(ws,'import')(nsig,ROOT.RooCmdArg())
	# nbkg = RooRealVar('nbkg%d'%i,'nbkg%d'%i,nEvents,0,2*nEvents)
	# getattr(ws,'import')(nbkg,ROOT.RooCmdArg())
	# ws.factory("SUM::model_pdf_%s_forFit%d((nsig%d / (nsig%d+nbkg%d))*sig_pdf_%s_%i, (nbkg%d / (nsig%d+nbkg%d))*bkgpdf_%s_forFit%d)"%(channel,i,i,i,i,channel,i, i,i,i,channel,i))
	ws.factory("SUM::model_pdf_%s_forFit%d(nsig%d[0,%d,%d]*sig_pdf_%s_%i,nbkg%d[%d,0,%d]*bkgpdf_%s_forFit%d)"%(channel,i,i,-1*nEvents,nEvents,channel,i,i,nEvents,2*nEvents,channel,i))

def getSigBkgModelFixed(ws,channel,nEvents,nSignal,i):
	# nsig = RooRealVar('nsig%d'%i,'nsig%d'%i,0,0,nEvents)
	# getattr(ws,'import')(nsig,ROOT.RooCmdArg())
	# nbkg = RooRealVar('nbkg%d'%i,'nbkg%d'%i,nEvents,0,2*nEvents)
	# getattr(ws,'import')(nbkg,ROOT.RooCmdArg())
	# ws.factory("SUM::model_pdf_%s_forFit%d((nsig%d / (nsig%d+nbkg%d))*sig_pdf_%s_%i, (nbkg%d / (nsig%d+nbkg%d))*bkgpdf_%s_forFit%d)"%(channel,i,i,i,i,channel,i, i,i,i,channel,i))
	ws.factory("SUM::model_pdf_%s_fixed%d(nsig_fixed%d[%d,%d,%d]*sig_pdf_%s_%i,nbkg_fixed%d[%d,%d,%d]*bkgpdf_%s)"%(channel,i,i,nSignal,nSignal,nSignal,channel,i,i,nEvents,nEvents,nEvents,channel))

def makePullHist(ws, hist, curve, norm):

	ymax = 0.

	out = ROOT.RooHist(hist.getFitRangeBinW())
	xstart = ROOT.Double(-999)
	xstop  = ROOT.Double(-999)
	y      = ROOT.Double(-999)
	curve.GetPoint(0,xstart,y) ;
	curve.GetPoint(curve.GetN()-1,xstop,y)

	for i in range(hist.GetN()):
		x  = ROOT.Double(-999)
		y  = ROOT.Double(-999)
		hist.GetPoint(i,x,y)

		if (x<xstart) or (x>xstop):
			continue

		yy = ROOT.Double(-999)
		eyl = hist.GetErrorYlow(i)
		eyh = hist.GetErrorYhigh(i)
		exl = hist.GetErrorXlow(i)
		exh = hist.GetErrorXhigh(i)
		if (exl<=0 ):  exl = hist.GetErrorX(i)
		if (exh<=0 ):  exh = hist.GetErrorX(i)
		if (exl<=0 ):  exl = 0.5*hist.getNominalBinWidth()
		if (exh<=0 ):  exh = 0.5*hist.getNominalBinWidth()

		ws.var("mass_%s"%channel).setRange(x-exl, x+exh)
		bkgpdf_integ_i = ws.pdf("bkgpdf_%s"%channel).createIntegral(
			ROOT.RooArgSet(ws.var("mass_%s"%channel)),
			ROOT.RooFit.NormSet(ROOT.RooArgSet(ws.var("mass_%s"%channel)))
		)
		bkgpdf_integ_i = bkgpdf_integ_i.getVal() * norm
		yy = y - bkgpdf_integ_i

		err = eyl if yy>0 else eyh
		if err == 0.:
			yy=0
			eyh=0
			eyl=0
		else:
			yy   = yy / err
			eyh  = eyh / err
			eyl  = eyl / err

		if abs(yy) > ymax:
			ymax = abs(yy)

		out.addBinWithError(x,yy,eyl,eyh)

	return out, ymax

def makeRatioHist(ws, hist, curve, norm):

	ymax = 0.

	out = ROOT.RooHist(hist.getFitRangeBinW())
	xstart = ROOT.Double(-999)
	xstop  = ROOT.Double(-999)
	y      = ROOT.Double(-999)
	curve.GetPoint(0,xstart,y) ;
	curve.GetPoint(curve.GetN()-1,xstop,y)

	for i in range(hist.GetN()):
		x  = ROOT.Double(-999)
		y  = ROOT.Double(-999)
		hist.GetPoint(i,x,y)

		if (x<xstart) or (x>xstop):
			out.addBinWithError(x,0,0,0)
			continue

		yy = ROOT.Double(-999)
		eyl = hist.GetErrorYlow(i)
		eyh = hist.GetErrorYhigh(i)
		exl = hist.GetErrorXlow(i)
		exh = hist.GetErrorXhigh(i)
		if (exl<=0 ):  exl = hist.GetErrorX(i)
		if (exh<=0 ):  exh = hist.GetErrorX(i)
		if (exl<=0 ):  exl = 0.5*hist.getNominalBinWidth()
		if (exh<=0 ):  exh = 0.5*hist.getNominalBinWidth()

		ws.var("mass_%s"%channel).setRange(x-exl, x+exh)
		bkgpdf_integ_i = ws.pdf("bkgpdf_%s"%channel).createIntegral(
			ROOT.RooArgSet(ws.var("mass_%s"%channel)),
			ROOT.RooFit.NormSet(ROOT.RooArgSet(ws.var("mass_%s"%channel)))
		)
		bkgpdf_integ_i = bkgpdf_integ_i.getVal() * norm

		yy = (float(y) / bkgpdf_integ_i) - 1.
		eyh = eyh / bkgpdf_integ_i
		eyl = eyl / bkgpdf_integ_i

		if abs(yy) > ymax:
			ymax = abs(yy)

		out.addBinWithError(x,yy,eyl,eyh)

	return out, ymax


if __name__ == '__main__':
	import ROOT
	from ROOT import RooFit
	from ROOT import RooRealVar, RooMCStudy, RooNumIntConfig
	from ROOT import TF1, TRandom3
	ROOT.gROOT.SetBatch(True)
	ROOT.TH1.AddDirectory(False)
	ROOT.TH1.SetDefaultSumw2(True)
	ROOT.RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)  # WARNING, ERROR, FATAL

	import argparse
	parser = argparse.ArgumentParser(description='Steering tool for Zprime -> ll analysis interpretation in combine')
	parser.add_argument("-ch","--channel",   dest = "channel", default='XXX',     type=str,   help="channel")
	parser.add_argument("-f", "--function",  dest = "func",    default='default', type=str,   help="function to be used for bias test")
	parser.add_argument(      "--nToys",     dest = "nToys",   default=200,       type=int,   help="")
	parser.add_argument(      "--nEvents",   dest = "nEvents", default=-1,        type=int,   help="")
	parser.add_argument(      "--nSigFac",   dest = "nSigFac", default=3,         type=int,   help="")
	parser.add_argument(      "--mass",      dest = "mass",    default= 200.,     type=float, help="")
	parser.add_argument(      "--jobId",     dest = "jobId",   default=0,         type=int,   help="")
	parser.add_argument(      "--fix",       dest = "fix",     default=1,         type=int, help="")
	parser.add_argument(      "--saveAll",   dest = "saveAll", action='store_true', default=False, help="")
	parser.add_argument(      "--binw",      dest = "binw",    default= 10.,      type=float, help="")
	parser.add_argument("-m", "--min",       dest = "min",     default=-1.,       type=float, help="")
	parser.add_argument("-M", "--max",       dest = "max",     default=-1.,       type=float, help="")

	args = parser.parse_args()

	# -- Setup -- #
	if "msoh" in os.popen('whoami').read():
		combDir  = "/u/user/msoh/ZprimeBias/CMSSW_10_2_13/src/ZprimeBias/Analysis/"
	else:
		combDir  = ""
	inputDir = combDir+"input/"
	funcDir  = combDir+"userfuncs/"
	#sys.path.append(combDir)
	sys.path.append(inputDir)
	sys.path.append(funcDir)
        for f in glob.glob("userfuncs/*.cxx"):
                gROOT.ProcessLine(".L "+f+"+")


	ROOT.gSystem.Load(funcDir+"PowFunc_cxx.so")
	ROOT.gSystem.Load(funcDir+"RooDCBShape_cxx.so")
	ROOT.gSystem.Load(funcDir+"ZPrimeMuonBkgPdf4_cxx.so")
	ROOT.gSystem.Load(funcDir+"ZPrimeMuonBkgPdf4Low_cxx.so")
	ROOT.gSystem.Load(funcDir+"ZPrimeMuonBkgPdfMBWZg_cxx.so")
	ROOT.gSystem.Load(funcDir+"ZPrimeMuonBkgPdfMRBW_cxx.so")

	ROOT.gSystem.ListLibraries()

	from tools import getMassRange

	channel     = args.channel
	channelName ="channelConfig_%s"%channel
	channelConfig =  __import__(channelName)
	dataFile = combDir+channelConfig.dataFile

	sigWidth = 0.1

	nToys       = args.nToys
	massVal     = args.mass
	nEvents = 0
	nEventsFix = 0
	if args.nEvents > 0 and args.min > 0 and args.max > 0:
		nEvents     = args.nEvents
		fitRange    = [args.min, args.max]
		binw        = args.binw
		nBinsForFit = int((args.max-args.min)/binw)
	else:
		gamma    = sigWidth
		sigma    = (channelConfig.getResolution(massVal)['res'])
		effWidth = gamma + sigma
		massLow, massHigh = getMassRange(massVal,100,effWidth,dataFile,150,4)
		fitRange = [massLow, massHigh]
		nBinsForFit = 40
		binw = (fitRange[1]-fitRange[0])/float(nBinsForFit)
		nEventsFix = 0
		with open(dataFile) as f:
			_masses = f.readlines()
		for m in _masses:
			if float(m) >= fitRange[0] and float(m) <= fitRange[1]:
				nEventsFix += 1
	drawRanges  = fitRange
	fix = args.fix

	print "Setting: ", args.func, nToys, nEvents, nEventsFix, massVal, fitRange

	canvdir = "canvases_%s_fix%d_nsigfac%d/" % (args.func, fix, args.nSigFac)
	if not os.path.isdir(canvdir):
		os.makedirs(canvdir)

	ws = ROOT.RooWorkspace("w")

	mass = RooRealVar('mass_%s'%channel,'mass_%s'%channel,(fitRange[1]-fitRange[0])/2 ,fitRange[0],fitRange[1])
	mass.setBins(nBinsForFit)
	getattr(ws,'import')(mass,ROOT.RooCmdArg())

	getTruthFunc(ws, func = args.func, channel=channel)

	bkgpdf_integ = ws.pdf("bkgpdf_%s"%channel).createIntegral(
		ROOT.RooArgSet(ws.var("mass_%s"%channel)),
		ROOT.RooFit.NormSet(ROOT.RooArgSet(ws.var("mass_%s"%channel)))
	)
	bkgpdf_integ = bkgpdf_integ.getVal()

	ratios = []
	values = []
	trueValues = trueValueList(ws,ws.pdf("bkgpdf_%s"%channel),fitRange[0],fitRange[1])

	t = ROOT.TTree( 't', 't' )

	nEv       = array( 'i', [-1] )
	nTo       = array( 'i', [-1] )
	t.Branch( 'nEv', nEv, 'nEv/I' )
	t.Branch( 'nTo', nTo, 'nTo/I' )

	isSucD    = array( 'i', [-1] )
	nSigsD    = array( 'f', [ -999. ] )
	eSigsD    = array( 'f', [ -999. ] )
	elSigsD   = array( 'f', [ -999. ] )
	ehSigsD   = array( 'f', [ -999. ] )
	nRelSigsD = array( 'f', [ -999. ] )
	t.Branch( 'isSucD', isSucD, 'isSucD/I' )
	t.Branch( 'nSigsD', nSigsD, 'nSigsD/F' )
	t.Branch( 'eSigsD', eSigsD, 'eSigsD/F' )
	t.Branch( 'elSigsD', elSigsD, 'elSigsD/F' )
	t.Branch( 'ehSigsD', ehSigsD, 'ehSigsD/F' )
	t.Branch( 'nRelSigsD', nRelSigsD, 'nRelSigsD/F' )

	isSucT    = array( 'i', [ -1 ] )
	nSigsT    = array( 'f', [ -999. ] )
	eSigsT    = array( 'f', [ -999. ] )
	elSigsT   = array( 'f', [ -999. ] )
	ehSigsT   = array( 'f', [ -999. ] )
	nRelSigsT = array( 'f', [ -999. ] )
	t.Branch( 'isSucT', isSucT, 'isSucT/I' )
	t.Branch( 'nSigsT', nSigsT, 'nSigsT/F' )
	t.Branch( 'eSigsT', eSigsT, 'eSigsT/F' )
	t.Branch( 'elSigsT', elSigsT, 'elSigsT/F' )
	t.Branch( 'ehSigsT', ehSigsT, 'ehSigsT/F' )
	t.Branch( 'nRelSigsT', nRelSigsT, 'nRelSigsT/F' )

	ROOT.gRandom.SetSeed(0)
	ROOT.RooRandom.randomGenerator().SetSeed(0)
	for i in range(1,nToys+1):

		nEvents = nEventsFix  # ROOT.gRandom.Poisson(nEventsFix)
		nSignal = int(math.sqrt(nEvents)*args.nSigFac)
		print "\t update nEvents for %dth toy from %d to %d. nSignal = %d, (seed=%d)" % (i, nEventsFix, nEvents, nSignal, ROOT.gRandom.GetSeed())
		print ""

		# -- PDFs -- #
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])
		getSigFunc(ws,massVal,sigWidth,channel,channelConfig,-1*i)
		getTestFunc(ws,"default",channel,-1*i,fix)  # fix
		getSigBkgModel(ws,channel,nEvents,-1*i)

		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])
		getSigFunc(ws,massVal,sigWidth,channel,channelConfig,i)
		getTestFunc(ws,args.func,channel,i,fix)
		getSigBkgModel(ws,channel,nEvents,i)
		getSigBkgModelFixed(ws,channel,nEvents,nSignal,i)

		# -- Generate  -- #
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])
		dataSet = ws.pdf("bkgpdf_%s"%channel).generate(ROOT.RooArgSet(ws.var("mass_%s"%channel)),nEvents)
		dataSet.append( ws.pdf("sig_pdf_%s_%d"%(channel,i)).generate(ROOT.RooArgSet(ws.var("mass_%s"%channel)), nSignal))
		dataSet.SetNameTitle("dataSet%d"%i, "dataSet%d"%i)

		dataHist = ROOT.RooDataHist("dataHist%d"%i, "dataHist%d"%i, ROOT.RooArgSet(ws.var("mass_%s"%channel)), dataSet)
		hist = dataHist.createHistogram("hist%d"%i,ws.var("mass_%s"%channel))
		sys.stdout.flush()
		gc.collect()

		# -- Fit  -- #
		nIter = 5
		ws.var("nsig%d"%i).setVal(nSignal)
		ws.var("nsig%d"%(-1*i)).setVal(nSignal)

		bkgConstraints = ROOT.RooArgSet(ws.pdf("constraint_bkg_a_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_b_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_c_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_e_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_b2_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_c2_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_d2_forFit%d_%s"%(-1*i,channel)),ws.pdf("constraint_bkg_thr_forFit%d_%s"%(-1*i,channel)))
		systConstraints = ROOT.RooArgSet(ws.pdf("constraint_beta_peak%d"%(-1*i)),ws.pdf("constraint_beta_res%d"%(-1*i)))
		#fitResultsD = ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(bkgConstraints,systConstraints)),ROOT.RooFit.Save(True))
		fitResultsD = ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(bkgConstraints,systConstraints)),ROOT.RooFit.Save(True))
		fitResultsT = ws.pdf("model_pdf_%s_forFit%d"%(channel,i) ).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(ws.pdf("constraint_beta_peak%d"%i),ws.pdf("constraint_beta_res%d"%i))),ROOT.RooFit.Save(True))
#		fitResultsD = ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.Save(True))
#		fitResultsT = ws.pdf("model_pdf_%s_forFit%d"%(channel,i) ).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.Save(True))

		for it in range(nIter):
			fitResultsD = ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(bkgConstraints,systConstraints)),ROOT.RooFit.Save(True))
#			fitResultsD = ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.Save(True))
			fitResultsT = ws.pdf("model_pdf_%s_forFit%d"%(channel,i) ).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(ws.pdf("constraint_beta_peak%d"%i),ws.pdf("constraint_beta_res%d"%i))),ROOT.RooFit.Save(True))


#			fitResultsT = ws.pdf("model_pdf_%s_forFit%d"%(channel,i) ).fitTo(dataSet,ROOT.RooFit.PrintLevel(0),ROOT.RooFit.Save(True))
		sys.stdout.flush()
		gc.collect()

		# -- Ratio -- #
		ratioD, ratioListD, valueListD = makeRatioGraph(ws,ws.pdf("bkgpdf_%s"%channel),ws.pdf("bkgpdf_%s_forFit%d"%(channel,-1*i)),fitRange[0],fitRange[1])
		ratioD.SetLineColor(ROOT.kRed)
		isFailD = False
		for key, rr in ratioListD.iteritems():
			if rr < 0.01 or rr > 2:
				print "Fit D failed: ", key, rr
				isFailD = True
				break

		ratio, ratioList, valueList    = makeRatioGraph(ws,ws.pdf("bkgpdf_%s"%channel),ws.pdf("bkgpdf_%s_forFit%d"%(channel,i)),fitRange[0],fitRange[1])
		isFail = False
		for key, rr in ratioList.iteritems():
			if rr < 0.01 or rr > 2:
				print "Fit failed: ", key, rr
				isFail = True
				break
		if isFail or isFailD:
			ratios.append([])
			values.append([])
		else:
			ratios.append(ratioList)
			values.append(valueList)

		# init
		nEv[0]       = -1
		nTo[0]       = -1
		isSucD[0]    = -1
		nSigsD[0]    = -99999.
		eSigsD[0]    = -99999.
		elSigsD[0]   = -99999.
		ehSigsD[0]   = -99999.
		nRelSigsD[0] = -99999.
		isSucT[0]    = -1
		nSigsT[0]    = -99999.
		eSigsT[0]    = -99999.
		elSigsT[0]   = -99999.
		ehSigsT[0]   = -99999.
		nRelSigsT[0] = -99999.

		nEv[0]       = nEvents
		nTo[0]       = (int(args.jobId)*nToys + i)
		isSucD[0]    = int(not isFailD)
		nSigsD[0]    = float(fitResultsD.floatParsFinal().find("nsig-%d"%i).getVal())
		eSigsD[0]    = float(fitResultsD.floatParsFinal().find("nsig-%d"%i).getError())
		elSigsD[0]   = float(fitResultsD.floatParsFinal().find("nsig-%d"%i).getErrorLo())
		ehSigsD[0]   = float(fitResultsD.floatParsFinal().find("nsig-%d"%i).getErrorHi())
		nRelSigsD[0] = float(fitResultsD.floatParsFinal().find("nsig-%d"%i).getVal() / (fitResultsD.floatParsFinal().find("nsig-%d"%i).getVal() + fitResultsD.floatParsFinal().find("nbkg-%d"%i).getVal()))
		isSucT[0]    = int(not isFail)
		nSigsT[0]    = float(fitResultsT.floatParsFinal().find("nsig%d"%i).getVal())
		eSigsT[0]    = float(fitResultsT.floatParsFinal().find("nsig%d"%i).getError())
		elSigsT[0]   = float(fitResultsT.floatParsFinal().find("nsig%d"%i).getErrorLo())
		ehSigsT[0]   = float(fitResultsT.floatParsFinal().find("nsig%d"%i).getErrorHi())
		nRelSigsT[0] = float(fitResultsT.floatParsFinal().find("nsig%d"%i).getVal() / (fitResultsT.floatParsFinal().find("nsig%d"%i).getVal() + fitResultsT.floatParsFinal().find("nbkg%d"%i).getVal()))

		t.Fill()

		c1 = ROOT.TCanvas("c1%d"%i,"c1",800,1200)
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.6,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.2)
		pullPad = ROOT.TPad("pullPad","pullPad",0,0.2,1,0.4)
		dratioPad = ROOT.TPad("dratioPad","dratioPad",0,0.4,1,0.6)
		ROOT.gStyle.SetOptStat(0)
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		pullPad.UseCurrentStyle()
		dratioPad.UseCurrentStyle()
		plotPad.Draw()
		ratioPad.Draw()
		pullPad.Draw()
		dratioPad.Draw()
		plotPad.cd()

		frame = ws.var("mass_%s"%channel).frame()
		dataHist.plotOn(frame, ROOT.RooFit.Name("dataHist%d"%i))
		# ws.pdf("bkgpdf_%s"%channel).plotOn(frame, RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("bkgpdf_%s"%channel))
		# ws.pdf("bkgpdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("bkgpdf_%s_forFit%d"%(channel,i)))
		ws.pdf("model_pdf_%s_fixed%d"%(channel,i)).plotOn(frame, RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("model_pdf_%s"%channel))
		ws.pdf("model_pdf_%s_fixed%d"%(channel,i)).plotOn(frame, RooFit.Components("bkgpdf_%s"%(channel)), RooFit.LineColor(ROOT.kBlue), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("bkgpdf_%s"%(channel)))
		ws.pdf("model_pdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("model_pdf_%s_forFit%d"%(channel,i)))
		ws.pdf("model_pdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.Components("bkgpdf_%s_forFit%d"%(channel,i)), RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("bkgpdf_%s_forFit%d"%(channel,i)))
		# ws.pdf("model_pdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.Components("sig_pdf_%s_%d"%(channel,i)), RooFit.LineColor(ROOT.kGreen), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("sig_pdf_%s"%(channel)))
	
		hist_integ = nEvents  # dataHist.sum(False)
		hpull, hpull_max =  makePullHist(
			ws = ws,
			hist = frame.findObject("dataHist%d"%i),
			curve = frame.findObject("model_pdf_%s_forFit%d"%(channel,i)),
			norm = hist_integ/bkgpdf_integ
		)
		hpull.SetName("hpull%d"%i);
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])

		hratio, hratio_max =  makeRatioHist(
			ws = ws,
			hist = frame.findObject("dataHist%d"%i),
			curve = frame.findObject("model_pdf_%s_forFit%d"%(channel,i)),
			norm = hist_integ/bkgpdf_integ
		)
		hratio.SetName("hratio%d"%i);
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])

		xtmp = ROOT.Double(-999)
		ymin   = ROOT.Double(-999)
		ymax   = ROOT.Double(-999)
		frame.findObject("bkgpdf_%s"%channel).GetPoint(1,xtmp,ymax)
		ymin = 0
		#ymax = 1e5
		plotPad.DrawFrame(fitRange[0],ymin,fitRange[1],1.2*ymax,";m_{#mu#mu} [GeV]; Events / %.1f GeV" % binw)
		# plotPad.DrawFrame(fitRange[0],0.8*ymin,fitRange[1],2*ymax,";m_{#mu#mu} [GeV]; Events / %.1f GeV" % binw)
		frame.Draw("same")
		plotPad.SetLogy(0)

		ratioPad.cd()
		ratioPad.SetLogy(0)
		fr_ratioPad = ratioPad.DrawFrame(fitRange[0],0.8,fitRange[1],1.2,";; PDF_{bkg} gen/test")  # (F%d%d)"% (isFailD,isFail))
		fr_ratioPad.GetXaxis().SetLabelSize(0.1)
		fr_ratioPad.GetYaxis().SetLabelSize(0.08)
		fr_ratioPad.GetYaxis().SetTitleSize(0.12)
		fr_ratioPad.GetYaxis().SetTitleOffset(0.4)
		fr_ratioPad.Draw()
		ratio.Draw("L same")
		ratioD.Draw("L same")

		dratioPad.cd()
		dratioPad.SetLogy(0)
		fr_dratioPad = dratioPad.DrawFrame(fitRange[0],-1.3*hratio_max,fitRange[1],1.3*hratio_max,";;data/PDF_{bkg} -1")
		fr_dratioPad.GetXaxis().SetLabelSize(0.1)
		fr_dratioPad.GetYaxis().SetLabelSize(0.08)
		fr_dratioPad.GetYaxis().SetTitleSize(0.12)
		fr_dratioPad.GetYaxis().SetTitleOffset(0.4)
		fr_dratioPad.Draw()
		hratio.Draw("samepe0")

		pullPad.cd()
		pullPad.SetLogy(0)
		fr_pullPad = pullPad.DrawFrame(fitRange[0],-1.5*hpull_max,fitRange[1],1.5*hpull_max,";;pull(PDF_{bkg})") 
		fr_pullPad.GetXaxis().SetLabelSize(0.1)
		fr_pullPad.GetYaxis().SetLabelSize(0.08)
		fr_pullPad.GetYaxis().SetTitleSize(0.12)
		fr_pullPad.GetYaxis().SetTitleOffset(0.4)
		fr_pullPad.Draw()
		hpull.Draw("samepe0")

		if args.saveAll or (i<10):
			ROOT.gErrorIgnoreLevel = ROOT.kFatal
			c1.SaveAs( canvdir+"plot_%s_fix%d_nsigfac%d_%s_M%.0f_%d_%d--Job%d.pdf" %(channel,fix,args.nSigFac,args.func,massVal,nEvents,i,args.jobId), "pdf")
			# outfile = ROOT.TFile(canvdir+"plot_%s_fix%d_%s_M%.0f_%d_%d.root"%(channel,fix,args.func,massVal,nEvents,i), "RECREATE")
			# outfile.cd()
			# c1.Write()
			# hist.Write()
			# hratio.Write()
			# hpull.Write()
			# frame.findObject("dataHist%d"%i).Write()
			# outfile.Close()
			ROOT.gErrorIgnoreLevel = ROOT.kPrint


		c2 = ROOT.TCanvas("c2%d"%i,"c2",800,1200)
		plotPad2 = ROOT.TPad("plotPad2","plotPad2",0,0.6,1,1)
		ratioPad2 = ROOT.TPad("ratioPad2","ratioPad2",0,0.,1,0.2)
		pullPad2 = ROOT.TPad("pullPad2","pullPad2",0,0.2,1,0.4)
		dratioPad2 = ROOT.TPad("dratioPad2","dratioPad2",0,0.4,1,0.6)
		ROOT.gStyle.SetOptStat(0)
		plotPad2.UseCurrentStyle()
		ratioPad2.UseCurrentStyle()
		pullPad2.UseCurrentStyle()
		dratioPad2.UseCurrentStyle()
		plotPad2.Draw()
		ratioPad2.Draw()
		pullPad2.Draw()
		dratioPad2.Draw()
		plotPad2.cd()

		frame2 = ws.var("mass_%s"%channel).frame()
		dataHist.plotOn(frame2, ROOT.RooFit.Name("dataHist%d"%i))
		# ws.pdf("bkgpdf_%s"%channel).plotOn(frame, RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("bkgpdf_%s"%channel))
		# ws.pdf("bkgpdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("bkgpdf_%s_forFit%d"%(channel,i)))
		ws.pdf("model_pdf_%s_fixed%d"%(channel,i)).plotOn(frame2, RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("model_pdf_%s"%channel))
		ws.pdf("model_pdf_%s_fixed%d"%(channel,i)).plotOn(frame2, RooFit.Components("bkgpdf_%s"%(channel)), RooFit.LineColor(ROOT.kBlue), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("bkgpdf_%s"%(channel)))
		ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).plotOn(frame2, RooFit.LineColor(ROOT.kGreen+2), ROOT.RooFit.Name("model_pdf_%s_forFit%d"%(channel,-1*i)))
		ws.pdf("model_pdf_%s_forFit%d"%(channel,-1*i)).plotOn(frame2, RooFit.Components("bkgpdf_%s_forFit%d"%(channel,-1*i)), RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("bkgpdf_%s_forFit%d"%(channel,-1*i)))
		# ws.pdf("model_pdf_%s_forFit%d"%(channel,i)).plotOn(frame, RooFit.Components("sig_pdf_%s_%d"%(channel,i)), RooFit.LineColor(ROOT.kGreen), RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("sig_pdf_%s"%(channel)))
	
		hist_integ = nEvents  # dataHist.sum(False)
		hpull2, hpull_max2 =  makePullHist(
			ws = ws,
			hist = frame2.findObject("dataHist%d"%i),
			curve = frame2.findObject("model_pdf_%s_forFit%d"%(channel,-1*i)),
			norm = hist_integ/bkgpdf_integ
		)
		hpull2.SetName("hpull%d"%-1*i);
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])

		hratio2, hratio_max2 =  makeRatioHist(
			ws = ws,
			hist = frame2.findObject("dataHist%d"%i),
			curve = frame2.findObject("model_pdf_%s_forFit%d"%(channel,-1*i)),
			norm = hist_integ/bkgpdf_integ
		)
		hratio2.SetName("hratio%d"%-1*i);
		ws.var("mass_%s"%channel).setBins(nBinsForFit)
		ws.var("mass_%s"%channel).setRange(fitRange[0],fitRange[1])

		xtmp = ROOT.Double(-999)
		ymin   = ROOT.Double(-999)
		ymax   = ROOT.Double(-999)
		frame2.findObject("bkgpdf_%s"%channel).GetPoint(1,xtmp,ymax)
		ymin = 0
		#ymax = 1e5
		plotPad2.DrawFrame(fitRange[0],ymin,fitRange[1],1.2*ymax,";m_{#mu#mu} [GeV]; Events / %.1f GeV" % binw)
		# plotPad.DrawFrame(fitRange[0],0.8*ymin,fitRange[1],2*ymax,";m_{#mu#mu} [GeV]; Events / %.1f GeV" % binw)
		frame2.Draw("same")
		plotPad2.SetLogy(0)

		ratioPad2.cd()
		ratioPad2.SetLogy(0)
		fr_ratioPad2 = ratioPad2.DrawFrame(fitRange[0],0.8,fitRange[1],1.2,";; PDF_{bkg} gen/default")  # (F%d%d)"% (isFailD,isFail))
		fr_ratioPad2.GetXaxis().SetLabelSize(0.1)
		fr_ratioPad2.GetYaxis().SetLabelSize(0.08)
		fr_ratioPad2.GetYaxis().SetTitleSize(0.12)
		fr_ratioPad2.GetYaxis().SetTitleOffset(0.4)
		fr_ratioPad2.Draw()
		ratio.Draw("L same")
		ratioD.Draw("L same")

		dratioPad2.cd()
		dratioPad2.SetLogy(0)
		fr_dratioPad2 = dratioPad2.DrawFrame(fitRange[0],-1.3*hratio_max,fitRange[1],1.3*hratio_max,";;data/PDF_{bkg} -1")
		fr_dratioPad2.GetXaxis().SetLabelSize(0.1)
		fr_dratioPad2.GetYaxis().SetLabelSize(0.08)
		fr_dratioPad2.GetYaxis().SetTitleSize(0.12)
		fr_dratioPad2.GetYaxis().SetTitleOffset(0.4)
		fr_dratioPad2.Draw()
		hratio2.Draw("samepe0")

		pullPad2.cd()
		pullPad2.SetLogy(0)
		fr_pullPad2 = pullPad2.DrawFrame(fitRange[0],-1.5*hpull_max,fitRange[1],1.5*hpull_max,";;pull(PDF_{bkg})") 
		fr_pullPad2.GetXaxis().SetLabelSize(0.1)
		fr_pullPad2.GetYaxis().SetLabelSize(0.08)
		fr_pullPad2.GetYaxis().SetTitleSize(0.12)
		fr_pullPad2.GetYaxis().SetTitleOffset(0.4)
		fr_pullPad2.Draw()
		hpull2.Draw("samepe0")

		if args.saveAll or (i<10):
			ROOT.gErrorIgnoreLevel = ROOT.kFatal
			c2.SaveAs( canvdir+"plot_%s_fix%d_nsigfac%d_%s_M%.0f_%d_%d--Job%d.pdf" %(channel,fix,args.nSigFac,"default",massVal,nEvents,i,args.jobId), "pdf")
			# outfile = ROOT.TFile(canvdir+"plot_%s_fix%d_%s_M%.0f_%d_%d.root"%(channel,fix,args.func,massVal,nEvents,i), "RECREATE")
			# outfile.cd()
			# c1.Write()
			# hist.Write()
			# hratio.Write()
			# hpull.Write()
			# frame.findObject("dataHist%d"%i).Write()
			# outfile.Close()
			ROOT.gErrorIgnoreLevel = ROOT.kPrint

		sys.stdout.flush()
		gc.collect()

	# -- Bias -- #
	rootdir = "results_%s/"%args.func
	if not os.path.isdir(rootdir):
		os.makedirs(rootdir)

	cname = "%s_fix%d_nsigfac%d_%s_M%.0f_nToys%d_nEv%d--Job%d" %(channel,fix,args.nSigFac,args.func,massVal,nToys,nEventsFix,args.jobId)
	f_bias = ROOT.TFile( rootdir+"/nsig_%s.root"%cname, 'recreate' )
	f_bias.cd()
	t.Write()
	f_bias.Close()
