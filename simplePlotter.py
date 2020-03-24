import os, sys
import gc
from array import array
import numpy as np
import math

def getVariable( t, tt ):
	if tt == 0:
		return ( t.nSigsT )
	elif tt==1:
		return (t.nSigsT - t.nSigsD)
	elif tt==2:
		return (t.nSigsT / t.eSigsT)
	elif tt==3:
		return ( (t.nSigsT - t.nSigsD) / math.sqrt( abs((t.eSigsT)**.2 - (t.eSigsD)**2.)) )
	elif tt==4:
		return ( (t.nSigsT - t.nSigsD) / math.sqrt( t.nEv ) )
	elif tt==5:
		return ( t.nSigsD )
	elif tt==6:
		return (t.nSigsD / t.eSigsD)
	else:
		return -999

def getTitle( tt ):
	if tt == 0:
		return "Nsig_{test}"
	elif tt==1:
		return "Nsig_{test} - Nsig_{default}"
	elif tt==2:
		return "Nsig_{test} / #sigma_{test}"
	elif tt==3:
		return "(Nsig_{test} - Nsig_{default}) / (#sigma_{test}^{2} - #sigma_{default}^{2})"
	elif tt==4:
		return "(Nsig_{test} - Nsig_{default}) / (stat. unc.)"
	elif tt==5:
		return "Nsig_{default}"
	elif tt==6:
		return "Nsig_{default} / #sigma_{default}"
	else:
		return "Not defined type: %d" % tt

def getXRange( tt, nsigfac = 0, nEv = -1 ):
	# xmax = 500 # 0.05*nEv if nEv > 0 else 200
	# xmaxnorm = 5 # 2*xmax / math.sqrt(nEv) if nEv > 0 else 5

	nsig = nsigfac * math.sqrt(nEv) if nEv > 0 else 0

	if tt == 0:
		return 100,-500+nsig,500+nsig
	elif tt==1:
		return 100,-200,200
	elif tt==2:
		return 100,-5+nsigfac,5+nsigfac
	elif tt==3:
		return 100,-5,5
	elif tt==4:
		return 100,-2,2
	elif tt==5:
		return 100,-500+nsig,500+nsig
	elif tt==6:
		return 100,-5+nsigfac,5+nsigfac
	else:
		return -999,-999,-999

def getYRange( tt, nsigfac = 0, nEv = -1 ):
	nsig = nsigfac * math.sqrt(nEv) if nEv > 0 else 0

	if tt == 0:
		return 50+nsig,50+nsig
	elif tt==1:
		return 50,50
	elif tt==2:
		return 1+nsigfac,1+nsigfac
	elif tt==3:
		return 10,10
	elif tt==4:
		return 10,10
	elif tt==5:
		return 50+nsig,50+nsig
	elif tt==6:
		return 1+nsigfac,1+nsigfac
	else:
		return -999,-999

def removeOutliers( biases, q = -1 ):
	out = []

	a = biases

	ap = []
	for i in a:
		ap.append(abs(i))
	ap.sort()

	if q > 0:
		qi = int(len(ap)*q)
		slicer = ap[ qi ]
	else:
		slicer = 1e99

	for b in a:
		if abs(b) < slicer:
			out.append(b)

	return out

if __name__ == '__main__':
	import ROOT
	from ROOT import RooFit
	from ROOT import RooRealVar, RooMCStudy, RooNumIntConfig
	from ROOT import TF1, TRandom3, TGraph, TH1F, TH1D, TGraph, TGraphAsymmErrors, TTree
	ROOT.gROOT.SetBatch(True)
	ROOT.TH1.AddDirectory(False)
	ROOT.TH1.SetDefaultSumw2(True)


	fixtag = 1  # 0 1
	fixstr = "parm. fixed" if fixtag == 1 else "parm. float"

	funcs = [
		("default",    fixtag, ROOT.kBlack,"Z' shape (%s)"%fixstr),
		# ("defaultlow", fixtag, ROOT.kBlack,"Z' shape(%s)"%fixstr),
		("MRBW",       fixtag, ROOT.kBlue,"MRBW (%s)"%fixstr),
		("MBWZg",      fixtag, ROOT.kRed,"MBWZg (%s)"%fixstr)
	]

	# -- 4x
	masses = [
		(200, 90591),
		(250, 97763),
		(300, 80939),
		(350, 53204),
		(400, 36330),
		(500, 18678)
		# (600, xxx)
	]
#	masses = [
#		(200, 20662),
#		(250, 11621),
#		(300, 6903),
#		(350, 4576),
#		(400, 3096),
#		(500, 1526)
		# (600, xxx)
	#]

	# -- 6x
	# masses = [
	# 	(200, 31687),
	# 	(250, 17929),
	# 	(300, 10740),
	# 	(350, 7081),
	# 	(400, 4843),
	# 	(500, 2428)
	# 	# (600, 1426)
	# ]
	nMasses = len(masses)

	types = [
		0,  # nSigT
		1,  # nSigT - nSigD
		2,  # nSigT / unc
		3,  # nSigT - nSigD / unc
		4,  # nSigT - nSigD / stat. unc.
		5,  # nSigD
		6   # nSigD / unc
	]

	nsigfacs = [0,3]  # [ 0, 3 ]

	for nsigfac in nsigfacs:


		for tt in types:

			title = getTitle(tt)

			c1 = ROOT.TCanvas("c%d"%tt,"c%d"%tt)
			plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
			ROOT.gStyle.SetOptStat(0)
			plotPad.UseCurrentStyle()
			plotPad.SetLogy(False)
			plotPad.Draw()
			plotPad.cd()

			g  = [TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses)]
			gm = [TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses),TGraphAsymmErrors(nMasses)]

			for iff, (ff, fix, cc, ss) in enumerate(funcs):
				g[iff].SetName(ff)
				g[iff].SetLineColor(cc)
				gm[iff].SetName("mean_%s"%ff)
				gm[iff].SetLineColor(cc)
				# medians = {}
				for i, (m, nev) in enumerate(masses):
					tag = "%s_fix%d_nsigfac%d_M%d_nEv%d_type%d" % (ff, fix, nsigfac, m, nev, tt)
					fn = "results/nsig_dimuon_Legacy2018_BB_fix%d_nsigfac%d_%s_M%d_nToys200_nEv%d.root" % (fix, nsigfac, ff, m, nev)
					f = ROOT.TFile( fn )
					# print f

					xrange = getXRange(tt, nsigfac, nev)
					hb = TH1D("hb_%s"%tag, "%s M%d nEv %d" % (ff, m, nev), xrange[0], xrange[1], xrange[2])
					biases = []

					t = f.Get("t").Clone()
					for toy in t:
						if t.isSucD and t.isSucT:
							var = getVariable( t, tt )
							biases.append( var )
							hb.Fill( var )
					biases.sort()
					n = len(biases)
					print n
					med  = np.median(biases)
					mean = np.mean(biases)  # removeOutliers(biases)
					std  = np.std(biases)
					mean_rms = mean / std
					nToys = len(biases)

					g[iff].SetPoint(i,m,med)
					gm[iff].SetPoint(i,m,mean)

					print "[%s, %d, %d]:" % (ff, m, nev), med, mean

					c2 = ROOT.TCanvas("c_%s"%tag,"")
					plotPad2 = ROOT.TPad("plotPad2","plotPad2",0,0,1,1)
					ROOT.gStyle.SetOptStat(0)
					plotPad2.UseCurrentStyle()
					plotPad2.SetLogy(False)
					plotPad2.Draw()
					plotPad2.cd()
					ymax = 1.5*hb.GetMaximum()
					fr_plotPad2 = plotPad2.DrawFrame(xrange[1],0,xrange[2],ymax,";%s; # datasets"%title)
					fr_plotPad2.GetYaxis().SetTitleOffset(1.3)
					fr_plotPad2.GetYaxis().SetLabelSize(0.025)
					fr_plotPad2.Draw()
					hb.Draw("SAME HIST")
					latex2 = ROOT.TLatex()
					# latex2.DrawLatexNDC(0.15, 0.86, '#bf{#scale[0.8]{Dataset PDF: default}}')
					latex2.DrawLatexNDC(0.15, 0.82, '#bf{#scale[0.8]{Fit(bkg) PDF: %s}}' % ss)
					latex2.DrawLatexNDC(0.15, 0.78, '#bf{#scale[0.8]{%d GeV}}' % m)
					latex2.DrawLatexNDC(0.15, 0.74, '#bf{#scale[0.8]{# datasets=%d, # bkg/dataset=%d, # inj. sig=%d}}' % (nToys, nev, int(nsigfac*math.sqrt(float(nev)))) )
					latex2.DrawLatexNDC(0.15, 0.70, '#bf{#scale[0.8]{Median = %.2lg, Mean = %.2lg}}' % (med, mean) )
					c2.Print("./dists/plot_%s.pdf"%tag)
					c2.Print("./dists/plot_%s.png"%tag)

			nsigstr = "%d*#sigma_{stat.}"%nsigfac if nsigfac > 0 else "0"

			plotPad.cd()
			ymax = getYRange(tt, nsigfac, masses[0][1])[0]
			if nsigfac > 0 and (tt == 0 or tt == 2 or tt == 5 or tt == 6):
				ymax = 2*ymax
				ymin = 0
			else:
				ymin = -1*ymax
			fr_plotPad = plotPad.DrawFrame(masses[0][0],ymin,masses[-1][0],ymax,";m_{#mu#mu} [GeV]; median of %s"%title)
			fr_plotPad.GetYaxis().SetTitleOffset(1.3)
			fr_plotPad.GetYaxis().SetLabelSize(0.025)
			fr_plotPad.SetLineWidth(0)
			fr_plotPad.SetLineColorAlpha(0,0)
			fr_plotPad.Draw()

			g[0].Draw("L SAME")
			g[1].Draw("L SAME")
			g[2].Draw("L SAME")
			# g[3].Draw("L SAME")

			leg=ROOT.TLegend(0.5,0.65,0.9,0.85)
			leg.SetBorderSize(0)
			leg.SetLineWidth(0)
			leg.SetLineStyle(0)
			leg.SetFillStyle(0)
			leg.SetLineColor(0)
			leg.AddEntry(g[0], funcs[0][3],"l")
			leg.AddEntry(g[1], funcs[1][3],"l")
			leg.AddEntry(g[2], funcs[2][3],"l")
			# leg.AddEntry(g[3], funcs[3][3],"l")
			leg.Draw()

			latex = ROOT.TLatex()
			latex.DrawLatexNDC(0.50, 0.86, "#bf{#scale[0.8]{Test PDF (dataset gen.):}}")
			latex.DrawLatexNDC(0.15, 0.86, "#bf{#scale[0.8]{Default PDF: Z' shape}}")
			latex.DrawLatexNDC(0.15, 0.81, "#bf{#scale[0.8]{# datasets / mass = 2000}}" )
			latex.DrawLatexNDC(0.15, 0.76, "#bf{#scale[0.8]{# inj. sig = %s}}"%nsigstr )

			c1.SaveAs("./plots/bias_median_default_fix%d_nsigfac%d_type%d.pdf"%(fixtag, nsigfac, tt),"pdf")


			plotPad.cd()
			ymax = getYRange(tt, nsigfac, masses[0][1])[1]
			if nsigfac > 0 and (tt == 0 or tt == 2 or tt == 5 or tt == 6):
				ymax = 2*ymax
				ymin = 0
			else:
				ymin = -1*ymax
			fr_plotPad = plotPad.DrawFrame(masses[0][0],ymin,masses[-1][0],ymax,";m_{#mu#mu} [GeV]; mean of %s"%title)
			fr_plotPad.GetYaxis().SetTitleOffset(1.3)
			fr_plotPad.GetYaxis().SetLabelSize(0.025)
			fr_plotPad.SetLineWidth(0)
			fr_plotPad.SetLineColorAlpha(0,0)
			fr_plotPad.Draw()

			gm[0].Draw("L SAME")
			gm[1].Draw("L SAME")
			gm[2].Draw("L SAME")
			# gm[3].Draw("L SAME")

			leg=ROOT.TLegend(0.5,0.65,0.9,0.85)
			leg.SetBorderSize(0)
			leg.SetLineWidth(0)
			leg.SetLineStyle(0)
			leg.SetFillStyle(0)
			leg.SetLineColor(0)
			leg.AddEntry(gm[0], funcs[0][3],"l")
			leg.AddEntry(gm[1], funcs[1][3],"l")
			leg.AddEntry(gm[2], funcs[2][3],"l")
			# leg.AddEntry(gm[3], funcs[3][3],"l")
			leg.Draw()

			latex = ROOT.TLatex()
			latex.DrawLatexNDC(0.50, 0.86, "#bf{#scale[0.8]{Test PDF (dataset gen.):}}")
			latex.DrawLatexNDC(0.15, 0.86, "#bf{#scale[0.8]{Default PDF: Z' shape}}")
			latex.DrawLatexNDC(0.15, 0.81, "#bf{#scale[0.8]{# datasets / mass = 2000}}" )
			latex.DrawLatexNDC(0.15, 0.76, "#bf{#scale[0.8]{# inj. sig = %s}}"%nsigstr )

			c1.SaveAs("./plots/bias_mean_default_fix%d_nsigfac%d_type%d.pdf"%(fixtag, nsigfac, tt),"pdf")
