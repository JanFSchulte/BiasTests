import ROOT
import math
#from ROOT import RooFit, RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData, RooArgSet
from ROOT import *
ROOT.gROOT.SetBatch(ROOT.kTRUE) 

import sys,getopt
from math import sqrt

#def getResolution(mass):
#   return res_scale*0.01*(sqrt(res_s*res_s/mass+res_n*res_n/mass/mass+res_c*res_c)+res_slope*mass)

def getResolution(mass):
   parameters={}
   parameters['alphaL']={}
   parameters['alphaR']={}
   parameters['nL']={}
   parameters['nR']={}
   parameters['scale'] ={}
   parameters['sigma'] ={}
   #####alphaL BB
   a_BB=1.67
   b_BB=-5e-05
   c_BB=1.86e-09
   parameters['alphaL']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass

   #####alphaL BE
   a_BE= 1.22
   b_BE= -5e-07
   c_BE=3.04e-09
   parameters['alphaL']['BE']= a_BE + b_BE*mass + c_BE*mass*mass

   #####alphaR BB
   a_BB= 2
   b_BB= 1e-07
   c_BB= -1e-08
   d_BB= 1.66e-13
   parameters['alphaR']['BB']=a_BB + b_BB*mass + c_BB*mass*mass + d_BB*mass**3

   #####alphaR BE
   a_BE= 1.27
   b_BE= 1e-05
   c_BE= -1e-10
   d_BE= 3.67e-13
   parameters['alphaR']['BE']=a_BE + b_BE*mass + c_BE*mass*mass + d_BE*mass**3

   #####nL BB
   a_BB=1.27
   b_BB=4.48e-05
   c_BB=-7.26e-09
   parameters['nL']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass

   #####nL BE
   a_BE= 2.48
   b_BE= -0.000354
   c_BE=6.16e-08
   parameters['nL']['BE']= a_BE + b_BE*mass + c_BE*mass*mass

   #####nR BB
   parameters['nR']['BB']=20

   #####nR BE
   parameters['nR']['BE']=20

   #####scale BB
   a_BB = -4.74e-05
   b_BB = -5.48e-06
   c_BB = 2.36e-09
   d_BB = 3.82e-17
   parameters['scale']['BB']= a_BB + b_BB*mass + c_BB*mass**2 + d_BB*mass**3

   #####scale BE
   a_BE = -0.000131
   b_BE = -8.22e-06
   c_BE = 2.92e-09
   d_BE = -5e-13
   parameters['scale']['BE']= a_BB + b_BB*mass + c_BB*mass**2 + d_BB*mass**3


   #####sigma BB
   a_BB=0.00606
   b_BB=3.41e-05
   c_BB=-1.33e-08
   d_BB=2.39e-12
   e_BB=-1.5e-16
   parameters['sigma']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass + d_BB*mass**3 + e_BB*mass**4

   #####sigma BE
   a_BE=0.0129
   b_BE=4.02e-05
   c_BE=-1.45e-08
   d_BE=2.49e-12
   e_BE=-1.48e-16
   parameters['sigma']['BE'] = a_BE + b_BE*mass + c_BE*mass*mass + d_BE*mass**3 + e_BE*mass**4


   return parameters




