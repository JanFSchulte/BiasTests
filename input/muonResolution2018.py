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
   a_BB=1.65
   b_BB=-5e-05
   c_BB=2.88e-09
   parameters['alphaL']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass

   #####alphaL BE
   a_BE= 1.25
   b_BE= -5e-07
   c_BE=2.19e-10
   parameters['alphaL']['BE']= a_BE + b_BE*mass + c_BE*mass*mass

   #####alphaR BB
   a_BB= 1.97
   b_BB= 1e-05
   c_BB= -1e-08
   d_BB= 1e-14
   parameters['alphaR']['BB']=a_BB + b_BB*mass + c_BB*mass*mass + d_BB*mass**3

   #####alphaR BE
   a_BE= 1.33
   b_BE= 1e-05
   c_BE= -1.97e-09
   d_BE= 1e-14
   parameters['alphaR']['BE']=a_BE + b_BE*mass + c_BE*mass*mass + d_BE*mass**3

   #####nL BB
   a_BB=1.2
   b_BB=0.000177
   c_BB=-2.79e-08
   parameters['nL']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass

   #####nL BE
   a_BE= 2.54
   b_BE= -0.00025
   c_BE=4.46e-08
   parameters['nL']['BE']= a_BE + b_BE*mass + c_BE*mass*mass

   #####nR BB
   parameters['nR']['BB']=20

   #####nR BE
   parameters['nR']['BE']=20

   #####scale BB
   a_BB = -0.000359
   b_BB = -5.3e-06
   c_BB = 2.29e-09
   d_BB = -4.89e-13
   e_BB = 3.87e-17
   parameters['scale']['BB']= a_BB + b_BB*mass + c_BB*mass**2 + d_BB*mass**3 + e_BB*mass**4

   #####scale BE
   a_BE = -0.00161
   b_BE = -6.37e-06
   c_BE = 2.18e-09
   d_BE = -4.37e-13
   e_BE = 3.83e-17
   parameters['scale']['BE']= a_BE + b_BE*mass + c_BE*mass**2 + d_BE*mass**3 + e_BE*mass**4


   #####sigma BB
   a_BB=0.00608
   b_BB=3.42e-05
   c_BB=-1.34e-08
   d_BB=2.4e-12
   e_BB=-1.5e-16
   parameters['sigma']['BB'] = a_BB + b_BB*mass + c_BB*mass*mass + d_BB*mass**3 + e_BB*mass**4

   #####sigma BE
   a_BE=0.0134
   b_BE=3.92e-05
   c_BE=-1.42e-08
   d_BE=2.46e-12
   e_BE=-1.5e-16
   parameters['sigma']['BE'] = a_BE + b_BE*mass + c_BE*mass*mass + d_BE*mass**3 + e_BE*mass**4


   return parameters




