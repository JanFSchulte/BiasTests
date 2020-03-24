import ROOT
import sys

def init_fun_sigma_BB(p0,p1,p2,p3,p4,p5,p6):
    fun= ROOT.TF1('sigma_BB', '(x<[0])*(sqrt(([1]*[1]/x)+([2]*[2]/(x*x))+[3]*[3])+[4]*x)+(x>=[0])*([5]+[6]*x)')
    fun_name='(x<p0)*(#sqrt{#frac{(p1)^{2}}{x}+#frac{(p2)^{2}}{x^{2}}+(p3)^{2}}+p4*x)+(x#geqp0)*(p5+p6*x)'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    fun.SetParameter(5,p5)
    fun.SetParameter(6,p6)
    return fun
def init_fun_sigma_BE(p0,p1,p2):
    fun= ROOT.TF1('sigma_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    return fun
def init_fun_pol1_pol2(p0,p1,p2,p3,p4,p5):
    fun= ROOT.TF1('mean_BE', '(x<[0])*([1]+[2]*x)+(x>=[0])*([3]+[4]*x+[5]*x*x)')
    fun_name='(x<p0)*(p1+p2*x)+(x#geqp0)*(p3+p4*x+p5*x^{2})'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    fun.SetParameter(5,p5)
    return fun
def init_fun_pol3_pol2_v1(p0,p1,p2,p3,p4,p5):
    fun= ROOT.TF1('mean_BE', '(x<[0])*([1]+[2]*x+[3]*x*x)+(x>=[0])*([4]+[5]*x)')
    fun_name='(x<p0)*(p1+p2*x+p3*x^{2})+(x#geqp0)*(p4+p5*x)'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    fun.SetParameter(5,p5)
    return fun
def init_fun_pol3_pol3(p0,p1,p2,p3,p4,p5,p6):
    fun= ROOT.TF1('line2' , '(x<[0])*([1]+[2]*x+[3]*x*x)+(x>=[0])*([4]+[5]*x+[6]*x*x)')
    fun_name='(x<p0)*(p1+p2*x+p3*x^{2})+(x#geqp0)*(p4+p5*x+p6*x^{2})'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    fun.SetParameter(5,p5)
    fun.SetParameter(6,p6)
    return fun
def init_fun_pol1(p0,p1):
    fun= ROOT.TF1('POL1' , '[0]+[1]*x')
    fun_name='p0+p1*x'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    return fun
def init_fun_PowerL_BE(p0,p1,p2,p3):
    fun= ROOT.TF1('PowerL_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])+[3]*x')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}+p3*x'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    return fun


BB_sigma =init_fun_sigma_BB(3900,1.44900e+01,5.16499e-02,1.17605e+00,1.75522e-05,7.07217e-01,1.42834e-04)
BE_sigma =init_fun_sigma_BE(-2.05926e+01,8.71924e+01,1.96286e+00)
BB_mean  =init_fun_pol1_pol2(4.00000e+02,9.94251e-01,6.84913e-06,9.97039e-01,3.70794e-07,-1.91138e-10)
BE_mean  =init_fun_pol1_pol2(7.00000e+02,9.94276e-01,6.61709e-06,9.98535e-01,1.90473e-07,4.87478e-11)
BB_PowerR=init_fun_pol3_pol2_v1(2.00000e+03,5.19673e+01,-4.60266e-02,2.31366e-05,9.14599e+01,-1.79042e-02)
BB_PowerL=init_fun_sigma_BE(-3.89294e+01,8.38611e-03,1.27215e+00)
BB_CutL  =init_fun_pol1_pol2(6.10000e+02,1.05013e+00,1.09063e-03,1.65652e+00,6.04812e-05,3.48820e-09)
BB_CutR  =init_fun_pol3_pol3(1.50000e+03,1.61420e+00,8.68022e-04,-4.17800e-07,1.21676e+00,7.03864e-04,-1.29203e-07)
BE_PowerL=init_fun_PowerL_BE(-1.82292e+01,3.97244e+02,1.43349e+00,-2.99232e-06)
BE_PowerR=init_fun_pol3_pol3(1.20000e+03,1.02902e+01,-6.73964e-02,8.86873e-05,1.17254e+02,-5.64383e-02,7.00598e-06)
BE_CutL  =init_fun_pol3_pol3(2.20000e+03,1.45486e+00,7.69878e-04,-2.99246e-07,-1.47734e-01,1.23381e-03,-1.76108e-07)
BE_CutR  =init_fun_pol1(2.75129e+00,-5.35424e-05)


class DCB_para:
    def __init__(self,name):
        self.name=name
        self.sigma =0 
        self.mean  =0 
        self.PowerL=0
        self.PowerR=0
        self.CutL  =0
        self.CutR  =0
    def get_value(self,mass,isBB):
        if isBB:
            self.sigma =BB_sigma .Eval(mass)/100
            self.mean  =BB_mean  .Eval(mass)-1
            self.PowerR=BB_PowerR.Eval(mass) if mass < 4750 else 7.08
            self.PowerL=BB_PowerL.Eval(mass)
            self.CutL  =BB_CutL  .Eval(mass)
            self.CutR  =BB_CutR  .Eval(mass)
        else:
            self.sigma =BE_sigma .Eval(mass)/100
            self.mean  =BE_mean  .Eval(mass)-1
            self.PowerR=BE_PowerR.Eval(mass)
            self.PowerL=BE_PowerL.Eval(mass)
            self.CutL  =BE_CutL  .Eval(mass)
            self.CutR  =BE_CutR  .Eval(mass)
######### main #############
dCB=DCB_para("dcb")
dCB.get_value(5000,True) ##True for BB, False for BE
print "mean=%f, sigma=%f, PowerL=%f, PowerR=%f, CutL=%f, CutR=%f"%(dCB.mean,dCB.sigma,dCB.PowerL,dCB.PowerR,dCB.CutL,dCB.CutR)
###### BCD parameter #############    
