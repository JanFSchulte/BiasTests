import ROOT
import sys
def init_fun_sigma_BB():
    fun= ROOT.TF1('sigma_BB', '(x<[0])*(sqrt(([1]*[1]/x)+([2]*[2]/(x*x))+[3]*[3])+[4]*x)+(x>=[0])*([5]+[6]*x+[7]*x*x)')
    fun_name='(x<p0)*(#sqrt{#frac{(p1)^{2}}{x}+#frac{(p2)^{2}}{x^{2}}+(p3)^{2}}+p4*x)+(x#geqp0)*(p5+p6*x+p7*x^{2})'
    fun.SetParameter(0,1200)
    fun.SetParameter(1,6.026)
    fun.SetParameter(2,136)
    fun.SetParameter(3,1.022)
    fun.SetParameter(4,0.0001271)
    fun.SetParameter(5,1.244)
    fun.SetParameter(6,-6.014e-5)
    fun.SetParameter(7,1.356e-8)
    return fun
def init_fun_sigma_BE():
    fun= ROOT.TF1('sigma_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])+[3]*x')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}+p3*x'
    fun.SetParameter(0,-20.25)
    fun.SetParameter(1,30.1)
    fun.SetParameter(2,1.514)
    fun.SetParameter(3,4.986e-5)
    return fun
def init_fun_PowerL_BE():
    fun= ROOT.TF1('PowerL_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])+[3]*x')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}+p3*x'
    fun.SetParameter(0,36.4)
    fun.SetParameter(1,-9.41)
    fun.SetParameter(2,1.186)
    fun.SetParameter(3,8.588e-5)
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
def init_fun_pol2(p0,p1,p2):
    fun= ROOT.TF1('POL2', '[0]+[1]*x+[2]*x*x)')
    fun_name='p0+p1*x+p2*x^{2}'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
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
def init_fun_pol3_pol0(p0,p1,p2,p3,p4):
    fun= ROOT.TF1('line2' , '(x<[0])*([1]+[2]*x+[3]*x*x)+(x>=[0])*([4])')
    fun_name='(x<p0)*(p1+p2*x+p3*x^{2})+(x#geqp0)*(p4)'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    return fun
BB_sigma =init_fun_sigma_BB()
BE_sigma =init_fun_sigma_BE()
BB_mean  =init_fun_pol1_pol2(400,9.94850e-01,4.50685e-06,9.97197e-01,-7.96237e-07,-1.28048e-10)
BE_mean  =init_fun_pol1_pol2(400,0.9874,1.71e-5,0.9935,1.849e-6,-2.28e-10)
BB_PowerR=init_fun_pol3_pol3(1300,-15.42,0.1783,-0.0001144,-17.3,0.03644,-4.719e-6)
BB_PowerL=init_fun_pol3_pol3(1250,3.778,-0.004699,2.629e-6,2.803,-0.0008787,2.021e-7)
BB_CutL  =init_fun_pol1_pol2(700,1.003,0.0008312,1.465,0.0002005,-7.185e-8)
BB_CutR  =init_fun_pol3_pol3(610,2.426,-0.003341,5.153e-6,2.233,0.0001188,-4.319e-8)
BE_PowerL=init_fun_PowerL_BE()
BE_PowerR=init_fun_pol3_pol0(7.20000e+02,-8.90230e+01,7.00910e-01,-7.91725e-04,1.65938e+01)
BE_CutL  =init_fun_pol3_pol3(2600,1.337,0.0006109,-2.029e-7,-0.6927,0.001364,-1.911e-7)
BE_CutR  =init_fun_pol3_pol3(1600,3.952,-0.006593,5.119e-6,15.74,-0.007326,1.012e-6)


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
            self.PowerR=BB_PowerR.Eval(mass)
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
dCB.get_value(5000,False) ##True for BB, False for BE
print "mean=%f, sigma=%f, PowerL=%f, PowerR=%f, CutL=%f, CutR=%f"%(dCB.mean,dCB.sigma,dCB.PowerL,dCB.PowerR,dCB.CutL,dCB.CutR)
###### BCD parameter #############    
