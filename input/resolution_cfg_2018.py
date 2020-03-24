import ROOT
import sys

def init_fun_sigma_BB(p0,p1,p2,p3,p4,p5):
    fun= ROOT.TF1('sigma_BB', '(x<[0])*(sqrt(([1]*[1]/x)+([2]*[2]/(x*x))+[3]*[3]))+(x>=[0])*([4]+[5]*x)')
    fun_name='(x<p0)*(#sqrt{#frac{(p1)^{2}}{x}+#frac{(p2)^{2}}{x^{2}}+(p3)^{2}})+(x#geqp0)*(p4+p5*x)'
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
def init_fun_sigma_BE(p0,p1,p2):
    fun= ROOT.TF1('sigma_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    return fun
def init_fun_PowerL_BE(p0,p1,p2,p3):
    fun= ROOT.TF1('PowerL_BE', 'sqrt(([0]*[0]/x)+([1]*[1]/(x*x))+[2]*[2])+[3]*x')
    fun_name='#sqrt{#frac{(p0)^{2}}{x}+#frac{(p1)^{2}}{x^{2}}+(p2)^{2}}+p3*x'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    return fun

# linear interpolation for each point
#            (x<[2]) *( (([3] -[1]) /([2] -[0])) *(x-[2])  + [3] ) +
# (x>=[2]  && x<[4]) *( (([5] -[3]) /([4] -[2])) *(x-[4])  + [5] ) +
# (x>=[4]  && x<[6]) *( (([7] -[5]) /([6] -[4])) *(x-[6])  + [7] ) +
# (x>=[6]  && x<[8]) *( (([9] -[7]) /([8] -[6])) *(x-[8])  + [9] ) +
# (x>=[8]  && x<[10])*( (([11]-[9]) /([10]-[8])) *(x-[10]) + [11] ) +
# (x>=[10] && x<[12])*( (([13]-[11])/([12]-[10]))*(x-[12]) + [13] ) +
# (x>=[12])          *( (([15]-[13])/([14]-[12]))*(x-[14]) + [15] )
def init_fun_interpolation(p0,p1, p2,p3, p4,p5, p6,p7, p8,p9, p10,p11, p12,p13, p14,p15):
    fun= ROOT.TF1('interpolation' , '(x<[2])*((([3]-[1])/([2]-[0]))*(x-[2])+[3])+(x>=[2]&&x<[4])*((([5]-[3])/([4]-[2]))*(x-[4])+[5])+(x>=[4]&&x<[6])*((([7]-[5])/([6]-[4]))*(x-[6])+[7])+(x>=[6]&&x<[8])*((([9]-[7])/([8]-[6]))*(x-[8])+[9])+(x>=[8]&&x<[10])*((([11]-[9])/([10]-[8]))*(x-[10])+[11])+(x>=[10]&&x<[12])*((([13]-[11])/([12]-[10]))*(x-[12])+[13])+(x>=[12])*((([15]-[13])/([14]-[12]))*(x-[14])+[15])')
    fun_name='(x<p2)*(((p3-p1)/(p2-p0))*(x-p2)+p3)+(x>=p2&&x<p4)*(((p5-p3)/(p4-p2))*(x-p4)+p5)+(x>=p4&&x<p6)*(((p7-p5)/(p6-p4))*(x-p6)+p7)+(x>=p6&&x<p8)*(((p9-p7)/(p8-p6))*(x-p8)+p9)+(x>=p8&&x<p10)*(((p11-p9)/(p10-p8))*(x-p10)+p11)+(x>=p10&&x<p12)*(((p13-p11)/(p12-p10))*(x-p12)+p13)+(x>=p12)*(((p15-p13)/(p14-p12))*(x-p14)+p15)'
    fun.SetParameter(0,p0)
    fun.SetParameter(1,p1)
    fun.SetParameter(2,p2)
    fun.SetParameter(3,p3)
    fun.SetParameter(4,p4)
    fun.SetParameter(5,p5)
    fun.SetParameter(6,p6)
    fun.SetParameter(7,p7)
    fun.SetParameter(8,p8)
    fun.SetParameter(9,p9)
    fun.SetParameter(10,p10)
    fun.SetParameter(11,p11)
    fun.SetParameter(12,p12)
    fun.SetParameter(13,p13)
    fun.SetParameter(14,p14)
    fun.SetParameter(15,p15)
    return fun

BB_sigma =init_fun_sigma_BB(3.90000e+03,1.39560e+01,-1.73043e-04,1.09326e+00,5.72535e-01,1.40630e-04)
BE_sigma =init_fun_sigma_BB(2.50000e+03,2.49398e+01,5.74472e-02,1.74731e+00,1.91796e+00,-4.17142e-05)
# BB_mean  =init_fun_pol3_pol3(1.11000e+03,9.89793e-01,1.66460e-05,-9.00994e-09,9.96821e-01,4.45426e-07,-1.70423e-10)
BB_mean  =init_fun_interpolation(160, 0.992173, 300, 0.994466, 600, 0.996192, 1100, 0.997227, 1850, 0.997059, 2900, 0.996742, 4000, 0.9957, 4750, 0.995135)
BE_mean  =init_fun_pol1_pol2(6.10000e+02,9.94091e-01,6.56980e-06,9.97773e-01,8.80234e-07,-4.28809e-11)
BB_PowerR=init_fun_pol1_pol2(1.84000e+03,1.08998e+01,2.13291e-02,1.29390e+02,-4.83340e-02,4.64367e-06)
BB_PowerL=init_fun_sigma_BE(-2.30005e+01,3.93681e+02,1.45509e+00)
BB_CutL  =init_fun_pol1_pol2(6.10000e+02,8.79732e-01,1.37066e-03,1.72308e+00,1.88687e-05,3.73202e-09)
BB_CutR  =init_fun_pol3_pol3(2.25000e+03,1.68136e+00,7.51854e-04,-3.23341e-07,4.08748e-02,1.06536e-03,-1.41032e-07)
BE_PowerL=init_fun_PowerL_BE(-6.92166e-03,4.96247e+02,1.40980e+00,7.72180e-05)
BE_PowerR=init_fun_pol1_pol2(6.10000e+02,-2.22560e+01,1.43042e-01,7.66779e+01,-2.56456e-02,2.40496e-06)
BE_CutL  =init_fun_pol3_pol3(1.50000e+03,1.39523e+00,1.13683e-03,-6.08260e-07,1.50300e+00,2.51233e-04,-4.68012e-08)
BE_CutR  =init_fun_sigma_BE(-4.00632e+01,-9.11385e-03,1.87796e+00)


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
dCB.get_value(2000,False) ##True for BB, False for BE
print "mean=%f, sigma=%f, PowerL=%f, PowerR=%f, CutL=%f, CutR=%f"%(dCB.mean,dCB.sigma,dCB.PowerL,dCB.PowerR,dCB.CutL,dCB.CutR)
###### BCD parameter #############    
