/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
o *****************************************************************************/

#ifndef ZPRIMEMUONBKGPDFTRIGDOWN
#define ZPRIMEMUONBKGPDFTRIGDOWN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class ZPrimeMuonBkgPdfTrigDown : public RooAbsPdf {
public:
  ZPrimeMuonBkgPdfTrigDown() {} ; 
  ZPrimeMuonBkgPdfTrigDown(const char *name, const char *title,
              RooAbsReal& _mass,
              RooAbsReal& _bkg_a,
              RooAbsReal& _bkg_b,
              RooAbsReal& _bkg_c,
              RooAbsReal& _bkg_d,
              RooAbsReal& _bkg_e,
              RooAbsReal& _bkg_syst_a,
              RooAbsReal& _bkg_syst_b,
	      RooAbsReal& _cat);
  ZPrimeMuonBkgPdfTrigDown(const ZPrimeMuonBkgPdfTrigDown& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ZPrimeMuonBkgPdfTrigDown(*this,newname); }
  inline virtual ~ZPrimeMuonBkgPdfTrigDown() { }

protected:

  RooRealProxy mass ;
  RooRealProxy bkg_a ;
  RooRealProxy  bkg_b ;
  RooRealProxy  bkg_c ;
  RooRealProxy  bkg_d ;
  RooRealProxy  bkg_e ;
  RooRealProxy  bkg_syst_a ;
  RooRealProxy  bkg_syst_b ;
  RooRealProxy  cat ;
  
  Double_t evaluate() const ;

private:

  ClassDef(ZPrimeMuonBkgPdfTrigDown,1) // Your description goes here...
};
 
#endif