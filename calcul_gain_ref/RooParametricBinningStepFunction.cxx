/*****************************************************************************
 * Authors: Xalbat Aguerre and Loic Labit
 * Modification of RooParametricStepFunction to allow modification of the
 * binning

 * Project: RooFit                                                           *
 * Package: RooFitBabar                                                      *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *    Aaron Roodman, Stanford Linear Accelerator Center, Stanford University *
 *                                                                           *
 * Copyright (c) 2004, Stanford University. All rights reserved.        *
 *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooParametricBinningStepFunction
  \ingroup Roofit

  The Parametric Step Function PDF is a binned distribution whose parameters
  are the heights of each bin.  This PDF was first used in BaBar's B0->pi0pi0
  paper BABAR Collaboration (B. Aubert et al.) Phys.Rev.Lett.91:241801,2003.

  This PDF may be used to describe oddly shaped distributions.  It differs
  from a RooKeysPdf or a RooHistPdf in that a RooParametricBinningStepFunction
  has free parameters.  In particular, any statistical uncertainty in
  sample used to model this PDF may be understood with these free parameters;
  this is not possible with non-parametric PDFs.

  The RooParametricBinningStepFunction has Nbins-1 free parameters. Note that
  the limits of the dependent variable must match the low and hi bin limits.

 */

#include "Riostream.h"
#include "TArrayD.h"
#include <math.h>

#include "RooParametricBinningStepFunction.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#include "TError.h"

using namespace std;

ClassImp(RooParametricBinningStepFunction);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooParametricBinningStepFunction::RooParametricBinningStepFunction(const char* name, const char* title,
		RooAbsReal& x, const RooArgList& coefList, RooArgList& limits, Int_t nBins) :
	RooAbsPdf(name, title),
	_x("x", "Dependent", this, x),
	_coefList("coefList","List of coefficients",this),
	_limits("limitList","List of bin limits",this),
	_nBins(nBins)
{

	// Check lowest order
	if (_nBins<0) {
		cout << "RooParametricBinningStepFunction::ctor(" << GetName()
			<< ") WARNING: nBins must be >=0, setting value to 0" << endl ;
		_nBins=0 ;
	}
//To make it a little more foolproof

  // Check that lengths of arrays are consistent
  if (_nBins +1 != limits.getSize()) {
    cout << "RooParametricBinningStepFunction::ctor(" << GetName()
    << ") ERROR: "<< limits.GetName()<<" RooArgList must have one more element than the number of bins nBins >=0" << endl ;
      R__ASSERT(0) ;
  }
  if (_nBins != coefList.getSize() ) {
    cout << "RSTEtrueParamPdf::ctor(" << GetName()
    << ") ERROR: coefficients list " << coefList.GetName() << " must have same number of elements as nBins" << endl ;
      R__ASSERT(0) ;
  }

	for (auto *coef : coefList) {
		if (!dynamic_cast<RooAbsReal*>(coef)) {
			cout << "RooParametricBinningStepFunction::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
				<< " is not of type RooAbsReal" << endl ;
			R__ASSERT(0) ;
		}
		_coefList.add(*coef) ;
	}

	// Bin limits
	//  limits.Copy(_limits);
	for (auto *lim : limits) {
		if (!dynamic_cast<RooAbsReal*>(lim)) {
			cout << "RooParametricBinningStepFunction::ctor(" << GetName() << ") ERROR: limit " << lim->GetName()
				<< " is not of type RooAbsReal" << endl ;
			R__ASSERT(0) ;
		}
		_limits.add(*lim) ;
	}


}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooParametricBinningStepFunction::RooParametricBinningStepFunction(const RooParametricBinningStepFunction& other, const char* name) :
	RooAbsPdf(other, name),
	_x("x", this, other._x),
	_coefList("coefList",this,other._coefList),
	_limits("limitList",this,other._limits),
	_nBins(other._nBins)
{
	//  (other._limits).Copy(_limits);
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooParametricBinningStepFunction::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
	if (matchArgs(allVars, analVars, _x)) return 1;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////

double RooParametricBinningStepFunction::analyticalIntegral(Int_t code, const char* rangeName) const
{
	R__ASSERT(code==1) ;

	// Case without range is trivial: p.d.f is by construction normalized
	if (!rangeName) {
		double sum = 0;
		for (int i = 0; i < _coefList.getSize(); i++) {
			sum+=((RooRealVar*)_coefList.at(i))->getVal()*(((RooRealVar*)_limits.at(i+1))->getVal()-((RooRealVar*)_limits.at(i))->getVal());
		}
		return sum;
	}

	// Case with ranges, calculate integral explicitly
	double xmin = _x.min(rangeName) ;
	double xmax = _x.max(rangeName) ;

	double sum=0 ;
	Int_t i ;
	for (i=1 ; i<=_nBins ; i++) {
		double binVal = (i<_nBins) ? (static_cast<RooRealVar*>(_coefList.at(i-1))->getVal()) : lastBinValue() ;
		if (((RooRealVar*)_limits.at(i-1))->getVal()>=xmin && ((RooRealVar*)_limits.at(i))->getVal()<=xmax) {
			// Bin fully in the integration domain
			sum += (((RooRealVar*)_limits.at(i))->getVal()-((RooRealVar*)_limits.at(i-1))->getVal())*binVal ;
		} else if (((RooRealVar*)_limits.at(i-1))->getVal()<xmin && ((RooRealVar*)_limits.at(i))->getVal()>xmax) {
			// Domain is fully contained in this bin
			sum += (xmax-xmin)*binVal ;
			// Exit here, this is the last bin to be processed by construction
			return sum ;
		} else if (((RooRealVar*)_limits.at(i-1))->getVal()<xmin && ((RooRealVar*)_limits.at(i))->getVal()<=xmax && ((RooRealVar*)_limits.at(i))->getVal()>xmin) {
			// Lower domain boundary is in bin
			sum +=  (((RooRealVar*)_limits.at(i))->getVal()-xmin)*binVal ;
		} else if (((RooRealVar*)_limits.at(i-1))->getVal()>=xmin && ((RooRealVar*)_limits.at(i))->getVal()>xmax && ((RooRealVar*)_limits.at(i-1))->getVal()<xmax) {
			sum +=  (xmax-((RooRealVar*)_limits.at(i-1))->getVal())*binVal ;
			// Upper domain boundary is in bin
			// Exit here, this is the last bin to be processed by construction
			return sum ;
		}
	}

	return sum;
}

////////////////////////////////////////////////////////////////////////////////

double RooParametricBinningStepFunction::lastBinValue() const
{
	double sum(0.);
	double binSize(0.);
	for (Int_t j=1;j<_nBins;j++){
		RooRealVar* tmp = (RooRealVar*) _coefList.at(j-1);
		binSize = ((RooRealVar*)_limits.at(j))->getVal() - ((RooRealVar*)_limits.at(j-1))->getVal();
		sum = sum + tmp->getVal()*binSize;
	}
	binSize = ((RooRealVar*)_limits.at(_nBins))->getVal() - ((RooRealVar*)_limits.at(_nBins-1))->getVal();
	return (1.0 - sum)/binSize;
}

////////////////////////////////////////////////////////////////////////////////

double RooParametricBinningStepFunction::evaluate() const
{
	double value(0.);
	if (_x >= ((RooRealVar*)_limits.at(0))->getVal() && _x < ((RooRealVar*)_limits.at(_nBins))->getVal()){

		for (Int_t i=1;i<=_nBins;i++){
			if (_x < ((RooRealVar*)_limits.at(i))->getVal()){
				// in Bin i-1 (starting with Bin 0)
//				cout<<"LABIT::DEBUG "<<"_x="<<_x<<" _limits["<<i<<"]="<<((RooRealVar*)_limits.at(i))->getVal()<<endl;
				if (i<_nBins) {

					// not in last Bin
					RooRealVar* tmp = (RooRealVar*) _coefList.at(i-1);
					value =  tmp->getVal();
//					cout<<"LABIT::DEBUG "<<"_i-1="<<i-1<<" __coefList["<<i-1<<"]="<<value<<endl;
					break;
				} else {
					// in last Bin
					double sum(0.);
					double binSize(0.);
					for (Int_t j=1;j<_nBins;j++){
						RooRealVar* tmp = (RooRealVar*) _coefList.at(j-1);
						binSize = ((RooRealVar*)_limits.at(j))->getVal() - ((RooRealVar*)_limits.at(j-1))->getVal();
						sum = sum + tmp->getVal()*binSize;
					}
					binSize = ((RooRealVar*)_limits.at(_nBins))->getVal() - ((RooRealVar*)_limits.at(_nBins-1))->getVal();
					value = (1.0 - sum)/binSize;
					if (value<=0.0){
						value = 0.000001;
						//       cout << "RooParametricBinningStepFunction: sum of values gt 1.0 -- beware!!" <<endl;
					}
					break;
				}
			}
		}

	}
	return value;

}

////////////////////////////////////////////////////////////////////////////////

Int_t RooParametricBinningStepFunction::getnBins(){
	return _nBins;
}

////////////////////////////////////////////////////////////////////////////////

double* RooParametricBinningStepFunction::getLimits(){

	TArrayD* Arraylimits = new TArrayD(_limits.getSize());
	for (size_t i = 0; i < _limits.getSize(); i++) {
		Arraylimits->SetAt(((RooRealVar*)_limits.at(i))->getVal(), i);
	}
	double* limoutput = Arraylimits->GetArray();
	return limoutput;
}
