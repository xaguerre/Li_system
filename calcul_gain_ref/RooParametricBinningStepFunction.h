/*****************************************************************************
* Authors: Xalbat Aguerre and Loic Labit
* Modification of RooParametricStepFunction to allow modification of the
* binning
*
* Project: RooFit                                                           *
* Package: RooFitModels                                                     *
*    File: $Id: RooParametricBinningStepFunction.h,v 1.5 2007/05/11 09:13:07 verkerke Exp $
* Authors:                                                                  *
*    Aaron Roodman, Stanford Linear Accelerator Center, Stanford University *
*                                                                           *
* Copyright (c) 2000-2005, Stanford University. All rights reserved.        *
*
*                                                                           *
* Redistribution and use in source and binary forms,                        *
* with or without modification, are permitted according to the terms        *
* listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
*****************************************************************************/
#ifndef ROO_PARAMETRIC_BINNING_STEP_FUNCTION
#define ROO_PARAMETRIC_BINNING_STEP_FUNCTION

#include "TArrayD.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooParametricBinningStepFunction : public RooAbsPdf {
public:

  RooParametricBinningStepFunction() : _nBins(0), _coefIter(0) {}

  RooParametricBinningStepFunction(const char *name, const char *title,
    RooAbsReal& x, const RooArgList& coefList, RooArgList& limits, Int_t nBins=1) ;

    RooParametricBinningStepFunction(const RooParametricBinningStepFunction& other, const char* name = 0);
    TObject* clone(const char* newname) const override { return new RooParametricBinningStepFunction(*this, newname); }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
    double analyticalIntegral(Int_t code, const char* rangeName=0) const override ;
    Int_t getnBins();
    double* getLimits();

  protected:

    double lastBinValue() const ;

    RooRealProxy _x;
    RooListProxy _coefList ;
    RooListProxy _limits;
    Int_t _nBins ;
    TIterator* _coefIter ;  //! do not persist

    double evaluate() const override;

    ClassDefOverride(RooParametricBinningStepFunction,1) // Parametric Step Function Pdf
  };

  #endif
