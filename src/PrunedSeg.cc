#include "CallSegmentor.h"
#include "GeneralFunctionsDeclarations.h"
#include "R.h"
#include "Rmath.h"

extern "C"
{

  void SegmentPoisson(int *Size, int *KMax, int *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorPoisson(Size, KMax, Data, DataComp, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentBinNeg(int *Size, int *KMax, double *theta, int *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorBinNeg(Size, KMax, theta, Data, DataComp, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentNormal(int *Size, int *KMax, double *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorNormal(Size, KMax, Data, DataComp, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentVariance(int *Size, int *KMax, double *mu, double *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorVariance(Size, KMax, mu, Data, DataComp, Breakpoints, Parameters, Likelihood);
    return;
  }
  
  void SegmentPoissonKeep(int *Size, int *KMax, int *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost, int *Pos)
  {
    CallSegmentorPoissonKeep(Size, KMax, Data, DataComp, Breakpoints, Parameters, Likelihood, Cost, Pos);
    return;
  }

  void SegmentBinNegKeep(int *Size, int *KMax, double *theta, int *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost, int *Pos)
  {
    CallSegmentorBinNegKeep(Size, KMax, theta, Data, DataComp, Breakpoints, Parameters, Likelihood, Cost, Pos);
    return;
  }

  void SegmentNormalKeep(int *Size, int *KMax, double *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost, int *Pos)
  {
    CallSegmentorNormalKeep(Size, KMax, Data, DataComp, Breakpoints, Parameters, Likelihood, Cost, Pos);
    return;
  }

  void SegmentVarianceKeep(int *Size, int *KMax, double *mu, double *Data, int *DataComp, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost, int *Pos)
  {
    CallSegmentorVarianceKeep(Size, KMax, mu, Data, DataComp, Breakpoints, Parameters, Likelihood, Cost, Pos);
    return;
  }
}
