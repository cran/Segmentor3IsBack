#include "CallSegmentor.h"
#include "GeneralFunctionsDeclarations.h"
#include "R.h"
#include "Rmath.h"

extern "C"
{

  void SegmentPoisson(int *Size, int *KMax, int *Data, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorPoisson(Size, KMax, Data, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentBinNeg(int *Size, int *KMax, double *theta, int *Data, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorBinNeg(Size, KMax, theta, Data, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentNormal(int *Size, int *KMax, double *Data, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorNormal(Size, KMax, Data, Breakpoints, Parameters, Likelihood);
    return;
  }

  void SegmentVariance(int *Size, int *KMax, double *mu, double *Data, int *Breakpoints, double *Parameters, double *Likelihood)
  {
    CallSegmentorVariance(Size, KMax, mu, Data, Breakpoints, Parameters, Likelihood);
    return;
  }
}
