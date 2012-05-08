
#include "BinNegative.h"
#include "Segment.h"
#include "Constants.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cmath>
#include "GeneralFunctionsDeclarations.h"



  // Negative Binomial(mu)= A - Y.ln(mu) - Theta.ln(1-mu)
  // New parametrisation: A - S.ln(mu) - T.ln(1-mu)

BinNegative::BinNegative()
{
  ResetMe();
}


BinNegative::BinNegative(double t, int Y)
{
  ResetMe(t, Y);
}

BinNegative::BinNegative(double a, double s, double t)
{
  ResetMe(a, s, t);
}

void BinNegative::ResetMe()
{
  A = 0;
  T = 0;
  S = 0;
  FirstElement = 0;
  FirstElementSpecified = true;
}


void BinNegative::ResetMe(double t, int Y)
{
  A = 0;
  T = t;
  S = Y;
  FirstElement = Y;
  FirstElementSpecified = true;
  for(double i=0; i<Y; i++)
     A+=log(i+1)-log(T+i);
}

void BinNegative::ResetMe(double a, double s, double t)
{
  A=a;
  T = t;
  S = s;
  FirstElement = s;
  FirstElementSpecified = true;
}

BinNegative *BinNegative::operator+(BinNegative &Other)
{
  BinNegative *Res = new BinNegative;
  Res->A = (*this).A + Other.A;
  Res->S = (*this).S + Other.S;
  Res->T = (*this).T + Other.T;
  return Res;
}

BinNegative *BinNegative::operator+(const double &C)
{
  BinNegative *Res = new BinNegative;
  Res->A = (*this).A + C ;
  Res->S = (*this).S;
  Res->T = (*this).T;
  (*Res).FirstElementSpecified = true;
  return Res;
}

void BinNegative::operator+=(const double &C)
{
  A += C;
  FirstElementSpecified = true;
}

void BinNegative::operator+=(BinNegative &Other)
{
  A += Other.A;
  S += Other.S;
  T += Other.T;
  FirstElementSpecified = true;
}


void BinNegative::SpecializeMe(int Y)
{
  A = 0;
  S = Y;
  FirstElementSpecified = true;
  if (Y!=0)
    for(int i=0; i<Y; i++)
      A+=log((double) (i+1))-log((double) (T+i));

}


double BinNegative::operator()(int y, double mu)
{
  double Ap = 0;
  if (y!=0)
    for(double i = 0; i < y; i++)
      Ap += log(i+1)-log(T+i);
  if (mu == 0)
    {
       if (y == 0)
         return Ap;
     else
         return PLUS_INFINITY;
    }
  else
    {
      if (mu ==1)
      {
	if (T == 0)
	  return Ap;
	else
	  return PLUS_INFINITY;
      }
      else
        return Ap - y * log(mu) - T * log(1-mu);
    }
}


double BinNegative::operator()(double mu)
{
    if (mu == 0)
    {
       if (S == 0)
         return A;
     else
         return PLUS_INFINITY;
    }
    else
    {
      if (mu ==1)
      {
	if (T == 0)
	  return A;
	else
	  return PLUS_INFINITY;
      }
      else
        return A - S * log(mu) - T * log(1-mu);
    }
}

double BinNegative::operator[](double mu)
{
    if (mu !=0)
      return (- S / mu + T /(1-mu) );
    return 0;
}

double BinNegative::Min(Segment &LS)
{
  double xmin = S/(S+T);
  return (*this)(xmin);
}

double BinNegative::ArgMin(Segment &LS)
{
  return S/(S+T);
}

double BinNegative::Min()
{
  double xmin = S/(S+T);
  return (*this)(xmin);
}

double BinNegative::ArgMin()
{
  return S/(S+T);
}


double BinNegative::Min(MultiSegment &MS)
{
  if (MS.Empty())
    return PLUS_INFINITY;
  double Answer = PLUS_INFINITY;
  for (MyVector<Segment>::iterator I = MS.GetMySegments().begin(); I != MS.GetMySegments().end(); I++)
  {
    Answer = std::min(Answer, (*this).Min(*I));
  }
  return Answer;
}

double BinNegative::ArgMin(MultiSegment &MS)
{
  if (MS.Empty())
    return PLUS_INFINITY;
  double Answer = PLUS_INFINITY;
  double min = PLUS_INFINITY;
  for (MyVector<Segment>::iterator I = MS.GetMySegments().begin(); I != MS.GetMySegments().end(); I++)
	if ((*this).Min(*I) < min)
	{
	  Answer = (*this).ArgMin(*I);
	  min = (*this).Min(*I);
	}

  return Answer;
}


MultiSegment *BinNegative::LowerThanZero(MultiSegment &MS)
{
  Segment I(MINUS_INFINITY,PLUS_INFINITY);
  if ((*this).S==0)
  {
    if ((*this).T==0)
    {
	if ((*this).A<=0)
		I.SetMe(MINUS_INFINITY,PLUS_INFINITY,false,false);
	else
		I.SetMe(PLUS_INFINITY,MINUS_INFINITY,false,false);
    }
    else
    {
	double R = 1 - exp((*this).A/(*this).T);
	I.SetMe(MINUS_INFINITY,R,false,true);
    }
  }
  else
  {
    if ((*this).T==0)
    {
	double R = exp((*this).A/(*this).S);
	I.SetMe(R, PLUS_INFINITY, true, false);
    }
    else
    {
	double xmin = S/(S+T);
	double TheMin = (*this)(xmin);
	if (abs(TheMin) < EPSILON)
		I.SetMe(xmin,xmin,true,true);
	else if (TheMin > 0)
		I.SetMe(PLUS_INFINITY,MINUS_INFINITY,false,false);
	else
	{
		double FirstRoot, SecondRoot;
		// Now computing the first root (near 0)
		double V = xmin, U;
		while ((*this)(V) < 0)
			V /= 2;
		U = 2 * V;
		while (abs(V - U) > EPSILON)
		{
			U = V;
			V = U - ((*this)(U)) / (T/(1-U) - S / U);
		}
		FirstRoot = V;

		// Now computing the second root (near 1)
		V = xmin;
		while ((*this)(V) < 0)
			V = (1 + V) / 2;
		U = 2 * V - 1;
		while (abs(V - U) > EPSILON)
		{
			U = V;
			V = U - ((*this)(U)) / (T/(1-U) - S / U);
		}
		SecondRoot = V;
		I.SetMe(FirstRoot, SecondRoot,true, true);
	}
    }
  }
  return MS.Intersect(I);

}


MultiSegment *BinNegative::IsLowerThan(MultiSegment &S, double C)
{
  A = A - C;
  MultiSegment *Answer = LowerThanZero(S);
  A = A + C;
  return Answer;
}




