/*
 *  Segmentor.h
 *  Segments
 *
 *  Created by Michel Koskas on 06/01/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */

#ifndef _Segmentor_H_
#define _Segmentor_H_

#include <string>
#include <fstream>
#include <iostream>

#include "Function.h"
#include "Sets.h"
#include <math.h>
#include "Constants.h"
#include "Observations.h"
#include <cmath>
#include "MyVector.h"
#include "GeneralFunctions.h"


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
class Segmentor
{
// Fields
private:
  int K;
  MyVector<DataTypeName> y;
  double **C;
  double **Par;
  int **M;
  int n;
	SetTypeName MySet;
  SumOfFunctionsTypeName g;
  FunctionTypeName gamma;


// Methods

	// setter and getter for all of them
public:
  //unsigned int sizeof(SumOfFunctionsTypeName *H, SetTypeName *S, SetTypeName *I, MyVector<int> Candidates, MyVector<int> CandidatesToRemove);
  int **GetM() const;
  int Getn() const;
  int GetK() const;
  double **GetC() const;
  double **GetPar() const;
  Segmentor();
  Segmentor(Observations<DataTypeName>  &yc, int Kc, FunctionTypeName Mg, FunctionTypeName Mgam, SetTypeName *MS);
  ~Segmentor();
  void Initialize(Observations<DataTypeName>  &yc, int Kc, FunctionTypeName Mg, FunctionTypeName Mgam, SetTypeName *MS);
  int ChoosekBreak(double rup);
private:
	// call only by the constructor
  void Initialize();
};


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::~Segmentor()
{
  for (int i = 0; i < K; i++)
    delete[] M[i];
  delete[] M;
  for (int j = 0; j < K; j++)
    delete[] C[j];
  delete[] C;
  for (int j = 0; j < K; j++)
    delete[] Par[j];
  delete[] Par;
	y.clear();
}



template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
int Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::ChoosekBreak(double rup)
{
  int Index = (K-1);
  double FirstSlope = -C[Index-1][n-1]+C[Index][n-1];
  double SecondSlope = -C[Index-2][n-1]+C[Index-1][n-1];


  while ((SecondSlope > (rup*FirstSlope)) && (Index>2))
    {
	Index --;
	FirstSlope = SecondSlope;
	SecondSlope = -C[Index-2][n-1]+C[Index-1][n-1];
    }

  return (Index + 1);
}


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
int **Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::GetM() const
{
  return M;
}

template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
double **Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::GetC() const
{
  return C;
}

template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
double **Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::GetPar() const
{
  return Par;
}

template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
int Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::Getn() const
{
  return n;
}

template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
int Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::GetK() const
{
  return K;
}


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::Segmentor()
{
  K = 0;
  C = NULL;
  Par = NULL;
  M = NULL;
  n = 0;
}


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::Segmentor(Observations<DataTypeName> &yc, int Kc, FunctionTypeName Mg, FunctionTypeName Mgam, SetTypeName *MS)
{
	Initialize(yc, Kc, Mg, Mgam, MS);
}

template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
void Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::Initialize(Observations<DataTypeName> &yc, int Kc, FunctionTypeName Mg, FunctionTypeName Mgam, SetTypeName *MS)
{
  K = Kc;
  n = yc.y.size();
	y = yc.y;
  C = new double *[K];
  for (int i = 0; i < K; i++)
    C[i] = new double[n];
  for (int i = 0; i < K; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = 0;
  Par = new double *[K];
  for (int i = 0; i < K; i++)
    Par[i] = new double[n];
  for (int i = 0; i < K; i++)
    for (int j = 0; j < n; j++)
      Par[i][j] = 0;
  M = new int *[K];
  for (int i = 0; i < K; i++)
    M[i] = new int[n];
  for (int i = 0; i < K; i++)
    for (int j = 0; j < n; j++)
      M[i][j] = 0;
  g = Mg;
  gamma = Mgam;
  MySet.SetMe(*MS);
   Initialize();
}


template <typename SumOfFunctionsTypeName, typename FunctionTypeName, typename AtomicSetElementTypeName, typename DataTypeName, typename SetTypeName>
void Segmentor<SumOfFunctionsTypeName, FunctionTypeName, AtomicSetElementTypeName, DataTypeName, SetTypeName>::Initialize()
{
  SumOfFunctionsTypeName SumOfGammas;
  gamma.SpecializeMe(y[0]);
  SumOfGammas = gamma;
  C[0][0] = (SumOfGammas.Min(MySet));
  Par[0][0] = (SumOfGammas.ArgMin(MySet));
  for (int t = 1; t < n; t++)
  {
    gamma.SpecializeMe(y[t]);
    SumOfGammas += gamma;
    C[0][t] = (SumOfGammas.Min(MySet));
    Par[0][t] = (SumOfGammas.ArgMin(MySet));
  }
  for (int t = 0; t < n; t++)
    M[0][t] = -1;

  // H is 2-dimensional because:
  //    H[k][k][k] is initialized as H[k][k]
  //    and H[k][tau][t] is computed thanks to H[k][tau][t - 1] and replaces it since it is no longer used.
	// In H: first index is the number of break-points, the secon is the value of the breakpoint.
  SumOfFunctionsTypeName *H;   // changed on aug 23rd 2011
  H = new SumOfFunctionsTypeName[n];

  SetTypeName *S = new SetTypeName [n]; // changed on aug 23rd 2011
  // Same remark than for H

  SetTypeName *I = new SetTypeName [n]; // changed on aug 23rd 2011
  // Same again for I



  MyVector<int> *Candidates =new MyVector<int>[K];
  for (int k = 1; k < K; k++)  // k is the number of breakpoints.
  {
    for (int j = 0; j < n; j++)
      S[j].SetMe(MySet);  // changed on aug 23rd 2011
    for (int j = 0; j < n; j++)
      I[j].SetMe(MySet); // Changed on august 23-rd.
    for (int i = 0; i < n; i++)
      H[i].ResetMe();
    Candidates[k].push_back(k-1);
    H[k - 1] = C[k-1][k-1];
    M[k][k] = k - 1;
    for (int t = k; t < n; t++)
    {
      MyVector<int> CandidatesToRemove;
      for (MyVector<int>::iterator Tau = Candidates[k].begin(); Tau != Candidates[k].end(); Tau++)
      {
        gamma.SpecializeMe(y[t]);
        H[*Tau] += gamma;
        SetTypeName *Aux = H[*Tau].IsLowerThan(MySet, C[k-1][t]);
        I[*Tau] = *Aux;
        delete Aux;
        S[*Tau].SelfIntersect(&(I[*Tau]));
        if (S[*Tau].AlmostEmpty())
	    CandidatesToRemove.push_back(*Tau);
      }
      C[k][t] = PLUS_INFINITY; 
      Par[k][t] = MINUS_INFINITY;
      for (MyVector<int>::iterator Tau = Candidates[k].begin(); Tau != Candidates[k].end(); Tau++)
	{
	  S[t].SelfIntersectWithComplementary(&(I[*Tau]));
	  SumOfFunctionsTypeName *Aux = H[*Tau]+g;
	  double mi = Aux->Min(MySet);
	  double par = Aux->ArgMin(MySet);
	  delete Aux;
	  if (mi < C[k][t])
	    {
	      C[k][t]= mi;
	      Par[k][t]= par;
	      M[k][t] = *Tau;
	    }
	}

      for (MyVector<int>::iterator CR = CandidatesToRemove.begin(); CR != CandidatesToRemove.end(); CR++)
        Candidates[k].remove(*CR);
      CandidatesToRemove.clear();
      if (!S[t].AlmostEmpty())
      {
        Candidates[k].push_back(t);
        H[t] += C[k - 1][t];
      }
    }
  }

  delete[] H;
  delete[] S;
  delete[] I;
  for (int i = 0; i < K; i++)
    Candidates[i].clear();
  delete[] Candidates;
	
}


#endif
