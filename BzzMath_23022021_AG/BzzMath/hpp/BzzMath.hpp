// BZZMATH: Release 7.0

// COMMENTARE le tre istruzioni con compilatori diversi da Visual C++ 6
// e poi ripristinare
#ifndef BZZ_COMPILER
#define BZZ_COMPILER 3
#endif
//////////////////////////////////////////////////////////////
// MODIFICARE SOLO PER CREARE LE LIBRERIE
// POI COMMENTARE
//#define BZZ_COMPILER 0 // default: Visual C++ 6 without OPENMP 32 bit
//#define BZZ_COMPILER 1 // Visual C++ 9 (Visual 2008) with OPENMP 32 bit
//#define BZZ_COMPILER 2 // Visual C++ 10 (Visual 2010) with OPENMP 64 bit
//#define BZZ_COMPILER 3 // Visual C++ 11 (Visual 2012) with OPENMP 64 bit
//#define BZZ_COMPILER 11 // INTEL with OPENMP 32 bit
//#define BZZ_COMPILER 101 // g++ LINUX con OPENMP 64 bit
//////////////////////////////////////////////////////////////

#if BZZ_COMPILER == 0
#define BZZ_OPENMP 0
#define BZZ_BLAS 0
#else
#define _CRT_SECURE_NO_DEPRECATE
#define BZZ_OPENMP 1
#define BZZ_BLAS 0
#endif

#if BZZ_COMPILER == 11
#define BZZ_BLAS 0
#endif

#if BZZ_OPENMP == 1
#include <omp.h>
#endif

// Header files for BzzMath library

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include <float.h>

#include "BzzUtility.hpp"

class BzzVectorInt;

class BzzComplex;
class BzzQuadraticEquation;

class BzzMatrixInt;
class BzzMatrixIntBand;
class BzzMatrixCoefficientsExistence;
class BzzMatrixSymmetricCoefficientsExistence;
class BzzMatricesExistence;
class BzzVectorIntArray;

// double precision

class BzzVector;
class BzzVectorSparse;
class BzzMatrix;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;
class BzzMatrixDiagonal;
class BzzMatrixSparse;
class BzzMatrixSparseLockedByRows;
class BzzMatrixSparseLockedByColumns;
class BzzMatrixSparseLocked;
class BzzMatrixSparseSymmetric;
class BzzMatrixBand;

class BzzFactorized;
class BzzFactorizedGauss;
class BzzFactorizedPLR;
class BzzFactorizedSymmetric;
class BzzFactorizedSparseCholesky;
class BzzFactorizedQRLQ;
class BzzFactorizedSVD;
class BzzFactorizedSparseGauss;
class BzzFactorizedBandGauss;
class BzzFactorizedSparseLQ;

class BzzInterpolation;
class BzzVectorArray;
class BzzCubicInterpolation;
class BzzBiCubicInterpolation;

class BzzMatrixTridiagonalBlocks;
class BzzFactorizedTridiagonalBlocksGauss;
class BzzMatrixDiagonalBlocks;
class BzzFactorizedDiagonalBlocksGauss;
class BzzFactorizedFourBlocksGauss;
class BzzMatrixStaircase;
class BzzFactorizedStaircaseGauss;
class BzzMatrixWithNonDecreasingBand;
class BzzMatrixBlocks;

class BzzFactorizedTridiagonalGauss;
class BzzFactorizedSparseLockedByRowsGauss;
class BzzFactorizedBandGaussAndMatrixLocked;
class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked;
class BzzString;
class BzzColorTransformations;

class FaSpaUnd;
class BzzFactorizedSparseSymmetric;
class BzzFactorizedMatrixBlocksGauss;
//class FactBlod;
class FactBlockLQd;
class DiaBloSd;
class LeastCod;
class BzzFactorizedSparseArrayGauss;
class BzzFactorizedSparseLockedGauss;
/*
class Poly5;
*/

class BzzSave;
class BzzLoad;

// Advanced
class BzzFunctionRoot;
class BzzFunctionRootMP;
class BzzFunctionRootRobust;

//Definite integral
class BzzIntegral;

// Mono Dimensional Minimization
class BzzMinimizationMono;
class BzzMinimizationMonoObject;
class BzzMinimizationMonoVeryRobust;
class BzzMinimizationMonoVeryRobustObject;

// Multi Dimensional Minimization
class BzzMinimizationQuasiNewton;
class BzzMinimizationSimplex;
class BzzMinimizationRobust;
class BzzMinimizationRobustObject;
class BzzMinimizationMultiVeryRobust;
class BzzMinimizationMultiVeryRobustObject;

// Numerical Differentiation
class BzzNumericalDifferentiation;

// Linear Regressions
class BzzLinearRegression;
class BzzLinearRegressionExperimentsSearch;

// Non Linear Regression
class BzzNonLinearRegression;

// CSTR Network
class BzzCSTRNetwork;

// Non Linear System
class BzzNonLinearSystem;
class BzzNonLinearSystemSparse;
class BzzNonLinearSystemObject;
class BzzNonLinearSystemSparseObject;

// Ordinary Differential Equations
class BzzOdeMultiValue;
class BzzOdeMultiValueObject;
class BzzOdeRungeKutta;

// Differential-Algebraic Equations
class BzzDae;
class BzzDaeObject;

// BVP
class BzzBVP;
/*

// Linear Programming
class LinProgrd;
*/

///////////////////////////////////////////////////////////////////////////////
#include "BzzVectorInt.hpp"
#include "BzzVector.hpp"

#include "BzzComplex.hpp"
#include "BzzQuadraticEquation.hpp"

#include "BzzMatrixInt.hpp"
#include "BzzMatrixIntBand.hpp"
#include "BzzMatrixCoefficientsExistence.hpp"
#include "BzzMatrixSymmetricCoefficientsExistence.hpp"
#include "BzzMatricesExistence.hpp"
#include "BzzVectorIntArray.hpp"

// double precision

#include "BzzMatrix.hpp"
#include "BzzMatrixLeft.hpp"
#include "BzzMatrixRight.hpp"
#include "BzzMatrixSymmetric.hpp"
#include "BzzMatrixDiagonal.hpp"
#include "BzzMatrixSparse.hpp"
#include "BzzMatrixSparseLockedByRows.hpp"
#include "BzzMatrixSparseLockedByColumns.hpp"
#include "BzzMatrixSparseSymmetric.hpp"
#include "BzzMatrixBand.hpp"

#include "BzzFactorized.hpp"
#include "BzzFactorizedPLR.hpp"
#include "BzzFactorizedSymmetric.hpp"
#include "BzzFactorizedQRLQ.hpp"
#include "BzzFactorizedSVD.hpp"
#include "BzzFactorizedSparseGauss.hpp"
#include "BzzFactorizedBandGauss.hpp"
#include "BzzFactorizedSparseLQ.hpp"

#include "BzzInterpolation.hpp"
#include "BzzCubicInterpolation.hpp"
#include "BzzVectorArray.hpp"
#include "BzzBiCubicInterpolation.hpp"

#include "BzzMatrixTridiagonalBlocks.hpp"
#include "BzzFactorizedTridiagonalBlocksGauss.hpp"
#include "BzzMatrixDiagonalBlocks.hpp"
#include "BzzFactorizedDiagonalBlocksGauss.hpp"
#include "BzzFactorizedFourBlocksGauss.hpp"
#include "BzzMatrixStaircase.hpp"
#include "BzzFactorizedStaircaseGauss.hpp"
#include "BzzMatrixWithNonDecreasingBand.hpp"
#include "BzzMatrixBlocks.hpp"

#include "BzzFactorizedTridiagonalGauss.hpp"
#include "BzzFactorizedSparseLockedByRowsGauss.hpp"
#include "BzzFactorizedBandGaussAndMatrixLocked.hpp"
#include "BzzFactorizedDiagonalBlocksGaussAndMatrixLocked.hpp"
#include "BzzString.hpp"
#include "BzzColorTransformations.hpp"

#include "FaSpaUnd.hpp"
#include "BzzFactorizedSparseSymmetric.hpp"
#include "BzzFactorizedMatrixBlocksGauss.hpp"
//#include "FactBlod.hpp"
#include "FactBlockLQd.hpp"
#include "DiaBloSd.hpp"
#include "LeastCod.hpp"
#include "BzzFactorizedSparseArrayGauss.hpp"
#include "BzzMatrixSparseLocked.hpp"
#include "BzzMatrixCoefficientsLocked.hpp"
#include "BzzMatrixSparseSymmetricLocked.hpp"
#include "BzzFactorizedSparseLockedGauss.hpp"
/*
#include "Poly5.hpp"
*/

#include "BzzSave.hpp"
#include "BzzLoad.hpp"

// Advanced
#include "BzzFunctionRoot.hpp"
#include "BzzFunctionRootRobust.hpp"

//Definite integral
#include "BzzIntegral.hpp"

// Mono Dimensional Minimization
#include "BzzMinimizationMono.hpp"
#include "BzzMinimizationMonoObject.hpp"
#include "BzzMinimizationMonoVeryRobustObject.hpp"
#include "BzzMinimizationMonoVeryRobust.hpp"

// Multi Dimensional Minimization
#include "BzzMinimizationTwoVeryRobust.hpp"
#include "BzzOptnov.hpp"
#include "BzzMinimizationQuasiNewton.hpp"
#include "BzzMinimizationLargeSparseNewton.hpp"
#include "BzzQuadraticProgramming.hpp"
#include "BzzMinimizationSimplex.hpp"
#include "BzzMinimizationRobust.hpp"
#include "BzzMinimizationRobustObject.hpp"
#include "BzzMinimizationMultiVeryRobust.hpp"
#include "BzzMinimizationMultiVeryRobustObject.hpp"
// Numerical Differentiation
#include "BzzNumericalDifferentiation.hpp"

// Linear Regressions
#include "BzzLinearRegression.hpp"
#include "BzzLinearRegressionExperimentsSearch.hpp"

// Non Linear Regression
#include "BzzNonLinearRegression.hpp"

// CSTR Network
#include "BzzCSTRNetwork.hpp"

// Non Linear System

#include "BzzNonLinearSystemUtilities.hpp"
#include "BzzNonLinearSystem.hpp"
#include "BzzNonLinearSystemSparse.hpp"
#include "BzzNonLinearSystemObject.hpp"
#include "BzzNonLinearSystemSparseObject.hpp"

// Ordinary Differential Equations
#include "BzzOdeMultiValue.hpp"
#include "BzzOdeMultiValueObject.hpp"
#include "BzzOdeRungeKutta.hpp"

// Differential-Algebraic Equations
#include "BzzDae.hpp"
#include "BzzDaeObject.hpp"

// BVP
#include "BzzBVP.hpp"

// Constrained Minimization
//#include "BzzLinearProgrammingAttic.hpp"
//#include "BzzQuadraticProgramming.hpp"

#include "BzzConstrainedMinimization.hpp"

/*

// Linear Programming
#include "LinProgrd.hpp"
*/