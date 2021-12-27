//#define WINDOWS 1  // 1 WINDOWS VS6
//#define WINDOWS 0  // 0 LINUX, INTEL 10.1

// BZZMATH: Release 7.0

// =============================< UTILITY.HPP >================================
// * Header File for UTILITY.CPP                                              *
// * Description: Dal Fortan al C++ (Capitolo 8)                              *
// *              by G. Buzzi-Ferraris                                        *
// *              Addison Wesley(1991)                                        *
// *                    and                                                   *
// *              Scientific C++ (Chapter 8)                                  *
// *              by G. Buzzi-Ferraris                                        *
// *              Addison Wesley(1993)                                        *
// *                                                                          *
// * Examples: c:\bzzmath\examples\BzzMathBaisc											*
// *					\utilities\utilities\utilities.cpp									*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
// 04-1991  Date Written
// 11-1992  English version
// 01-1994  Modified Max, Min, MaxAbs, MinAbs.
// 03-1994  Added BzzBaseClass for BzzPrint and BzzMessage.
// 04-1994  Added BzzMessageAndGetchar
// 02-1995  Added the Bzz prefix to the names of the classes.
// 03-1995  Added #define BZZ_LONG_DOUBLE and BZZ_INT.
// 04-1995  Added BzzPause function.
// 10-1995  Added BzzClock and BzzDiffClock functions.
// 03-1996  Added template Middle function.
// 03-1996  Abs, Max, Min, MaxAbs, MinAbs modified in template.
// 03-1996  Sort and Swap modified in template.
// 04-1996  Modified BzzError.
// 04-1997  Added BzzWarning.
// 06-1997  Added BzzGetUserTime, BzzGetKernelTime, BzzGetCpuTime.
// 10-1998  Added BzzSetFloatingPointExceptions.
// 10-1998  Added BzzPow2, BzzPow3, ..., BzzPow10.

//////////	Version 4.0
//	01-2000	Bug fixed in BzzGetCpuTime.
//	01-2000	Added RelativeEfficience.
//	02-2000	Added functions for Statistics: BzzLogGamma,BzzGamma,BzzLogFactorial,
//					BzzBinomial,BzzBeta,BzzIncompleteGammaP,BzzIncompleteGammaQ,
//					BzzErf,BzzErfc,BzzIncompleteBeta.
//	03-2000	Added BzzBiasedFactorForStandardDeviation.
//	03-2000	Added BzzOneWayTTest, BzzOneWayFTest, BzzTwoWayFTest.
//	02-2001	Added BzzBinomialProbabilityDistribution.
//	02-2001	Added BzzChiSquareLowerThanFixed and BzzChiSquareGreaterThanFixed.
//	02-2001	Added BzzPoissonProbabilityDistribution.
// 12-2002	Added inline BzzPow2,3,4,5,6,7,8,9,10.

//////////	Version 5.0
//	09-2003	Added RemoveWarningWindow, RemoveWarningWindowAndWarningPrint functions.
//	09-2003	Added SetWarningWindow function.
//	10-2004	Added GetBzzVersionMessage,GetBzzVersionMajor,GetBzzVersionMinor,GetBzzVersionRevision.
//	07-2010	Added RemoveErrorWindow function.
//	07-2010	Added SetErrorWindow function.

//////////	Version 6.0
//	01-2012	Added BzzFibonacci functions.

#ifndef BZZ_UTILITY_HPP
#define BZZ_UTILITY_HPP

#include <stdio.h>
#include <stdlib.h>
#include <math.h> // for inline functions
#include <time.h>
//#include <eh.h>

extern int bzzMajor;
extern int bzzMinor;
extern int bzzReview;

class BzzVector;

void BzzSetFloatingPointExceptions(int exception = 0);
void BzzMessage(const char* myFormat, ...);
double RelativeEfficience(void);

class BzzExceptionHandling
{
	int startEx, endEx;
public:
	BzzExceptionHandling(int sEx, int eEx)
	{
		startEx = sEx;
		endEx = eEx;
		::BzzSetFloatingPointExceptions(startEx);
	}
	BzzExceptionHandling(void)
	{
		startEx = 0;
		endEx = 0;
		::BzzSetFloatingPointExceptions(0);
	}
	~BzzExceptionHandling(void)
	{
		::BzzSetFloatingPointExceptions(endEx);
	}
};

class BzzPrintClockClass;
class BzzClockClass;

// ============================================================================
// ==============================   GLOBAL   ==================================
// ============================================================================
//extern BzzPrintClockClass BzzPrintClock;
//extern BzzClockClass BzzClock;
enum SplineType
{
	NATURAL,
	ASSIGNED_SECOND_DERIVATIVE,
	SECOND_DERIVATIVE_ADJACENT,
	ASSIGNED_FIRST_DERIVATIVE,
	PERIODIC
};

extern BzzExceptionHandling	bzzExceptionHandling;
extern FILE* bzzFileOut;
extern FILE* bzzFileIn;
extern int bzzOpenMP;

extern char bzzMessageYesNo;
extern char bzzStop;
extern char bzzUnfeasible;
extern char bzzWarningWindow;
extern char bzzErrorWindow;

extern char bzzSqrtError;
extern char bzzLogError;
extern char bzzLog10Error;
extern char bzzOverflowError;
inline void BzzSetSqrtAbsolute(void) { bzzSqrtError = 1; }
inline void BzzSetLogAbsolute(void) { bzzLogError = 1; }
inline void BzzSetLog10Absolute(void) { bzzLog10Error = 1; }
inline void BzzSetOverflowBigValue(void) { bzzOverflowError = 1; }

//extern char bzzE;

extern const float BZZ_MACH_EPS_FLOAT;
extern const double BZZ_MACH_EPS_DOUBLE;

extern const double BZZ_PI_GRECO;

extern const double BZZ_RADIX;
extern const double BZZ_BIG;
extern const float BZZ_BIG_FLOAT;
extern const float BZZ_TINY_FLOAT;
extern const double BZZ_BIG_DOUBLE;
extern const double BZZ_TINY_DOUBLE;
extern const long double BZZ_BIG_LONG_DOUBLE;
extern const long double BZZ_TINY_LONG_DOUBLE;
extern const char BZZ_BIG_CHAR;
extern const char BZZ_TINY_CHAR;
extern const unsigned char BZZ_BIG_UNSIGNED_CHAR;
extern const unsigned char BZZ_TINY_UNSIGNED_CHAR;
extern const int BZZ_BIG_INT;
extern const int BZZ_TINY_INT;
extern const unsigned int BZZ_BIG_UNSIGNED_INT;
extern const unsigned int BZZ_TINY_UNSIGNED_INT;
extern const long int BZZ_BIG_LONG;
extern const signed long int BZZ_TINY_LONG;
extern const unsigned long int BZZ_BIG_UNSIGNED_LONG;
extern const unsigned long int BZZ_TINY_UNSIGNED_LONG;

// ============================================================================
// ===============================  ERRORS   ==================================
// ============================================================================

//========== Error Type ==============
extern const char* const BZZ_ERR_RANGE;
extern const char* const BZZ_ERR_SPACE;
extern const char* const BZZ_ERR_OPEN_FILE;
extern const char* const BZZ_ERR_CLOSE_FILE;
extern const char* const BZZ_ERR_READING_FILE;
extern const char* const BZZ_ERR_WRITING_FILE;
extern const char* const BZZ_ERR_CHECK_DIMENSION;
extern const char* const BZZ_ERR_OUT_OF_SCOPE;
extern const char* const BZZ_ERR_FACTORIZED;
extern const char* const BZZ_ERR_IMPLEMENTATION;
extern const char* const BZZ_ERR_E;

//============ Function Type =============
extern const char* const BZZ_ERR_FUNCTION;
extern const char* const BZZ_ERR_CONSTRUCTOR;
extern const char* const BZZ_ERR_OPERATOR;

//============ Version Type =============
extern const char* const BZZ_BZZMATH_VERSION;

// ============================================================================
// ===========================< Template Functions >===========================
// ============================================================================
// *******************************< Swap >*************************************
// * Purpose: Having a unique function for Swap                               *
// * Example: Swap(&x,&y);                                                    *
// ****************************************************************************
template <class T>
void Swap(T* x, T* y)
{
	T temp = *x;
	*x = *y;
	*y = temp;
}

// *********************************< Abs >************************************
// * Purpose: To have a unique function for absolute value                    *
// * Description: Overload for float, int, double                             *
// * Example: float x = Abs(y);                                               *
// ****************************************************************************
template <class T>
inline T Abs(T a)
{
	return a < 0 ? -a : a;
}

// *********************< Max, Min, MaxAbs, MinAbs >***************************
// * Purpose: Maximum and minimum for a pair of numbers                       *
// * Description: Functions overload float, int, double                       *
// * Example: x = Max(a,b);i = Min(l,k); x = MaxAbs(a,b);                     *
// ****************************************************************************
template <class T>
inline T Max(T a, T b)
{
	return (a > b ? a : b);
}
template <class T>
inline T Min(T a, T b)
{
	return (a < b ? a : b);
}
template <class T>
inline T MaxAbs(T a, T b)
{
	return (Abs(a) > Abs(b) ? Abs(a) : Abs(b));
}
template <class T>
inline T MinAbs(T a, T b)
{
	return (Abs(a) < Abs(b) ? Abs(a) : Abs(b));
}

// *********************< Max, Min, MaxAbs, MinAbs >***************************
// * Purpose: Max, Min, MaxAbs, MinAbs of an array                            *
// * Description: It returns the Max of an array.                             *
// *              If the index im is put in the argument                      *
// *              it also returns the position.                               *
// * Example: y = Max(10,x,&im);                                              *
// ****************************************************************************
template <class T>
T Max(int n, T* x, int* im)
{
	if (n < 0) return x[0];
	T temp;
	temp = x[0];
	if (im != 0) *im = 0;
	for (int i = 1;i < n;i++)
		if (temp < x[i])
		{
			temp = x[i];
			if (im != 0)
				*im = i;
		}
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T Max(int n, T* x)
{
	if (n < 0) return x[0];
	T temp;
	temp = x[0];
	for (int i = 1;i < n;i++)
		if (temp < x[i])
			temp = x[i];
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T MaxAbs(int n, T* x, int* im)
{
	T temp = Abs(x[0]);
	if (n < 0) return temp;
	if (im != 0) *im = 0;
	for (int i = 1;i < n;i++)
		if (temp < Abs(x[i]))
		{
			temp = Abs(x[i]);
			if (im != 0)
				*im = i;
		}
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T MaxAbs(int n, T* x)
{
	T temp = Abs(x[0]);
	if (n < 0) return temp;
	for (int i = 1;i < n;i++)
		if (temp < Abs(x[i]))
			temp = Abs(x[i]);
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T Min(int n, T* x, int* im)
{
	if (n < 0) return x[0];
	T temp = x[0];
	if (im != 0) *im = 0;
	for (int i = 1;i < n;i++)
		if (temp > x[i])
		{
			temp = x[i];
			if (im != 0)
				*im = i;
		}
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T Min(int n, T* x)
{
	if (n < 0) return x[0];
	T temp = x[0];
	for (int i = 1;i < n;i++)
		if (temp > x[i])
			temp = x[i];
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T MinAbs(int n, T* x, int* im)
{
	T temp = Abs(x[0]);
	if (n < 0) return temp;
	if (im != 0) *im = 0;
	for (int i = 1;i < n;i++)
		if (temp > Abs(x[i]))
		{
			temp = Abs(x[i]);
			if (im != 0)
				*im = i;
		}
	return temp;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T MinAbs(int n, T* x)
{
	T temp = Abs(x[0]);
	if (n < 0) return temp;
	for (int i = 1;i < n;i++)
		if (temp > Abs(x[i]))
			temp = Abs(x[i]);
	return temp;
}

// *****************************< Middle >*************************************
// * Purpose: Retriving the middle of three objects                           *
// * Description: Given three objects a, b, c this function returns           *
// * midle object. It is necessary to overload the operator <=                *
// * Example: x = Middle(a,b,c);                                              *
// ****************************************************************************
template <class T>
inline T Middle(T a, T b, T c)
{
	return(a <= b ? (b <= c ? b : Middle(a, c, b)) : Middle(b, a, c));
}

// *******************************< Sort >*************************************
// * Purpose: To sort an array                                                *
// * Description: Chapter 5 Scientific C++.                                   *
// * It require the functions BzzLE, BzzGT, BzzGE and Swap.                   *
// * Example:                                                                 *
// *        float x[5] = {3.,2.,5.,1.,4.};                                    *
// *        Sort(5,x);                                                        *
// ****************************************************************************
template <class T>
inline unsigned char BzzLE(T a, T b)
{
	return (a <= b);
}

template <class T>
inline unsigned char BzzLT(T a, T b)
{
	return (a < b);
}

template <class T>
inline unsigned char BzzGE(T a, T b)
{
	return (a >= b);
}

template <class T>
inline unsigned char BzzGT(T a, T b)
{
	return (a > b);
}

template <class T>
void Sort(int n, T* x)
{
	if (n <= 1)return;
	int node, i, j, k, ik, jk;
	for (node = 1;node < n;node++)
	{
		i = node;
		j = ((i + 1) / 2) - 1;
		while (i != 0 && BzzLE(x[j], x[i]))
		{
			Swap(x + j, x + i);
			i = j;
			j = ((i + 1) / 2) - 1;
		}
	}
	for (i = n - 1;i >= 1;i--)
	{
		Swap(x + i, x);
		k = i - 1;
		ik = 0;
		jk = 1;
		if (k >= 2 && BzzGT(x[2], x[1]))
			jk = 2;
		while (jk <= k && BzzGT(x[jk], x[ik]))
		{
			Swap(x + jk, x + ik);
			ik = jk;
			jk = (2 * (ik + 1)) - 1;
			if (jk + 1 <= k)
				if (BzzGT(x[jk + 1], x[jk]))
					jk++;
		}
	}
}

// *******************************< Sort >*************************************
// * Purpose: To sort an array                                                *
// * Description: Chapter 5 Scientific C++.                                   *
// * It require the functions BzzLE, BzzGT, BzzGE and Swap.                   *
// * Example:                                                                 *
// *        float x[5] = {3.,2.,5.,1.,4.};                                    *
//	*			int iS[5] = {1,2,3,4,5};														*
// *        Sort(5,x,iS);                                                     *
// ****************************************************************************
template <class T>
void Sort(int n, T* x, int* iS)
{
	if (n <= 1)return;
	int node, i, j, k, ik, jk;
	for (node = 1;node < n;node++)
	{
		i = node;
		j = ((i + 1) / 2) - 1;
		while (i != 0 && BzzLE(x[j], x[i]))
		{
			Swap(x + j, x + i);
			Swap(iS + j, iS + i);
			i = j;
			j = ((i + 1) / 2) - 1;
		}
	}
	for (i = n - 1;i >= 1;i--)
	{
		Swap(x + i, x);
		Swap(iS + i, iS);
		k = i - 1;
		ik = 0;
		jk = 1;
		if (k >= 2 && BzzGT(x[2], x[1]))
			jk = 2;
		while (jk <= k && BzzGT(x[jk], x[ik]))
		{
			Swap(iS + jk, iS + ik);
			Swap(x + jk, x + ik);
			ik = jk;
			jk = (2 * (ik + 1)) - 1;
			if (jk + 1 <= k)
				if (BzzGT(x[jk + 1], x[jk]))
					jk++;
		}
	}
}

// ============================================================================
// =============================< PROTOTYPES >=================================
// ============================================================================
double BzzClock(void);
double BzzDiffClock(void);
double BzzGetUserTime(void);
double BzzGetKernelTime(void);
double BzzGetCpuTime(void);

float MachEpsFloat(void);
double MachEps(void);

void BzzMathVersionMessage(void);
int BzzMathVersionMajor(void);
int BzzMathVersionMinor(void);
int BzzMathVersionRevision(void);

void BzzMessageAndGetchar(const char* myFormat, ...);
void BzzPrint(const char* myFormat, ...);
//void BzzError(char *strFirst,char *strSecond);
//void BzzError(char *strFirst);
void BzzError(const char* myFormat, ...);
void BzzError(void);
void BzzWarning(const char* myFormat, ...);
void BzzWarning(void);

float Control(double value);
double Control(long double value);

void Sum(int n, float* lval, float* rval, float* result);
void Sum(int n, float* lvalAndResult, float* rval);
void Sum(int n, float* lvalRvalAndResult);
void Sum(int n, double* lval, double* rval, double* result);
void Sum(int n, double* lvalAndResult, double* rval);
void Sum(int n, double* lvalRvalAndResult);

void Difference(int n, float* lval, float* rval, float* result);
void Difference(int n, float* lvalAndResult, float* rval);
void Difference(int n, float* lval, float* rvalAndResult, int);
void Difference(int n, double* lval, double* rval, double* result);
void Difference(int n, double* lvalAndResult, double* rval);
void Difference(int n, double* lval, double* rvalAndResult, int);

float Dot(int n, float* lval, float* rval);
double Dot(int n, double* lval, double* rval);

void Product(int n, double lval, float* rval, float* result);
void Product(int n, double lval, float* rvalAndResult);
void Product(int n, double lval, double* rval, double* result);
void Product(int n, double lval, double* rvalAndResult);

void Division(int n, float* lval, double rval, float* result);
void Division(int n, float* lvalAndResult, double rval);
void Division(int n, double* lval, double rval, double* result);
void Division(int n, double* lvalAndResult, double rval);

float SqrtSumSqr(int n, float* x);
double SqrtSumSqr(int n, double* x);

inline void RemoveWarningWindow(void) { bzzWarningWindow = 0; }
inline void RemoveWarningWindowAndWarningPrint(void) { bzzWarningWindow = -1; }
inline void SetWarningWindow(void) { bzzWarningWindow = 1; }

inline void RemoveErrorWindow(void) { bzzErrorWindow = 0; }
inline void SetErrorWindow(void) { bzzErrorWindow = 1; }

double BzzPowInt(double x, int n);
inline double BzzPow2(double x)
{
	return x * x;
}
inline double BzzPow3(double x)
{
	return x * x * x;
}

inline double BzzPow4(double x)
{
	double x2 = x * x;
	return x2 * x2;
}

inline double BzzPow5(double x)
{
	double x2 = x * x;
	return x2 * x2 * x;
}

inline double BzzPow6(double x)
{
	double x2 = x * x;
	return x2 * x2 * x2;
}

inline double BzzPow7(double x)
{
	double x2 = x * x;
	return x2 * x2 * x2 * x;
}

inline double BzzPow8(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	return x4 * x4;
}

inline double BzzPow9(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	return x * x4 * x4;
}

inline double BzzPow10(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	double x5 = x4 * x;
	return x5 * x5;
}

double BzzCubicRoot(double x);
void Flinks(void);

char BulSto0(int np, float* x, float* y, float* prev, float* dprev);
char BulSto0(int np, double* x, double* y, double* prev, double* dprev);

void BzzPause(void);

double BzzLogGamma(double xx);
double BzzGamma(double xx);
double BzzLogFactorial(int n);
double BzzBinomial(int n, int k);
double BzzBinomialProbabilityDistribution(int n, int k, double p);
double BzzPoissonProbabilityDistribution(double lambda, double t, int k);
double BzzBeta(double z, double w);
double BzzIncompleteGammaP(double a, double x);
double BzzIncompleteGammaQ(double a, double x);
double BzzNormalGreaterThanFixed(double u);
double BzzChiSquareLowerThanFixed(int fd, double chiSquare);
double BzzChiSquareGreaterThanFixed(int fd, double chiSquare);
double BzzErf(double x);
double BzzErfc(double x);
double BzzIncompleteBeta(double a, double b, double x);
double BzzBiasedFactorForStandardDeviation(int nn);
double BzzOneWayTTest(int fd, double t);
double BzzOneWayFTest(int fdNum, int fdDen, double F);
double BzzTwoWayFTest(int fdNum, int fdDen, double F);
void BzzFibonacci(double L1, double delta, int N, BzzVector* Li);
double BzzFibonacci(double Li, double delta, int N);
int BzzFibonacci(double Li, double delta, double* Lf);
double BzzFibonacciGetL1GivenL3(double L3, double delta, int N);
/*
// ******************************< BzzTempFile >*******************************
// * Class for creating temporary files                                       *
// ****************************************************************************

// ============================================================================
// =============================< class BzzTempFile >==========================
// ============================================================================
class BzzTempFile
   {
private:
   static int countBzzTempFile;
   char *nameBzzTempFile;
public:
   // default constructor BzzTempFile ff[10];
   BzzTempFile(void){nameBzzTempFile = 0;}

   // constructor BzzTempFile f1("D:"),f2("c:\\temp\\");
   BzzTempFile(char *directoryTemp)
	  {NewFileName(directoryTemp);}

   // destructor
   ~BzzTempFile(void){delete nameBzzTempFile;}

   // provides the unique name f1.FileName
   char *FileName(void){return nameBzzTempFile;}

   // for using after default ff[0].NewFileName("D:");
   void NewFileName(char *directoryTemp);
   };
*/

// ***************************< BzzBaseClass >*********************************
// * Base Class for BzzPrint and BzzMessage                                   *
// ****************************************************************************

// ============================================================================
// ==========================< class BzzBaseClass >============================
// ============================================================================
class BzzBaseClass
{
public:
	// default constructor
	BzzBaseClass(void) {}

	// destructor
	~BzzBaseClass(void) {}
	void BzzPrint(void);
	void BzzPrint(const char* myFormat, ...);
	void BzzMessage(void)
	{
		if (bzzMessageYesNo == 'Y')BzzPrint();
	}
	void BzzMessage(const char* myFormat, ...);
	virtual void ObjectBzzPrint(void) = 0;
};

// *************************< BzzCloseFiles >**********************************
// * Class for closing files                                                  *
// ****************************************************************************

// ============================================================================
// =======================< class BzzCloseFiles >==============================
// ============================================================================

class BzzCloseFiles
{
public:
	// default constructor
	BzzCloseFiles(void)
	{
		if ((bzzFileOut = fopen("BzzFile.txt", "w")) == NULL)
			BzzError("%s", BZZ_ERR_OPEN_FILE);
	}

	void operator()(char* file);

	// destructor
 //   ~BzzCloseFiles(void){_fcloseall();}
};
extern BzzCloseFiles bzzFilePrint;

#endif // BZZ_UTILITY_HPP