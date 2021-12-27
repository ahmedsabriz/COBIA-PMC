// BZZMATH: Release 7.0

//	========================< BzzIntegral.hpp >===========================
//	* BzzIntegral and BzzIntegralGauss classes							*
//	* for numerical quadrature																	*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 27, 28)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	*																									*
//	* Examples: BzzMath\Examples\BzzMathAdvanced\DefiniteIntegral\					*
//	*				Integral\Integral.cpp										*
//	* Tests: BzzMath\Examples\BzzMathAdvanced\DefiniteIntegral\						*
//	*			IntegralTests\IntegralTests.cpp								*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1993	Date Written.
//	10-2002	Added class BzzIntegralGaussBF.

#ifndef INTEGRAD_HPP
#define INTEGRAD_HPP
//	============================================================================
//	***************************< BzzIntegral class >**********************
//	============================================================================
//	****** Constructors:																			*
//	* BzzIntegral o1(BzzIntegralFun);	// or o1(BzzIntegralFun,maxFun) *
//	* BzzIntegral o1(BzzIntegralFun,errAbs,errRel,maxFun,hMin);		*
//	****************************************************************************
//	***** Access function:																		*
//	* double I = o1(tA,tB);																		*
//	* int functionCallsCounter = o1.FunctionCallsCounter();							*
//	* double hMinUsed = o1.HMinUsed();														*
//	****************************************************************************
//	* Function prototype:																		*
//	* double BzzIntegralFun(double t);												*
//	****************************************************************************
class BzzIntegral
{
private:
	static const char* const BZZ_ERROR;
	double tA, tB, H,
		hMin,
		errAbs,
		errRel,
		integraltAtB,
		hMinUsed;
	double integralAbs;
	double T[6], M[4];
	double hh[6], prevision[6], dprevision[6];
	int functionCallsCounter,
		maxFunction;
	double (*funzInt)(double t);

	void ControlData(int mxFz);
	double Int12(double tI, double tF, double* f);
	double IntLeftRight(double tI, double tF);
public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
	BzzIntegral(double (*ptrInt)(double t),
		int mxFz = 100000);

	BzzIntegral(double (*ptrInt)(double t), double errA, double errR,
		int mxFz = 100000, double hMi = 100. * BZZ_TINY_FLOAT);

	//	============================================================================
	//	****************************< Functions >***********************************
	//	============================================================================
	void SetMaxFuncions(int numMax);
	void SetTolAbs(double tollA);
	void SetTolRel(double tollR);
	void SetMinH(double hm);
	int FunctionCallsCounter(void) { return functionCallsCounter; }
	double HMinUsed(void) { return hMinUsed; }
	double operator ()(double tI, double tF);
};

//	============================================================================
//	***************************< BzzIntegralGauss class >*****************
//	============================================================================

//	****** Constructors:																			*
//	* BzzIntegral o1(BzzIntegralFun);	// or o1(BzzIntegralFun,maxFun) *
//	* BzzIntegral o1(BzzIntegralFun,errAbs,errRel,maxFun,hMin);		*
//	****************************************************************************
//	***** Access function:																		*
//	* double I = o1(tA,tB);																		*
//	* int functionCallsCounter = o1.FunctionCallsCounter();							*
//	* double hMinUsed = o1.HMinUsed();														*
//	****************************************************************************
//	* Function prototype:																		*
//	* double BzzIntegralFun(double t);												*
//	****************************************************************************
class BzzIntegralGauss
{
private:
	static const char* const BZZ_ERROR;
	double tA, tB, H,
		hMin,
		errAbs,
		errRel,
		integraltAtB,
		hMinUsed;
	double integralAbs;
	int functionCallsCounter,
		maxFunction;
	double (*funzInt)(double t);

	void ControlData(int mxFz);
	double GaussKronrod(double tI, double tF);
public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
	BzzIntegralGauss(double (*ptrInt)(double t),
		int mxFz = 100000);

	BzzIntegralGauss(double (*ptrInt)(double t), double errA, double errR,
		int mxFz = 100000, double hMi = 100. * BZZ_TINY_FLOAT);

	//	============================================================================
	//	****************************< Functions >***********************************
	//	============================================================================
	void SetMaxFuncions(int numMax);
	void SetTolAbs(double tollA);
	void SetTolRel(double tollR);
	void SetMinH(double hm);
	int FunctionCallsCounter(void) { return functionCallsCounter; }
	double HMinUsed(void) { return hMinUsed; }
	double operator ()(double tI, double tF);
};

//	============================================================================
//	**************************< BzzIntegralGaussBF class >****************
//	============================================================================

//	****** Constructors:																			*
//	* BzzIntegral o1(BzzIntegralFun);	// or o1(BzzIntegralFun,maxFun) *
//	* BzzIntegral o1(BzzIntegralFun,errAbs,errRel,maxFun,hMin);		*
//	****************************************************************************
//	***** Access function:																		*
//	* double I = o1(tA,tB);																		*
//	* int functionCallsCounter = o1.FunctionCallsCounter();							*
//	* double hMinUsed = o1.HMinUsed();														*
//	****************************************************************************
//	* Function prototype:																		*
//	* double BzzIntegralFun(double t);												*
//	****************************************************************************
class BzzIntegralGaussBF
{
private:
	static const char* const BZZ_ERROR;
	double	tA, tB, H,
		hMin,
		errAbs,
		errRel,
		integraltAtB,
		hMinUsed;
	double integralAbs;
	int functionCallsCounter,
		maxFunction;

	//	BzzVector z,tt;
	//	BzzMatrix Q,D;
	double (*funzInt)(double t);
	BzzVector int7;
	int numThreads;
	void ControlData(int mxFz);
	double GaussBase(double tI, double tF);
	double GaussBF(double i7, double tI, double tF);
	double GaussFirst(double tI, double tF);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
	BzzIntegralGaussBF(double (*ptrInt)(double t),
		int mxFz = 100000);

	BzzIntegralGaussBF(double (*ptrInt)(double t), double errA, double errR,
		int mxFz = 100000, double hMi = 100. * BZZ_TINY_FLOAT);

	//	============================================================================
	//	****************************< Functions >***********************************
	//	============================================================================
	void SetMaxFuncions(int numMax);
	void SetTolAbs(double tollA);
	void SetTolRel(double tollR);
	void SetMinH(double hm);
	int FunctionCallsCounter(void) { return functionCallsCounter; }
	double HMinUsed(void) { return hMinUsed; }
	double operator ()(double tI, double tF);
};

#endif // INTEGRAD_HPP