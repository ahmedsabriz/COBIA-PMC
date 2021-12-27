// BZZMATH: Release 7.0

//	========================< BzzIntegral.hpp >===========================
//	* for numerical quadrature																	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1993	Date Written.
//	10-2002	Added class BzzIntegralGaussBF.
//	07-2013	Added class BzzIntegralGaussLobatto.
//	07-2013	Mofified classes BzzIntegralGauss,BzzIntegral and BzzIntegralGaussBF

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
class BzzIntegral : public BzzBaseClass
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
		maxFunction,
		numThreads;
	double (*funzInt)(double t);
	BzzVector ttt, fff;
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

	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
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
class BzzIntegralGauss : public BzzBaseClass
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
		maxFunction,
		numThreads;
	double (*funzInt)(double t);
	BzzVector ttt, fff;
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

	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
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
class BzzIntegralGaussBF : public BzzBaseClass
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
		maxFunction,
		numThreads;

	//	BzzVector z,tt;
	//	BzzMatrix Q,D;
	double (*funzInt)(double t);
	BzzVector int7;
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

	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
};

//	============================================================================
//	========================< class BzzGaussLobatto >===========================
//	============================================================================

class BzzIntegralGaussLobatto : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const double BZZ_LOG10;
	static const double BZZ_ABS_TOL;
	static const double BZZ_REL_TOL;
	static int count; // for whoAmI
	int whoAmI;
	double	tA,
		tB,
		H,
		integralBase,
		tLeft, fLeft, integralLeft, // tA tLeft
		tRight, fRight, integralRight, // tLeft tRight
		integralCentral, // tRight tB
		integral,
		absTol, relTol;
	int functionCallsCounter,
		maxFunction,
		numThreads;
	double (*funzInt)(double t);
	BzzVector int6L, int6C, int6R;

	double IntegralLeft(double t1, double tm, double fm, double t2, double f2);
	double IntegralCentral(double t1, double f1, double t2, double f2);
	//	double IntegralCentral(double I,double t1,double f1,
	//		double tm,double fm,double t2,double f2);
	double IntegralCentralLeft(double t1, double f1, double t2, double f2);
	double IntegralCentralRight(double t1, double f1, double t2, double f2);
	//	double IntegralCentral5(double t1,double f1,double t2,double f2);
	//	double IntegralCentral5(double I,double t1,double f1,
	//		double tm,double fm,double t2,double f2);
	//	double IntegralRight(double t1,double f1,double t2);
	double IntegralRight(double t1, double f1, double tm, double fm, double t2);

public:
	void SetAbsTol(double aT)
	{
		absTol = aT;
	}
	void SetRelTol(double aR)
	{
		relTol = aR;
	}
	int FunctionCallsCounter(void) { return functionCallsCounter; }
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
	BzzIntegralGaussLobatto(double (*ptrInt)(double t),
		int mxFz = 100000)
	{
		absTol = BzzIntegralGaussLobatto::BZZ_ABS_TOL;
		relTol = BzzIntegralGaussLobatto::BZZ_REL_TOL;
		funzInt = ptrInt;
		maxFunction = mxFz;
		integralBase = 0.;
		integral = 0.;
		functionCallsCounter = 0;
	}
	double operator ()(double tI, double tF);

	//	============================================================================
	//	====================< Non-modifying functions >=============================
	//	============================================================================
	double GetIntegral(void)
	{
		return integral;
	}
	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void) {}
};

#endif // INTEGRAD_HPP