// BZZMATH: Release 7.0

//	=========================< OBzzOdeRungeKutta.hpp >==========================
//	* Header FILE for Runge Kutta	methods													*
//	* Description:																		 			*
//	*					Metodi Numerici e Software in C++ (Capitoli 29, 30)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathAdvanced\Ode\OdeRungeKutta\		*
// *				OdeRungeKutta.hpp														*
// * Tests: BzzMath\Examples\BzzMathAdvanced\Ode\OdeNonStiffTests\		*
// *				OdeNonStiffTests.hpp													*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-1994	Date Written

//	============================================================================
//	* These classes allow to solve initial value systems of first order ODEs:	*
//	*					dy/dt = f(y,t)																*
//	****************************************************************************
//	****** Constructors for BzzOdeRK:												*
//	* BzzOdeRK o(y0,t0,Sys);															*
//	****** Constructors for BzzOdeRKF:												*
//	* BzzOdeRKF o(y0,t0,Sys);															*
//	****** Constructors for BzzOdeRKM:												*
//	* BzzOdeRKM o(y0,t0,Sys);															*
//	****** Constructors for BzzOdeRKS:												*
//	* BzzOdeRKS o(y0,t0,Sys);															*
//	****************************************************************************
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
//	****************************************************************************
//	* Before defining an object of these classes it is necessary to declare		*
// * the function in which the components of the f vector are calculated.		*
//	* The prototype of the function is:														*
//	*		void <name>(BzzVector &y,double t,BzzVector &f);			*
//	* This function has to be defined elsewhere and the vector						*
// * f[1], f[2], ..., f[N] has to be calculated in function of the				*
//	* scalar value t and vector elements y[1], y[2], ..., y[N].						*
//	* Example:																						*
//	* void MySystem(BzzVector &y,double t,BzzVector &f)				*
//	*	{																								*
//	*	f[1] = y[1]*t + y[2];																	*
//	*	f[2] = -2.*y[1]*y[2];																	*
//	*	}																								*
// * NOTE:																							*
// * - The user must not modify the vector y.											*
// * - The coefficients of the vectors y and f are accessible using [] or ().	*
// * - [] is more efficient but the user is responsible of the correct range.	*
//	****************************************************************************
//	****************************************************************************
//	***** Use:																						*
//	* Before using an object of these classes it is necessary to define			*
// * the independent and dependent variables at starting point and				*
//	* the vector y, solution of the problem.												*
//	* Example:																						*
//	* double t0 = 0.;																				*
//	* BzzVector y0(3,0.,1.,2.);														*
//	* BzzVector y;																			*
// * NOTE:																							*
// * - It is not necessary to give the correct dimension for the vector y.		*
//	****************************************************************************
//	***** 1. Normal computation.																*
//	*	y = o(tOut);																				*
//	* Output values of y(t) at t = tOut														*
//	****************************************************************************
//	***** 2. Fixed integration step															*
//	* y = o(tOut,numSteps);																		*
//	****************************************************************************
//	****************************************************************************
//	* Here following is a list of the optionals											*
//	****************************************************************************
// * The tollerances tollAbs and tollRel are used by the code in a local		*
// * error test at each step which require:												*
// * 		abs(local error in y(i)) < tollRel(i)*abs(y(i)) + tollAbs(i)			*
// * The default values of tollRel and tollAbs are:									*
// * tollRel = 100.*MachEps(); tollAbs = 1.e-10;										*
// * The user can modify these default values inserting by setting the			*
// * appropriate values.																		*
//	* Example:																						*
//	* o.SetTolAbs(tollAbs);																	*
//	* o.SetTolRel(tollRel);																	*
// * tollAbs and tollRel can be either scalar (double) or							*
// * vectors (BzzVector).																*
// * The user could modify the values of tollAbs and tollRel also during the	*
// * integration of the system by using the functions:								*
//	*		o.SetTolAbs(tollAbs);																*
//	*		o.SetTolRel(tollRel);																*
//	* Also in this case tollAbs and tollRel can be either scalar (double)		*
// * or vectors (BzzVector)															*
//	****************************************************************************
//	*		o.SetHMin(hMin); (double hMin)													*
//	*		o.SetHMax(hMax); (double hMax)													*
//	****************************************************************************
//	*		o.SetMaxStep(maxStep);																*
//	****************************************************************************
//	****************************************************************************
//	***** Get Information.																		*
// * You can use the functions BzzPrint or BzzMessage for a summary of the		*
// * situation after the last call of the solver.										*
//	* Example:																						*
//	* y = o(tOut);																					*
//	* o.BzzPrint("Results");																	*
// * You can also get separate information by using the following functions:	*
//	*																									*
// * int numStep = GetNumStep();																*
// * int numFunction = GetNumFunction();													*
// * int numErrorTestFailures = CountOk();												*
// * int numErrorTestFailures = CountOk();												*
//	****************************************************************************

#ifndef BZZODERK_HPP
#define BZZODERK_HPP

class BzzOdeRungeKutta : public BzzBaseClass
{
protected:
	static const char* const BZZ_ERROR;
	static const double 	DEFAULT_TOL_ABS,
		DEFAULT_TOL_REL;

	int numVariables;
	BzzVector	y0,
		y,
		yOut,
		dydt,
		errAbs,
		errRel;
	double	t0,
		t,
		h, H;
	int	numSys,
		numOK,
		numNoOK;

	int maxStep;

	void (*sysDiff)(BzzVector& y, double t, BzzVector& f);

	virtual void StepControl(double tF) = 0;
	virtual void Step(double tF, int nStep) = 0;

	void ControlData(void);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default
	BzzOdeRungeKutta(void)
	{
		BzzError("Default constructor not implemented");
	}

	BzzOdeRungeKutta(BzzVector& y0, double t0,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================

	//	============================================================================
	//	===============================< functions >================================
	//	============================================================================
	int GetNumFunction(void) { return numSys; }
	int CountOK(void) { return numOK; }
	int CountNoOK(void) { return numNoOK; }
	virtual void ObjectBzzPrint(void);
	void SetMaxStep(int maxS) { maxStep = maxS; }
	void SetTolAbs(double tA);
	void SetTolAbs(BzzVector* tA);
	void SetTolAbs(const BzzVector& tA);
	void SetTolRel(double tR);
	void SetTolRel(BzzVector* tR);
	void SetTolRel(const BzzVector& tR);

	BzzVector operator ()(double tF);
	BzzVector operator ()(double tF, int nStep);
};

//	============================================================================
//	* BzzOdeRK class (Runge - Kutta)															*
//	============================================================================
class BzzOdeRK : public BzzOdeRungeKutta
{
private:
	BzzVector	va, vb,
		dySave, yTemp;

	void RungeKutta(const BzzVector& yin, double tin,
		const BzzVector& dydtin, BzzVector* yout);

	void RungeKuttaStepControl(void);
	void StepControl(double tF);
	void Step(double tF, int nStep);

public:

	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
	BzzOdeRK(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeRungeKutta(y00, t00, ptrSys),
		va(y00.Size()),
		vb(y00.Size()),
		dySave(y00.Size()),
		yTemp(y00.Size())
	{}
};

//	============================================================================
//	* BzzOdeRKF class (Runge - Kutta - Fehlberg)											*
//	============================================================================
class BzzOdeRKF : public BzzOdeRungeKutta
{
private:
	BzzVector	y4,
		k1, k2, k3, k4, k5, k6,
		va, vb;

	void RungeKuttaFehlberg(const BzzVector& yin, double tin,
		const BzzVector& dydtin, BzzVector* yout);

	void RungeKuttaFehlbergStepControl(void);
	void StepControl(double tF);
	void Step(double tF, int nStep);

public:

	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
	BzzOdeRKF(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeRungeKutta(y00, t00, ptrSys),
		y4(y00.Size()),
		k1(y00.Size()),
		k2(y00.Size()),
		k3(y00.Size()),
		k4(y00.Size()),
		k5(y00.Size()),
		k6(y00.Size()),
		va(y00.Size()),
		vb(y00.Size())
	{}
};
//	============================================================================
//	* BzzOdeRKM class (Runge - Kutta - Merson)									*
//	============================================================================
class BzzOdeRKM : public BzzOdeRungeKutta
{
private:
	BzzVector	y4,
		k1, k2, k3, k4, k5, k6,
		va, vb;

	void RungeKuttaMerson(const BzzVector& yin, double tin,
		const BzzVector& dydtin, BzzVector* yout);

	void RungeKuttaMersonStepControl(void);
	void StepControl(double tF);
	void Step(double tF, int nStep);

public:

	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
	BzzOdeRKM(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeRungeKutta(y00, t00, ptrSys),
		y4(y00.Size()),
		k1(y00.Size()),
		k2(y00.Size()),
		k3(y00.Size()),
		k4(y00.Size()),
		k5(y00.Size()),
		k6(y00.Size()),
		va(y00.Size()),
		vb(y00.Size())
	{}
};

//	============================================================================
//	* BzzOdeRKS class (Runge - Kutta - Sarafyan)									*
//	============================================================================
class BzzOdeRKS : public BzzOdeRungeKutta
{
private:
	BzzVector	y4,
		k1, k2, k3, k4, k5, k6,
		va, vb;

	void RungeKuttaSarafyan(const BzzVector& yin, double tin,
		const BzzVector& dydtin, BzzVector* yout);

	void RungeKuttaSarafyanStepControl(void);
	void StepControl(double tF);
	void Step(double tF, int nStep);

public:

	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
	BzzOdeRKS(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeRungeKutta(y00, t00, ptrSys),
		y4(y00.Size()),
		k1(y00.Size()),
		k2(y00.Size()),
		k3(y00.Size()),
		k4(y00.Size()),
		k5(y00.Size()),
		k6(y00.Size()),
		va(y00.Size()),
		vb(y00.Size())
	{}
};

#endif // BZZODERK_HPP