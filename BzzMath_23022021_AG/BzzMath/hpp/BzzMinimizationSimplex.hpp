// BZZMATH: Release 7.0

//	================< BzzMinimizationSimplex.hpp >========================
//	* The class BzzMinimizationSimplex is used for unconstrained			*
//	* optimization																					*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 23, 24)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: BzzMath\examples\BzzMathAdvanced\Minimization\Unconstrained\	*
// *			MinimizationSimplex\MinimizationSimplex.cpp				*
// * Tests: BzzMath\examples\BzzMathAdvanced\Minimization\Unconstrained\		*
// *			MinimizationSimplexTest\MinimizationSimplexTest.cpp	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-1994	Date Written.
// Add Save functions.
//	11-1998	Bug fixed in GetBzzMinimumF function.
//	11-1998	Add GetMethodUsed.

/////// Release 4.0
//	02-2000	Added SetTasksPrint, SetSubTasksPrint.
//	02-2000	Simplex in Quasi Newtom more robust search.

/////// Release 5.0

//	============================================================================
//	******	constructors:																		*
//	* BzzMinimizationSimplex m; // default											*
//	* m(x0,FunMin);																				*
//	* BzzMinimizationSimplex m(x0,FunMin);											*
//	****************************************************************************
//	****** Functions for minimization:														*
//	* char control = m();																		*
//	* char control = m(niter);																	*
//	****************************************************************************
//	* Other access functions:																	*
//	* int numIterations = m.IterationCounter();											*
//	* double F = m.GetBzzMinimumF();															*
//	* double F = m.Solution(&x);																*
//	* m.BzzPrint("Comments");																	*
//	* m.BzzMessage("Comments");																*
//	****************************************************************************
//	* Functions for assigning tollerances:													*
//	* SetXRelativeTolerance(gT);																		*
//	* SetXAbsoluteTolerance(nT);																			*
//	* SetFunctionMinimumTolerance(fT);																*
//	****************************************************************************
//	* Functions for save partial results:													*
//	* Save("SAVE.TXT");																			*
//	****************************************************************************
//	* Restart(BzzVector x0);															*
//	****************************************************************************

#ifndef MINIMUM_MULTI_DOUBLE_SIMPLEX_HPP
#define MINIMUM_MULTI_DOUBLE_SIMPLEX_HPP

//	============================================================================
//	======================< class BzzMinimizationSimplex >=======================
//	============================================================================
class BzzMinimizationSimplex : public BzzBaseClass
{
private:
	enum MethodType
	{
		SIMPLEX_METHOD,
		QUASI_NEWTON_METHOD
	}methodType;

	enum MethodUsed
	{
		QUASI_NEWTON_USED = 1,
		SIMPLEX_USED = 2,
		QUASI_NEWTON_SIMPLEX_USED = 3
	}methodUsed;

	enum GradientHessianType
	{
		GRAD_AND_HESS_NUMERICAL,
		GRAD_ANALYTICAL_HESS_NUMERICAL,
		GRAD_AND_HESS_ANALYTICAL
	}gradHessType;

	enum GradientStatus
	{
		GRADIENT_AVAILABLE,
		GRADIENT_NOT_KNOWN
	}gradientStatus;

	enum MethodStatusBzzMinimumQuasiNewton
	{
		QUASI_NEWTON_INITIALIZE,
		QUASI_NEWTON_UPDATE,
		QUASI_NEWTON,
		NEWTON,
		MONODIMENSIONAL_FIRST,
		MONODIMENSIONAL_SECOND,
		NEGATIVE,
		GRADIENT,
		GRADIENT_MONO_FIRST,
		GRADIENT_MONO_SECOND,
		DOG_LEG
	}methodStatusBzzMinimumQuasiNewton;

	enum SimplexStatus
	{
		REFLECTION,
		EXPANSION,
		CONTRACTION,
		REVERSE_CONTRACTION,
		SIMPLEX_CONTRACTION,
		REINITIALIZATION
	}simplexStatus;

	static const char* const BZZ_ERROR;

	FILE* bzzFileSave;
	char saveResults;

	int	numVariables,
		niter, iter, iterTotal,
		iterGrad,
		iterHess,
		iterSolve,
		iterPositive,
		iterQuasi,
		iQuasi,
		nQuasi,
		reinitialize,
		iterationStatus,
		iInitialize,
		jInitialize,
		iContraction,
		jContraction,
		printTasks,
		printSubTasks;

	char	control,
		controlHess,
		stop;

	BzzVector	xi,
		xStart,
		gi,
		pi,
		xi1,
		dg,
		xDimensions,
		h,
		aux1, aux2,
		dgi,	// gi1 - gi
		dxi,	// xi1 - xi
		xiBase,
		giBase,
		qqq;

	BzzVector	vS,	 // the best
		vL,	 // largest
		vSL,	// second largest
		* v,	 // others
		vB,	 // base
		vR,	 // reflection
		vE,	 // expansion
		vC,	 // contraction
		x,
		xd;

	double	fi,
		fStart,
		fi1,
		fDimensions,
		giTpi,
		giTgi,
		giTdg,
		giNorm2,
		pMax,
		piNorm2,
		gradToll,
		newToll,
		fToll,
		u, w, z, fz, fu,
		alfa, beta,
		giTdxi, dgiTdxi;

	double	fS,	 // the best
		fL,	 // largest
		fSL,	// second largest
		* f,	 // others
		fB,	 // base
		fR,	 // reflection
		fE,	 // expansion
		fC,	 // contraction
		fI,	 // reinitialization
		fx,
		xTollAbs,
		xTollRel,
		step1, step2;

	BzzFactorizedGillMurray Gi;
	BzzMatrixSymmetric S;

	double (*ptrFun)(BzzVector& x);
	void (*ptrGradient)(BzzVector& xx, BzzVector& gg);
	void (*ptrHessian)(BzzVector& xx, BzzMatrixSymmetric& GG);

	void InitializeQuasiNewton(void);
	char ControlBzzMinimum(void);
	char ControlPrevision(void);
	char ControlGradient(void);
	char ControlMono(void);
	void SetBzzMinimum(double fm, BzzVector& xm);
	void StepH(double eta); // calcola h
	void FindStep(void);
	void GradientForward(void);
	void GradientCentral(void);
	void Hessian(void);
	void GradientAndHessianFirst(void);
	void GradientAndHessianSecond(void);
	char ControlQuasiNewton(void);
	double FindAlfa(void);
	void GetGradient(void);
	void GetHessian(void);
	char QuasiNewtonMethod(int ni);
	char SimplexMethod(int ni);

	char InitializeSimplex(void);
	void Baricenter(void);
	void Insert(BzzVector& x, double fx);
	void FindSL(void);
	char SimplexContraction(void);
	char ControlBzzSimplex(void);
	char MinimumSolve(int ni);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default constructor
	BzzMinimizationSimplex(void);

	// copy constructor
	BzzMinimizationSimplex
	(const BzzMinimizationSimplex& rval);

	// numerical gradient and Hessian; Simplex as base method
	BzzMinimizationSimplex
	(BzzVector x00, double (*ptr)(BzzVector& x));

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzMinimizationSimplex(void);

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================
	double GetSolution(BzzVector* xx);
	double GetBzzMinimumF(void);
	double GetMinimumF(void) { return GetBzzMinimumF(); }
	char GetMethodUsed(void) { return methodUsed; }
	int IterationCounter(void) { return iterTotal; }
	int GetIterTotal(void)
	{
		return iterTotal;
	}
	virtual void ObjectBzzPrint(void);
	void Save(char* sv);
	void StopTaskPrint(void) { printTasks = 0; }
	void StopSubTasksPrint(void) { printSubTasks = 0; }
	void SetTasksPrint(void) { printTasks = 1; }
	void SetSubTasksPrint(int psb = 1)
	{
		if (psb > 0)
			printSubTasks = psb;
		else
			printSubTasks = 1;
	}

	//	============================================================================
	//	**********************< Modifying functions >*******************************
	//	============================================================================
	void Restart(BzzVector x0);
	void operator()(BzzVector x0, double (*ptr)(BzzVector& x));
	void SetXRelativeTolerance(double gT);
	void SetXAbsoluteTolerance(double nT);
	void SetFunctionMinimumTolerance(double fT);
	char operator ()(void)
	{
		return (*this)(0);
	}
	char operator ()(int ni);

	// return < 0 Search failed before reaching the solution
	// return = 0 Found problems
	// return = 1 Maximum number of iterations
	// return = 2 grad <= GRAD_TOL
	// return = 3 fSolution <= fToll
	// return = 4 dNew < NEW_TOL
};

#endif // MINIMUM_MULTI_DOUBLE_SIMPLEX_HPP