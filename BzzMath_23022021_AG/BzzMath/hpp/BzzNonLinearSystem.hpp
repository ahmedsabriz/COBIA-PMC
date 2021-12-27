// BZZMATH: Release 7.0

//	=======================< BzzNonLinearSystem.hpp >=====================
//	* BzzNonLinearSystem class															*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 25, 26)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\BzzAdavanced\										*
// *					NonLinearSystems\NonLinearSystem.cpp						*
// * Tests: c:\bzzmath\examples\BzzAdavanced\											*
// *					NonLinearSystems\NonLinearSystemTests.cpp					*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1995	Date Written
// 07-1996	Introduction of upper and lower variables' constraints.
//	05-1998	Bug fixed in numFunctionsForGradient initialization.
//	06-1998	Introduction of new functions Restart.
//	10-1998	Bug fixed in Controlxi function with CONSTRAINTS.
//	10-1998	Bug fixed in Controlxi1 function with MINIMUM_CONSTRAINTS.
//	05-1999	Control in dimensions in the function SetConstraints.
//	01-2000	Added SetTasksPrint, SetSubTasksPrint.
//	03-2001	SetWeights.

////////////////// Release 5.0
// 09-2003	Added operator () for restarting.
// 11-2004	Added SetMaxNewtonCallsNumber function.

////////////////// Release 6.0
// 02-2009	Added GetWeights function.
// 07-2009	Added OpenMP parallelization.
// 11-2011	Added GetCalculationState function.
// 11-2011	Added GetMaximumResidual function.
// 11-2012	Added BzzNonLinearUnderdimensionedSystem class.

#ifndef NON_LINEAR_SYSTEM_DOUBLE_HPP
#define NON_LINEAR_SYSTEM_DOUBLE_HPP

//	============================================================================
//	=================< class BzzNonLinearSystemJacobian >=======================
//	============================================================================

class BzzNonLinearSystemJacobian
{
private:
	double zerDer, eta2, hj, xh, xdh, hInv;
	BzzVector aaa, bbb;
public:
	void (*ptrSystem)(BzzVector& xx, BzzVector& ff);
	void GetJacobian(int jStart, int jEnd, BzzVector& x, BzzVector& xDimensions,
		BzzVector& fi, BzzMatrix& Ji);
};

//	============================================================================
//	=======================< class BzzNonLinearSystem >===================
//	============================================================================

class BzzNonLinearSystem : public BzzBaseClass
{
	friend double MonSystem(double t);
protected:
	enum BzzNonLinearSystemCalculationState
	{
		INITIALIZATION_STATE = 0,
		CONTINUATION_STATE = 1,
		NEWTON_OK_STATE = 2,
		QUASI_NEWTON_OK_STATE = 3,
		GRADIENT_OK_STATE = 4,
		PHINEW_OK_STATE = 5,
		PHIW_OK_STATE = 6,
		DUBIOUS_STATE = 7,
		MAX_NEWTON_CALLS = 8,

		EXCESSIVE_WORK_STATE = -1,
		STOP_FOUND = -2,
		NOT_INITIALIZED = -3,
		STOP_RESTART = -4,
		STOP_FOR_BAD_CONVERGENCE = -5
	}nlsCalculationState;

	enum BzzNonLinearSystemJacobianType
	{
		JACOBIAN_NUMERICAL,
		JACOBIAN_ANALYTICAL
	}nlsJacobianType;

	enum BzzNonLinearSystemFactorizationMethod
	{
		FACTORIZATION_GAUSS,
		FACTORIZATION_LQ
	}nlsFactorizationMethod;

	enum BzzNonLinearSystemTollerance
	{
		SCALAR_A_SCALAR_R,
		ARRAY_A_SCALAR_R,
		SCALAR_A_ARRAY_R,
		ARRAY_A_ARRAY_R
	}nlsTollerance;

	enum BzzNonLinearSystemConstraintsState
	{
		NO_CONSTRAINTS,
		MINIMUM_CONSTRAINTS,
		MAXIMUM_CONSTRAINTS,
		MIN_MAX_CONSTRAINTS,
	}nlsConstraintsState;

	enum BzzNonLinearSystemMethodState
	{
		NEWTON,
		QUASI_NEWTON,
		MONODIMENSIONAL_FIRST,
		MONODIMENSIONAL_SECOND,
		GRADIENT,
		GRADIENT_MONODIMENSIONAL_FIRST,
		GRADIENT_MONODIMENSIONAL_SECOND,
		DOG_LEG
	}nlsMethodState;

	enum BzzNonLinearSystemWeightsType
	{
		AUTOMATIC_WEIGHTS,
		FIXED_WEIGHTS
	}nlsWeightsType;

	static const char* const BZZ_ERROR;

	static const double	DEFAULT_ABSOLUTE_TOLERANCE,
		DEFAULT_RELATIVE_TOLERANCE,
		NEWTON_TOLERANCE,
		DEFAULT_MAX_CONDITION,
		DEFAULT_PHI_RELATIVE_TOLERANCE,
		DEFAULT_PHI_ABSOLUTE_TOLERANCE,
		DEFAULT_MAX_W,
		BETA,
		GAMMA,
		MU;

	static const int	MAX_GRADIENT,
		MAX_BROYDEN_REGULAR,
		DEFAULT_MAX_FUNCTIONS_IN_THE_CALL,
		DEFAULT_MAX_FUNCTIONS_TOTAL;

	char stopForBadConvergence;

	int	numVariables,
		maxNumFunctionsInTheCall,
		maxNumFunctionsTotal,
		numFunctionsInTheCall,
		numFunctionsTotal,
		numQuasiNewtons,
		numNewtons,
		numGradients,
		numGradientsTrials,
		numFunctionsForGradient,
		numMonodimensional,
		numAnalyticalJacobians,
		numNumericalJacobians,
		numGaussFactorizations,
		numLQFactorizations,
		numLinearSolutions,
		numTrials,
		printTasks,
		printSubTasks,
		maxNewtonCallsNumber,
		numThread;

	BzzNonLinearSystemJacobian* getJac;

	char	//control,
		jacobianUpatating,
		singular,
		stepReduced,
		settingTollerance,
		stop,
		stopRestart,
		noJacobian;

	double	phiW,
		phiNew,
		phiLQ,
		phi1W,
		phi1New,
		phi1LQ,
		phiNewAbsoluteTollerance,
		phiNewRelativeTollerance,
		phiWAbsoluteTollerance,
		phiWAbsoluteTolleranceInStartingPoint,
		phiWRelativeTollerance,
		gradientToll,
		newtonToll,
		tollA,
		tollR,
		phiAbsToll,
		phiRelToll,
		xiNorm2,
		dxiNorm2,
		giNorm2,
		condition,
		sqrtInvSize,
		u, v, z, fz, fu;

	BzzVector	x0v, x0, f0v, f0,
		xi, fi,
		xi1, fi1,
		xMin, xMax, xDimensions,
		pi,
		pi1,
		dxi, dfi,
		gi,
		dLQ,
		weights,
		errorWeightBzzVector,
		tollAbs,
		tollRel,
		aux, aux1;

	BzzFactorizedGauss JiGauss;
	BzzFactorizedLQ JiLQ;
	BzzMatrix Ji;

	void (*ptrSystem)(BzzVector& xx, BzzVector& ff);
	void (*ptrJacobian)(BzzVector& xx, BzzMatrix& JJ);

	void Initialize(void);
	void ReInitialize(void);
	void BzzPrintErrorState(void);
	void ControlObjectiveFunctions(void);
	char ControlPrevision(double gamma);
	double ControlNewton(void);
	void Controlxi(void);
	void Controlxi1(void);
	char Controlxi1Newton(void);
	void PhiWInxi1(void);
	void PhiNewInxi1(void);
	void Weights(void);
	void ErrorWeightBzzVector(void);
	void GetJacobian(void);
	void NumericalJacobian(void);
	void NewtonPrevision(void);
	void QuasiNewtonPrevision(void);
	char NonLinearSolve(int nS);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================

		// default constructor
	BzzNonLinearSystem(void);

	// copy constructor
	BzzNonLinearSystem
	(const BzzNonLinearSystem& rval);

	// Numerical Jacobian
	BzzNonLinearSystem
	(const BzzVector& x0,
		void (*ptrS)(BzzVector& xx, BzzVector& ff));

	// Analytical Jacobian
	BzzNonLinearSystem
	(const BzzVector& x0,
		void (*ptrS)(BzzVector& xx, BzzVector& ff),
		void (*ptrJ)(BzzVector& xx, BzzMatrix& JJ));

	//	============================================================================
	//	*********************************< destructor >*****************************
	//	============================================================================
	~BzzNonLinearSystem(void)
	{
		if (numThread != 0)
#if BZZ_COMPILER != 0
			delete[] getJac;
#else
			//delete[numThread] getJac;
			delete[] getJac;
#endif
	}

	//	============================================================================
	//	******************************< Functions >*********************************
	//	============================================================================
	int GetSingular(void) { return singular; }
	int IterationCounter(void) { return numFunctionsTotal; }
	int GetCalculationState(void)
	{
		if (fi.MaxAbs() > 1.e-3 && nlsCalculationState != DUBIOUS_STATE)
			return 7;
		else
			return nlsCalculationState;
	}
	double GetMaximumResidual(void)
	{
		return fi.MaxAbs();
	}
	int NumNewtons(void) { return numNewtons; }
	int NumQuasiNewtons(void) { return numQuasiNewtons; }
	int NumGradients(void) { return numGradients; }
	int NumFunctionsForGradient(void) { return numFunctionsForGradient; }
	int NumFunctionsForNumericalJacobian(void) { return numNumericalJacobians * numVariables; }
	int NumMonodimensional(void) { return numMonodimensional; }
	int NumAnalyticalJacobians(void) { return numAnalyticalJacobians; }
	int NumNumericalJacobians(void) { return numNumericalJacobians; }
	int NumGaussFactorizations(void) { return numGaussFactorizations; }
	int NumLQFactorizations(void) { return numLQFactorizations; }
	int NumLinearSolutions(void) { return numLinearSolutions; }

	virtual void ObjectBzzPrint(void);
	void SetConstraints(BzzVector& xmi, BzzVector& xma);
	void SetMinimumConstraints(BzzVector& yMi);
	void SetMaximumConstraints(BzzVector& yMa);
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
	void SetTolerance(double tA, double tR);
	void SetTolerance(double tA, BzzVector& tR);
	void SetTolerance(BzzVector& tA, double tR);
	void SetTolerance(BzzVector& tA, BzzVector& tR);
	void SetObjectiveFunctionToleranceInStartingPoint(double ta)
	{
		phiWAbsoluteTolleranceInStartingPoint = ta;
	}
	void SetWeights(BzzVector& w);
	void Restart(BzzVector& x00);
	void Restart(BzzVector& x00,
		void (*ptrS)(BzzVector& xx, BzzVector& ff));
	void Restart(BzzVector& x00,
		void (*ptrS)(BzzVector& xx, BzzVector& ff),
		void (*ptrJ)(BzzVector& xx, BzzMatrix& JJ));
	void Restart(BzzVector& x00,
		void (*ptrS)(BzzVector& xx, BzzVector& ff),
		BzzMatrix& jac, char stopR);
	void GetSolution(BzzVector* xx, BzzVector* ff,
		double* pW = 0, double* pN = 0);
	void GetWeights(BzzVector* w)
	{
		*w = weights;
	}

	char operator ()(void) { return (*this)(DEFAULT_MAX_FUNCTIONS_IN_THE_CALL); }
	char operator ()(int nS);
	char HyperSolve(int ni);
	void operator ()(BzzVector& x00);
	void operator ()(BzzVector& x00,
		void (*ptrS)(BzzVector& xx, BzzVector& ff));
	void SetMaxNewtonCallsNumber(int s)
	{
		maxNewtonCallsNumber = s;
	}
	void SetStopForBadConvergence(char s)
	{
		if (s > 0)stopForBadConvergence = s;
	}
};

#endif // NON_LINEAR_SYSTEM_DOUBLE_HPP