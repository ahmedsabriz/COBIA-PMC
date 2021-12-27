// BZZMATH: Release 7.0

//	===============================< BzzDae.hpp >===============================
//	* BzzDaeMultiValue class: base class for DAE solution with				*
//	* multi - value algorithms in double precision										*
//	* BzzDae class: Gear algorithm for stiff problems							*
//	* BzzDaeSparse class: Gear algorithm for stiff problems					*
//	*											and sparse Jacobian								*
//	* Description:																		 			*
//	*					Metodi Numerici e Software in C++ Vol. II							*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathAdvanced\Dae\Dae\Dae.cpp	*
// * Examples: BzzMath\Examples\BzzMathAdvanced\Dae\DaeSparse\				*
// *				Dae.cpp.cpp																*
// * Tests: BzzMath\Examples\BzzMathAdvanced\Dae\DaeTests\					*
// *			DaeTests.cpp																*
// * Tests: BzzMath\Examples\BzzMathAdvanced\Dae\DaeSparseTests\			*
// *			DaeSparseTests.cpp														*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-1999	Date Written.
//	11-1999	Added SetTasksPrint, SetSubTasksPrint.

/////// Release 4.0
//	02-2000	Added GetInitAndEndTimeStep function.
//	02-2001	Bug fixed in constraints control.

/////// Release 5.0
//	08-2003	Added operator() equal to SetInitialConditions.
//	09-2003	Added operator() equal to SetInitialConditions for sparse Jacobian.
//	05-2008	Added GetTimeInMeshPoint function.
//	05-2008	Added GetYInMeshPoint function.
//	05-2008	Added GetY1InMeshPoint function.

////////////////// Release 6.0
// 07-2009	Added OpenMP parallelization.
// 02-2010	Added StepDebug and StopStepDebug functions.
//	05-2010	Bug fixed in Linear System Solution.
// 07-2012	Added GetInitialConditions function.

//	============================================================================
//	* These classes allow to solve initial value systems of first order DAEs:	*
//	*					dy/dt = f(y,t)																*
//	*					0		  = f(y,t)															*
//	****************************************************************************
//	****** Constructors for BzzDae class											*
//	* TODO: BzzDae o(y0,t0,A,Sys); // Numerical Jacobian						*
//	* BzzDae o(y0,t0,iDer,Sys); // Numerical Jacobian							*
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
//	* A: BzzMatrixSparse																	*
//	* iDer: BzzVectorInt = 1 if exist y'													*
//	****************************************************************************
//	****** Constructors for BzzDaeSparse class									*
//	* BzzDaeSparse o(y0,t0,A,Sys,&Je); // Numerical Jacobian					*
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
// * BzzMatrixCoefficientsExistence Je: coefficients in Jacobian (see below)	*
//	* A: BzzMatrixSparse																	*
//	* iDer: BzzVectorInt = 1 if exist y'													*
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
// * In the class	BzzDaeSparse he must use the constructor					*
//	* BzzDaeSparse o(y0,t0,Sys,&Je);													*
// * In this case it is necessary to define the											*
// * BzzMatrixCoefficientsExistence Je that contains the locations of the		*
// * variables in the functions.																*
//	* For example given the system:															*
// * f[1] = f1(y1,y5,y9);																		*
// * f[2] = f2(y3,y10);																			*
// * ..........																					*
//	* BzzMatrixCoefficientsExistence Je(size,size);										*
//	* Je(1,1),Je(1,5),Je(1,9),Je(2,3),Je(2,10),.....									*
//	****************************************************************************
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
//	* (by overshooting and interpolating)													*
//	****************************************************************************
//	***** 2. Computation without overshooting tCritic									*
//	* y = o(tOut,tCritic);																		*
//	* Output values of y(t) at t = tOut but without overshooting tCritic			*
//	* tCritic may be equal to or beyond tOut in the direction of integration	*
//	* tCritic may not take a value between t0 and tOut									*
//	* If the solver reaches tCritic uses tCritic as last mesh point				*
//	* This option is useful if the problem has a singularity beyond t = tCritic*
// * NOTE:																							*
// * - The solver can evaluate the functions f in t = tCritic						*
//	****************************************************************************
//	****************************************************************************
//	* Here following is a list of the optionals											*
//	****************************************************************************
//	* The step size used on first step is automatically determined					*
// * by the solver.																				*
//	* The user can modify the step size (double h0) on the first step				*
// * by using the following function.														*
//	*		o.SetH0(h0);																			*
// * NOTE:																							*
//	* - It is possible to set the step size only at start.							*
//	****************************************************************************
// * The tollerances tollAbs and tollRel are used by the code in a local		*
// * error test at each step which require:												*
// * 		abs(local error in y(i)) < tollRel(i)*abs(y(i)) + tollAbs(i)			*
// * The default values of tollRel and tollAbs are:									*
// * tollRel = 100.*MachEpsFloat(); tollAbs = 1.e-10;									*
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
// * The default maximum order is 5 for stiff and 12 for non stiff problems.	*
// * The user can reduce these values by using the following function:			*
//	*		o.SetMaxOrder(maxOrder);															*
//	****************************************************************************
//	*		o.SetMaxStep(maxStep);																*
//	****************************************************************************
// * If you know that some y can not assume values smaller or bigger to			*
// * an assigned constant (for example the solution will always be non			*
// * negative) it may help to set a constraint by using the functions:			*
// *		o.SetMinimumConstraints(&yMin);													*
// *		o.SetMaximumConstraints(&yMax);													*
// * or:																								*
// *		o.SetMinimumConstraints(yMin);													*
// *		o.SetMaximumConstraints(yMax);													*
// * yMin and yMax are BzzVector. 													*
// * If you pass in the previuos arguments a pointer the solver swaps the		*
// * vectors yMin or yMax with vectors dimensioned zero.								*
// * If you	pass a vector the solver creates a copy of yMin or yMax and			*
// * uses this	copy.																				*
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
// * int numAnalyticalJacobian = GetNumAnalyticalJacobian();						*
// * int numNumericalJacobian = GetNumNumericalJacobian();							*
// * int numFactorization = GetNumFactorization();										*
// * int numSolution = GetNumSolution();													*
// * double hUsedInPreviousStep = GetHUsed();											*
// * double hInNextStep GetHInNextStep();													*
// * int orderUsed = GetOrderUsed();														*
// * int orderInNextStep = GetOrderInNextStep();										*
// * int daeCalculationState = GetCalculationState();									*
// * int componentWithLargestError = GetComponentWithLargestError();				*
// * BzzVector estimatedErrors = GetEstimatedErrors();						*
//	****************************************************************************

#ifndef DAEMULTI_HPP
#define DAEMULTI_HPP

//	============================================================================
//	======================< class BzzDaeJacobian >==============================
//	============================================================================
class BzzDaeJacobian
{
private:
	double ZERO_DER, ETA2, BETA;
	double hJ, hJf, yh, hInv, hf;
	BzzVector aaa, bbb;
public:
	void (*sysDiff)(BzzVector& y, double t, BzzVector& f);
	void GetJacobianNoConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
	void GetJacobianConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& yMax, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
};

class BzzDaeMultiValue : public BzzBaseClass
{
	friend void bzzNlsDae(BzzVector& x, BzzVector& g);
	friend class BzzMyNonLinearSystemObjectDae;
	friend class BzzMyNonLinearSystemObjectSparseDae;
protected:
	BzzMatrix J, JJ;
	BzzMatrixSparse JS, JJS;
	BzzFactorizedSparseUnspecified S;
	BzzMatrixBand JB, JJB;
	BzzMatricesExistence Jm;
	BzzMatrixBlocks Bm;
	BzzMatrixTridiagonalBlocks Tm;
	//	BzzMatrixBlocks BBm;
	BzzMatrixDiagonalBlocks A11, A12, A21, A22;
	BzzMatrixDiagonalBlocks M12, M21, M22;
	BzzMatrixTridiagonalBlocks T22;
	BzzMatrixBand B22;
	BzzMatrixTridiagonalBlocks T11;
	BzzMatrixBand B11;
	BzzMatrixSparse S12, S21, MS12, MS21;
	BzzMatrix D22, MD22;
	//	BzzMatrixSparseLockedByRows SLh;
	//	BzzNonLinearSystemObject nlsDense;

	enum BzzDaeSparseJacobyanType
	{
		FULL,
		SPARSE,
		BAND,
		BLOCK,
		TRIDIAGONAL_BLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_DIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_TRIDIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_BAND,
		FOUR_A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,
		FOUR_A11_FACTORIZED_BAND_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE
	}daeJacobianType;

	enum BzzDaeCalculationState
	{
		INITIALIZATION_STATE = 1,
		CONTINUATION_STATE = 2,
		EXCESSIVE_WORK_STATE = 12,
		TOO_MUCH_ACCURACY_REQUESTED_STATE = -2,
		ILLEGAL_VALUE_OF_TOUT_STATE = -3,
		REPEATED_ERROR_TEST_FAILURE_STATE = -4,
		CONVERGENCE_TEST_FAILURE_STATE = -5,
		H_EQUAL_HMIN_STATE = -6,
		YOU_MUST_USE_TCRITIC_STATE = -7,
		TCRITIC_IS_BEHIND_TOUT = -8,
		ILLEGAL_CONSTRAINTS = -9,
		EXCEPTION_HANDLING_STOP = -10,
		DEINITIALIZE_STATE = -11,
		YOU_CANNOT_OVERSHOOT_TCRITIC = -12
	}daeCalculationState;

	enum BzzDaeTollerance
	{
		SCALAR_A_SCALAR_R,
		ARRAY_A_SCALAR_R,
		SCALAR_A_ARRAY_R,
		ARRAY_A_ARRAY_R
	}daeTollerance;

	enum BzzDaeTypeOfComputation
	{
		NORMAL_WITH_OVERSHOOTING,
		WITH_TCRITIC
	}daeTypeOfComputation;

	enum BzzDaeDirection
	{
		POSITIVE_DIRECTION,
		NEGATIVE_DIRECTION
	}daeDirection;

	enum BzzDaeHState
	{
		H_DECREASED,
		H_CONST,
		H_INCREASED
	}daeHState;

	enum BzzDaeOrderState
	{
		ORDER_DECREASED,
		ORDER_CONST,
		ORDER_INCREASED
	}daeOrderState;

	enum BzzDaeJacobianState
	{
		JAC_HAS_TO_BE_CHANGED,
		JAC_MODIFIED,
		JAC_OK
	}daeJacobianState;

	enum BzzDaeFactorizationState
	{
		MATRIX_HAS_TO_BE_FACTORIZED,
		MATRIX_FACTORIZED,
	}daeFactorizationState;

	enum BzzDaeConvergenceState
	{
		CONVERGENCE_FAILURE,
		CONVERGENCE_OK
	}daeConvergenceState;

	enum BzzDaeErrorState
	{
		ERROR_FAILURE,
		ERROR_OK
	}daeErrorState;

	enum BzzDaeConstraints
	{
		NO_CONSTRAINTS,
		MINIMUM_CONSTRAINTS,
		MAXIMUM_CONSTRAINTS,
		MIN_MAX_CONSTRAINTS,
	}daeConstraintsState;

	enum BzzDaeTypeDae
	{
		MATRIX_A,
		VECTOR_IDER
	}daeTypeDae;

	static const char* const BZZ_ERROR;
	static const double 	DEFAULT_TOL_ABS,
		DEFAULT_TOL_REL,
		DEFAULT_HSCALE_MAX1,
		DEFAULT_HSCALE_MAX2,
		DEFAULT_HSCALE_MAX3;

	static const unsigned int 	DEFAULT_MAX_STEP;

	static const int 	MAX_CONVERGENCE_FAILURE,
		MAX_ERROR_FAILURE;

	FILE* bzzFileSave;
	FILE* bzzFileDebug;

	double 	h0,
		h,
		hMax,
		hMin,
		hUsedInPreviousStep,
		hInNextStep,
		hUsedInNordsieck,
		hMaxUsed,
		hMinUsed,
		hScale,
		hScaleMax,
		hScalePM1,
		hScaleP,
		hScalePP1,
		hr0,
		t0,
		t,
		tInMeshPoint,
		tOut,
		tCritic,
		tStabilize,
		tollA,
		tollR,
		tollSafe,
		correctionNorm,
		errorNorm,
		sqrtInvSize,
		parOrderPM1,
		//parOrderP,
		parOrderPP1,
		* E,
		* parOrderP,
		timeSystem,
		startDebug,
		endDebug,
		timeFactorization;

	int	numVariables,
		maxOrder,
		minOrder,
		maxOrderUsed,
		//maxStep,
		maxConvergenceIterations,
		iterOrder,
		orderUsed,
		orderInNextStep,
		componentWithLargestError,
		iterConvergence,
		iterConvergenceFailure,
		iterErrorFailure,
		iterMaxOrder,
		numThread,
		stabilize;

	BzzDaeJacobian* getJac;
	//	virtual void SolveStartingNonLinearSystem(void) = 0;
	virtual int SolveStartingNonLinearSystem(void);
	virtual int SolveNonLinearSystem(void) = 0;
	virtual int SolveStabilizeNonLinearSystem(void) = 0;

	unsigned int	numStep,
		numFunction,
		numFunctionForJacobian,
		numNumericalJacobian,
		numAnalyticalJacobian,
		numFactorization,
		numSolution,
		numHIncreased,
		numHDecreased,
		numConvergenceFailure,
		numConvergenceFailureForOrderMax,
		numConvergenceSuccess,
		numErrorCheckFailure,
		numErrorCheckSuccess,
		numNoChange,
		numHConst,
		numOrderMin,
		maxStep,
		numNonLinearSystem,
		numNonLinearSystemSuccess,
		numNonLinearSystemFailure;

	char	fixedOrderMin,
		saveResults,
		saveDebug,
		first,
		printTasks,
		printSubTasks;

	BzzVectorInt varSaved;

	double tInMinH;
	BzzVector yInMinH;

	BzzVector	y0,
		y, yy,
		yMin, yMax,
		f,
		tollAbs,
		tollRel,
		errorWeightBzzVector,
		errorWeightSystemBzzVector,
		* z,
		* v,
		* r,
		b,
		db,
		va, vb, vc;

	BzzMatrixSparse A;
	BzzVectorInt iDer;

	void (*sysDiff)(BzzVector& y, double t, BzzVector& f);
	void (*stepPrintOut)(BzzVector& y, double t);

	void Initialize(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void BzzDaeSolve(void);
	void DriverStep(void);
	void ErrorWeightBzzVector(BzzVector& w);
	void ErrorWeightSystemBzzVector(BzzVector& w);
	double ErrorControl(BzzVector& w);
	double ErrorSystemControl(BzzVector& w);
	void RescaleTollerance(void);
	void MultiValueStep(void);
	void NewOrderNewH(void);
	void SetDirection(void);
	void ControlDirection(void);
	void InitialStepSize(void);
	void Interpolation(void);
	void PerformsOneStep(void);
	virtual void FindCorrection(void) = 0;
	virtual void WhatToDoInCaseOfConvergenceFailure(void) = 0;
	virtual void ConvergenceRate(void) = 0;
	void ObjectBzzPrint(void);
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void SetInitialConditions(BzzVector& y00, double t00);

public:
	double tt;
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default
	BzzDaeMultiValue(void);

	BzzDaeMultiValue(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	void Deinitialize(void);
	~BzzDaeMultiValue(void);

	//	============================================================================
	//	==============================< Functions > ================================
	//	============================================================================
	double GetSolutionInStopPoint(BzzVector* ys)
	{
		if (daeCalculationState == H_EQUAL_HMIN_STATE)
		{
			*ys = yInMinH;
			return tInMinH;
		}
		return 0.;
	}

	int GetNumStep(void) { return numStep; }
	int GetNumFunction(void) { return numFunction; }
	int GetNumFunctionForJacobian(void) { return numFunctionForJacobian; }
	int GetNumAnalyticalJacobian(void) { return numAnalyticalJacobian; }
	int GetNumNumericalJacobian(void) { return numNumericalJacobian; }
	int GetNumFactorization(void) { return numFactorization; }
	int GetNumSolution(void) { return numSolution; }
	double GetHUsed(void) { return hUsedInPreviousStep; }
	double GetHInNextStep(void) { return hInNextStep; }
	int GetOrderUsed(void) { return orderUsed; }
	int GetOrderInNextStep(void) { return orderInNextStep; }
	int GetCalculationState(void) { return daeCalculationState; }
	int GetComponentWithLargestError(void);
	BzzVector GetEstimatedErrors(void);
	BzzVector GetTollRel(void);
	BzzVector GetTollAbs(void);
	void GetInitAndEndTimeStep(double* tInitStep, double* tEndStep)
	{
		*tInitStep = t - hUsedInPreviousStep;*tEndStep = t;
	}
	double GetTimeInMeshPoint(void) { return t; }
	BzzVector GetYInMeshPoint(void) { return z[0]; }
	BzzVector GetY1InMeshPoint(void) { Division(z[1], hInNextStep, &va);return va; }
	BzzVector GetInitialConditions(void) { return y0; }

	void BzzPrintErrorState(void);
	void StepPrint(char* sv, const char* pr = "Labels");
	void StepPrint(char* sv, BzzVectorInt& isv, const char* pr = "Labels");
	void StepPrint(void (*stepPrintOut)(BzzVector& y, double t));
	void StepDebug(char* sd, char deb = 1, double st = 0., double et = 0.);
	void StopStepDebug(void)
	{
		if (saveDebug != 0)
		{
			saveDebug = 0;fclose(bzzFileDebug);bzzFileDebug = 0;
		}
	}

	void SetH0(double h00);
	void SetHMin(double hm);
	void SetMaxStep(int maxS);
	void SetTolAbs(double tA);
	void SetTolAbs(BzzVector* tA);
	void SetTolAbs(const BzzVector& tA);
	void SetTolRel(double tR);
	void SetTolRel(BzzVector* tR);
	void SetTolRel(const BzzVector& tR);
	void SetMinimumConstraints(BzzVector* yMi);
	void SetMinimumConstraints(BzzVector& yMi);
	void SetMaximumConstraints(BzzVector* yMa);
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

	BzzVector operator() (double tF);
	// ELIMINATO BzzVector Operator(double tF){return (*this)(tF);}
	BzzVector Operator(double tF);
	//{BzzVector result; result = (*this)(tF); return result;}
	BzzVector operator() (double tF, double tC);
	//ELIMINATO BzzVector Operator(double tF,double tC){return (*this)(tF,tC);}
	BzzVector Operator(double tF, double tC);
	//{BzzVector result; result = (*this)(tF, tC); return result;}
};

class BzzDaeStiffBase : public BzzDaeMultiValue
{
protected:
	static const int 	MAX_ORDER,
		DEFAULT_MAX_CONVERGENCE_ITER,
		MAX_ITERATIONS_FACTORIZATION;
	static const unsigned int	MAX_ITERATIONS_JACOBIAN;

	static const double	PAR_ORDER_PM1,
		PAR_ORDER_P,
		PAR_ORDER_PP1;

	enum JacobianType
	{
		JAC_CONST,	 // BzzDae(y0,t0,ptrSys,J);
		JAC_ANALYTIC, // BzzDae(y0,t0,ptrSys,ptrJac);
		JAC_NUMERICAL // BzzDae(y0,t0,ptrSys);
	}jacType;

	unsigned int	numStepFact,
		numStepJac;

	//	BzzMatrixSparse A;
	//	BzzVectorInt iDer;

	int	maxIterationsJacobian,
		iterConvergenceRate;

	virtual void FindCorrection(void);
	virtual void WhatToDoInCaseOfConvergenceFailure(void);
	virtual void ConvergenceRate(void) = 0;
	virtual void Jacobian(void) = 0;
	virtual void AnalyticalJacobian(void) = 0;
	virtual void BuildBzzMatrixG(void) = 0;
	virtual void SolveLinearSystem(void) = 0;
	//	virtual int SolveStartingNonLinearSystem(void) = 0;
	virtual int SolveNonLinearSystem(void) = 0;
	virtual int SolveStabilizeNonLinearSystem(void) = 0;
	void SetParameters(void);

public:
	void SetMaxOrder(int maxO);

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzDaeStiffBase(void)
		: BzzDaeMultiValue() {};

	BzzDaeStiffBase(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));

	BzzDaeStiffBase(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
};

//	============================================================================
//	*********************< BzzDae for stiff problems >******************
//	============================================================================

class BzzDae : public BzzDaeStiffBase
{
private:
	BzzFactorizedGauss G;
	BzzVector dy, df;

	char weMustChangeJ,
		changedJ;

	void (*sysJac)(BzzVector& yy, double tt, BzzMatrix& JJ);
	//	virtual void FindCorrection(void);
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void);
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);

	friend class BzzMyNonLinearSystemObjectDae;
	virtual int SolveNonLinearSystem(void);
	virtual int SolveStabilizeNonLinearSystem(void)
	{
		BzzError("Stabilize Non Linear System");
		return 1;
	}

public:
	BzzVector operator () (double tF);
	//		{
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		// ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzDaeMultiValue::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
	}

	//	void SetInitialConditions(BzzVector &y00,double t00,
	//		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		 BzzMatrix &JJ);

	//	void SetInitialConditions(BzzVector &y00,double t00,
	//		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		 void (*ptrJac)(BzzVector &yy,double tt,BzzMatrix &JJ));

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzDae(void)
		:BzzDaeStiffBase() {};

	BzzDae(BzzVector& y00, double t00, BzzMatrix& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrix& JJ);

	BzzDae(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, AA, ptrSys);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, AA, ptrSys);
	}

	BzzDae(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrJac)(BzzVector& yy, double tt, BzzMatrix& JJ));

	BzzDae(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, iDer, ptrSys);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, iDer, ptrSys);
	}
};

//	============================================================================
//	***********< BzzDaeSparse for stiff and sparse problems >********
//	============================================================================

class BzzDaeSparse : public BzzDaeStiffBase
{
private:
	BzzMatrixCoefficientsExistence Je;
	BzzFactorizedGauss G;
	BzzFactorizedBandGauss B;
	BzzFactorizedMatrixBlocksGauss Fm;
	BzzFactorizedTridiagonalBlocksGauss Tf;
	BzzVector d;

	int lowerBand, upperBand;
	int blockDimensions, numDiagonalMatrices;

	// four blocks
	int	n11, // A11 dimensions
		n22, // A22 dimensions
		nd1, // block dimensions in A11
		nd2; // block dimensions in A22
	BzzFactorizedDiagonalBlocksGauss FA11;
	// Four TridiagonalGauss-Sparse-Sparse-Dense
	// Four BandGauss-Sparse-Sparse-Dense
	int nS12; // numColumns S12 and numRows S21
	int numBlocks;
	BzzMatrixCoefficientsExistence Je12, Je21;

	BzzFactorizedFourBlocksGauss Ffour;
	void ProductFourBlocksDiagonalTridiagonal(void);
	void ProductFourBlocksDiagonalBand(void);
	void ProductFourBlocksDiagonalDiagonal(void);
	void ProductFourBlocksTridiagonalSparseSparseDense(void);
	void ProductFourBlocksBandSparseSparseDense(void);
	BzzVector ax1, ax2, bx1, bx2, cx1, cx2;

	void (*sysJac)(BzzVector& yy, double tt, BzzMatrixSparse& JJ);
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void);
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);
	virtual int SolveNonLinearSystem(void);
	virtual int SolveStabilizeNonLinearSystem(void)
	{
		BzzError("Stabilize Non Linear System sparse DAE");
		return 1;
	}

public:
	BzzVector operator () (double tF);
	//		{
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		// ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzDaeMultiValue::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
	}

	// Jacobian const
	//void SetInitialConditions(BzzVector &y00,
	//		double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		BzzMatrixSparse &JJJ);
	// Jacobian analytic
	//void SetInitialConditions(BzzVector &y00,double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		void (*ptrJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ));

	// Jacobian Block Tridiagonal

	// Jacobian banded

	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
	// Generic
	// (....,Je);

	// Generic block
	// (....,Jb);

	// Tridiagonal Block
	// (....,int dimBlock);

	// BandGeneric
	// (....,int lowerBand,int upperBand);

	// Four Blocks Diag,diag,diag,tridiag
	// (....,int numRows,int dimBlock1,int dimBlock2);

	// Four Blocks Diag,diag,diag,diag
	// (....,int numRows11,int numRows22,int dimBlock1,int dimBlock2);

	// Four Blocks Diag,diag,diag,band
	// (....,numRows,int dimBlock1,int dimBlock2,int lowerBand,int upperBand);

	// Four Blocks TridagonalBlock,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);

	// Four Blocks Band,sparse,sparse,dense
	// (....,int lB,uB,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);

	void operator()(BzzVector& y00, double t00)
	{
		SetInitialConditions(y00, t00);
	}
	BzzDaeSparse(void)
		:BzzDaeStiffBase() {};

	// Jacobian numerical
		// sparse generic
	///////////////////////////////////////////////////////////////////////////////
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, JE);
	}
	// TODO
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, JE);
	}
	///////////////////////////////////////////////////////////////////////////////
	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, JE);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, JE);
	}
	///////////////////////////////////////////////////////////////////////////////
	// Jacobian Block Matrices
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, JBm);
	}
	// TODO
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, JBm);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, JBm);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, JBm);
	}
	///////////////////////////////////////////////////////////////////////////////
	// block tridiagonal
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, dim);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, dim);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, dim);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, dim);
	}
	///////////////////////////////////////////////////////////////////////////////
	// band
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, low, up);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, low, up);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal tridiagonal
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, nV, nb1, nb2);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, nV, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal diagonal
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, n1, n2, nb1, nb2);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, n1, n2, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal band
	BzzDaeSparse(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrSys, nV, nb1, nb2, low, up);
	}

	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, nV, nb1, nb2, low, up);
	}

	// Four TridiagonalBlocksGauss,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, dimBlock, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, dimBlock, J12, J21);
	}

	// Four BandGauss,sparse,sparse,dense
	// (....,int lB,int uB,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzDaeSparse(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, lB, uB, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrSys, lB, uB, J12, J21);
	}

	void PrintSelectedMethod(void);
};

//	============================================================================
//	********************< BzzDaeDiagonalBlocks >***************************
//	============================================================================
class BzzDaeDiagonalBlocks : public BzzDaeStiffBase
{
private:
	int	numMat,
		numVarBlock;
	//	BzzVector d;
	BzzMatrix* E, EE;
	BzzFactorizedDiagonalBlocksAndSparse EF;
	BzzMatrixSparse* ES, SS;
	BzzVector en;
	void (*sysDiff2)(BzzVector& y, double t, BzzVector& f);

	//	int lowerBand,upperBand;

	//	void (*sysJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ);
	//	virtual void ConvergenceRate(void){}
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void) {};
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);
	//	virtual void SolveStartingNonLinearSystem(void)
	//		{
	//		BzzError("Starting Non Linear system diagonal block DAE");
	//		}

	virtual int SolveNonLinearSystem(void)
	{
		BzzError("Generic Non Linear System diagonal block DAE");
		return 1;
	}

	virtual int SolveStabilizeNonLinearSystem(void)
	{
		BzzError("Stabilize Non Linear System diagonal block ");
		return 1;
	}

public:
	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzDaeMultiValue::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
	}

	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrSys2)(BzzVector& y, double t, BzzVector& f),
		int numB, BzzMatrixSparse* SSS);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrSys2)(BzzVector& y, double t, BzzVector& f),
		int numB, BzzMatrixSparse* SSS)
	{
		SetInitialConditions(y00, t00, ptrSys, ptrSys2, numB, SSS);
	}

	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
	BzzDaeDiagonalBlocks(void)
		:BzzDaeStiffBase() {};

	BzzDaeDiagonalBlocks(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrSys2)(BzzVector& y, double t, BzzVector& f),
		int numB, BzzMatrixSparse* SSS);
	BzzDaeDiagonalBlocks(BzzVector& y00, double t00, BzzVectorInt& iDer,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrSys2)(BzzVector& y, double t, BzzVector& f),
		int numB, BzzMatrixSparse* SSS);
};

//=============================================================================================
class BzzMyNonLinearSystemObjectDae : public BzzMyNonLinearSystemObject
{
private:
	BzzDaeMultiValue* ptrNlsDae;
	BzzVector vb;
public:
	BzzMyNonLinearSystemObjectDae(void) {};
	BzzMyNonLinearSystemObjectDae(BzzDaeMultiValue* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	void operator()(BzzDaeMultiValue* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

//=============================================================================================
class BzzMyNonLinearSystemObjectSparseDae : public BzzMyNonLinearSystemSparseObject
{
private:
	BzzDaeMultiValue* ptrNlsDae;
	BzzVector vb;
public:
	BzzMyNonLinearSystemObjectSparseDae(void) {}
	BzzMyNonLinearSystemObjectSparseDae(BzzDaeMultiValue* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	void operator()(BzzDaeMultiValue* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

#endif // DAEMULTI_HPP