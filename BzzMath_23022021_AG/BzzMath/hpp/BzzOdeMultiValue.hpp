// BZZMATH: Release 7.0

//	==========================< BzzOdeMultiValue.hpp >==========================
//	* BzzOdeMultiValue class: base class for ODE solution with				*
//	* multi - value algorithms in double precision										*
//	* BzzOdeNonStiff class: Adams - Moulton algorithm							*
//	* BzzOdeStiff class: Gear algorithm for stiff problems					*
//	* BzzOdeSparseStiff class: Gear algorithm for stiff problems			*
//	*											and sparse Jacobian								*
//	* Description:																		 			*
//	*					Metodi Numerici e Software in C++ (Capitoli 29, 30)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathAdvanced\Ode\OdeNonStiff\			*
// *				OdeNonStiff.hpp														*
// * Examples: BzzMath\Examples\BzzMathAdvanced\Ode\OdeStiff\				*
// *				OdeStiff.hpp															*
// * Tests: BzzMath\Examples\BzzMathAdvanced\Ode\OdeNonStiffTests\		*
// *				OdeNonStiffTests.hpp													*
// * Tests: BzzMath\Examples\BzzMathAdvanced\Ode\OdeStiffTests\			*
// *				OdeStiffTests.hpp														*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-1994	Date Written
//	12-1994	Added BzzOdeSparseStiff class.
//	01-1995	Added constraints.
//	02-1995	Added the Bzz prefix to the names of the classes.
//	11-1996	Added the StepPrint functions.
//	05-1997	Added the SetInitialConditions and Deinitialize functions.
//	05-1997	Added the BzzExceptionHandling.
//	11-1999	Added SetTasksPrint, SetSubTasksPrint.
/////// Release 4.0
//	01-2000	Added control for time consuming factorization.
//	02-2000	Added GetInitAndEndTimeStep function.
//	02-2001	Bug fixed in constraints control.

/////// Release 5.0
//	08-2003	Added operator() equal to SetInitialConditions.
//	09-2003	Added operator() equal to SetInitialConditions for sparse Jacobian.
//	12-2004	Added GetOdeCalculationState function.
//	12-2004	Added GetTimeInMeshPoint function.
//	12-2004	Added GetYInMeshPoint function.
//	12-2004	Added GetY1InMeshPoint function.
//	01-2005	StopIntegrationWhenSumAbsY1IsLessThan.

////////////////// Release 6.0
// 07-2009	Added OpenMP parallelization.
// 02-2010	Added StepDebug and StopStepDebug functions.

//	============================================================================
//	* These classes allow to solve initial value systems of first order ODEs:	*
//	*					dy/dt = f(y,t)																*
//	****************************************************************************
//	****** Constructors for BzzOdeNonStiff class:								*
//	* BzzOdeNonStiff o(y0,t0,Sys);													*
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
//	****************************************************************************
//	****** Constructors for BzzOdeStiff class										*
//	* BzzOdeStiff o(y0,t0,Sys); // Numerical Jacobian							*
//	* BzzOdeStiff o(y0,t0,Sys,Jac); // Analytical Jacobian					*
//	* BzzOdeStiff o(y0,t0,Sys,J);	// Constant Jacobian							*
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
//	* Jac: Name of function of the analytical Jacobian (see below)					*
//	* BzzMatrix J: constant Jacobian													*
//	****************************************************************************
//	****** Constructors for BzzOdeSparseStiff class								*
//	* BzzOdeSparseStiff o(y0,t0,Sys,&Je); // Numerical Jacobian				*
//	* The arguments are as following:														*
//	* BzzVector y0: the vector of initial values									*
//	* double t0: the initial value of the independent variable						*
//	* Sys: Name of function for right-hand side vector f (see below)				*
// * BzzMatrixCoefficientsExistence Je: coefficients in Jacobian (see below)	*
//	* Jac: Name of function of the analytical Jacobian (see below)					*
//	* BzzMatrixSparse J: constant Jacobian											*
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
// * In the class BzzOdeStiff it is possible to use the constructor		*
//	* BzzOdeStiff o(y0,t0,Sys,Jac); 													*
// * that requires the analytical Jacobian. In this case it is necessary		*
// * to declare the function in which the Jacobian is calculated.					*
//	* The prototype of the function is:														*
//	*		void <name>(BzzVector &y,double t,BzzMatrix &J);			*
//	* In the definition of this function the components J[i][j] of					*
//	* the Jacobian as function of the scalar t and vector								*
//	* y[1], y[2], ..., y[N] are calculated.												*
//	* Example:																						*
//	* void MyJacobian(BzzVector &y,double t,BzzMatrix &J)				*
//	*	{																								*
//	*	J[1][1] = t;																				*
//	*	J[1][2] = 1.;																				*
//	*	J[2][1] = -2.*y[2];																		*
//	*	J[2][2] = -2.*y[1];																		*
//	*	}																								*
// * NOTE:																							*
// * - The user must not modify the vector y.											*
// * - The coefficients of the matrix J are accessible using [][] or (,).		*
// * - [][] is more efficient but the user is responsible of the correct		*
// * 		range.																					*
// * - The coefficient J[i][j] or J(i,j) is the derivative of function i		*
// *		respect variable j.																	*
//	****************************************************************************
// * If the user wants a numerical Jacobian in the class								*
// * BzzOdeSparseStiff he must use the constructor								*
//	* BzzOdeSparseStiff o(y0,t0,Sys,&Je);											*
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
// * The tolerances tollAbs and tollRel are used by the code in a local			*
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
// * int numStep = o.GetNumStep();															*
// * int numFunction = o.GetNumFunction();												*
// * int numAnalyticalJacobian = o.GetNumAnalyticalJacobian();						*
// * int numNumericalJacobian = o.GetNumNumericalJacobian();						*
// * int numFactorization = o.getNumFactorization();									*
// * int numSolution = o.GetNumSolution();												*
// * double hUsedInPreviousStep = o.GetHUsed();											*
// * double hInNextStep = o.GetHInNextStep();											*
// * int orderUsed = o.GetOrderUsed();														*
// * int orderInNextStep = o.GetOrderInNextStep();										*
// * int odeCalculationState = o.GetCalculationState();								*
// * int componentWithLargestError = o.GetComponentWithLargestError();			*
// * BzzVector estimatedErrors = o.GetEstimatedErrors();						*
// * double tInMeshPoint =	o.GetTimeInMeshPoint()										*
// * BzzVector yInMeshPoint =	o.GetYInMeshPoint()								*
// * BzzVector y1InMeshPoint =	o.GetY1InMeshPoint()							*
// * int status = GetOdeCalculationState();												*
//	****************************************************************************

#ifndef ODEMULTI_HPP
#define ODEMULTI_HPP

//	============================================================================
//	======================< class BzzOdeStiffJacobian >=========================
//	============================================================================
class BzzOdeStiffJacobian
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

class BzzOdeMultiValue : public BzzBaseClass
{
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

	enum BzzOdeSparseJacobyanType
	{
		FULL,
		SPARSE,
		SPARSE_BAND,
		SPARSE_FULL,
		BAND,
		BLOCK,
		TRIDIAGONAL_BLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_DIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_TRIDIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_BAND,
		FOUR_A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,
		FOUR_A11_FACTORIZED_BAND_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_FIRST,
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_SECOND,
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_THIRD,
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_FOURTH,
		JACOBIAN_LOCKED_GENERIC,
		JACOBIAN_LOCKED_BAND
	}odeJacobianType;

	enum BzzOdeCalculationState
	{
		INITIALIZATION_STATE = 1,
		CONTINUATION_STATE = 2,
		INTEGRATION_STOPPED_BEFORE_RECALCULATING_JACOBIAN = 10,
		INTEGRATION_STOPPED_WHEN_SUM_ABS_Y1_IS_LESS_THAN = 11,
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
		YOU_CANNOT_OVERSHOOT_TCRITIC = -12,
	}odeCalculationState;

	enum BzzOdeTollerance
	{
		SCALAR_A_SCALAR_R,
		ARRAY_A_SCALAR_R,
		SCALAR_A_ARRAY_R,
		ARRAY_A_ARRAY_R
	}odeTollerance;

	enum BzzOdeTypeOfComputation
	{
		NORMAL_WITH_OVERSHOOTING,
		WITH_TCRITIC
	}odeTypeOfComputation;

	enum BzzOdeDirection
	{
		POSITIVE_DIRECTION,
		NEGATIVE_DIRECTION
	}odeDirection;

	enum BzzOdeHState
	{
		H_DECREASED,
		H_CONST,
		H_INCREASED
	}odeHState;

	enum BzzOdeOrderState
	{
		ORDER_DECREASED,
		ORDER_CONST,
		ORDER_INCREASED
	}odeOrderState;

	enum BzzOdeJacobianState
	{
		JAC_HAS_TO_BE_CHANGED,
		JAC_MODIFIED,
		JAC_OK
	}odeJacobianState;

	enum BzzOdeFactorizationState
	{
		MATRIX_HAS_TO_BE_FACTORIZED,
		MATRIX_FACTORIZED,
	}odeFactorizationState;

	enum BzzOdeConvergenceState
	{
		CONVERGENCE_FAILURE,
		CONVERGENCE_OK
	}odeConvergenceState;

	enum BzzOdeErrorState
	{
		ERROR_FAILURE,
		ERROR_OK
	}odeErrorState;

	enum BzzOdeConstraints
	{
		NO_CONSTRAINTS,
		MINIMUM_CONSTRAINTS,
		MAXIMUM_CONSTRAINTS,
		MIN_MAX_CONSTRAINTS,
	}odeConstraintsState;

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
		timeFactorization,
		startDebug,
		endDebug,
		sumAbsY1;

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
		printTasks,
		printSubTasks,
		stopIntegrationBeforeRecalcuatingJacobian,
		stopIntegrationWhenSumAbsY1IsLessThan,
		numThread,
		numLockedGroup;

	BzzVectorInt jaux, iaux;
	BzzVector aux;

	BzzOdeStiffJacobian* getJac;
	BzzMatrixSparseLocked AD, AA;

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
		numHConstForLargeFactorizationTime,
		numOrderMin,
		maxStep;

	char	fixedOrderMin,
		saveResults,
		saveDebug,
		first;

	BzzVectorInt varSaved;

	BzzVector	y0,
		y,
		yMin, yMax,
		f,
		tollAbs,
		tollRel,
		errorWeightBzzVector,
		* z,
		* v,
		* r,
		b,
		db,
		va, vb;

	void (*sysDiff)(BzzVector& y, double t, BzzVector& f);
	void (*stepPrintOut)(BzzVector& y, double t);
	// generic new
// TODO
	/////////////////////////////////////////////////////////////////////////////////
	int jacDifferent; // se la funzione in cui calcolare lo Jacobiano è diversa = 1
	int lockedMatrix; // = 1 se usa locked
	void (*ptrFunRCV)(BzzVector& y, BzzVector& f); // Autonomous system
	void (*ptrJacRCV)(BzzVector& y, BzzVector& f); // Autonomous system
	/////////////////////////////////////////////////////////////////////////////////

	void Initialize(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void Initialize(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, BzzVector& f),
		void (*ptrSysJac)(BzzVector& y, BzzVector& f)); // TODO con Jac!!!!!!!!!!!!!!!!!!!!!!1
	void BzzOdeSolve(void);
	void DriverStep(void);
	void ErrorWeightBzzVector(BzzVector& w);
	double ErrorControl(BzzVector& w);
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
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default
	BzzOdeMultiValue(void);

	BzzOdeMultiValue(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	void Deinitialize(void);
	~BzzOdeMultiValue(void);

	//	============================================================================
	//	==============================< Functions > ================================
	//	============================================================================
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
	int GetCalculationState(void) { return odeCalculationState; }
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
	int GetOdeCalculationState(void) { return odeCalculationState; }
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
	void SetHMax(double hm) { hMax = hm; }
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
	void StopIntegrationBeforeRecalcuatingJacobian(int nJ = 1)
	{
		stopIntegrationBeforeRecalcuatingJacobian = nJ;
	}
	void StopIntegrationWhenSumAbsY1IsLessThan(double y1v)
	{
		stopIntegrationWhenSumAbsY1IsLessThan = 1;sumAbsY1 = y1v;
	}
	BzzVector operator() (double tF);
	// ELIMINATO
	//BzzVector Operator(double tF){return (*this)(tF);}
	BzzVector Operator(double tF);//{BzzVector result; result = (*this)(tF); return result;}
	BzzVector operator() (double tF, double tC);
	// ELIMINATO
	//BzzVector Operator(double tF,double tC){return (*this)(tF,tC);}
	BzzVector Operator(double tF, double tC);//{BzzVector result; result = (*this)(tF,tC); return result;}
};

//	============================================================================
//	***********< BzzOdeNonStiff class: Adams Moulton algorithm >*************
//	============================================================================

class BzzOdeNonStiff : public BzzOdeMultiValue
{
private:
	static const int 	MAX_ORDER,
		DEFAULT_MAX_CONVERGENCE_ITER;

	static const double	PAR_ORDER_PM1,
		PAR_ORDER_P,
		PAR_ORDER_PP1;
	virtual void FindCorrection(void);
	virtual void WhatToDoInCaseOfConvergenceFailure(void);
	virtual void ConvergenceRate(void) {};
	void SetParameters(void);

public:
	void SetMaxOrder(int maxO);
	BzzVector operator () (double tF);
	//		{
	//		//ELIMINATOreturn Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		//ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF,tC); return result;
	//		}

	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00, ptrSys);
		SetParameters();
	}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		maxConvergenceIterations = DEFAULT_MAX_CONVERGENCE_ITER;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
	}

	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
	BzzOdeNonStiff(void)
		: BzzOdeMultiValue() {}

	BzzOdeNonStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeMultiValue(y00, t00, ptrSys)
	{
		SetParameters();
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00, ptrSys);
		SetParameters();
	}
};

class BzzOdeStiffBase : public BzzOdeMultiValue
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
		JAC_CONST,	 // BzzOdeStiff(y0,t0,ptrSys,J);
		JAC_ANALYTIC, // BzzOdeStiff(y0,t0,ptrSys,ptrJac);
		JAC_NUMERICAL // BzzOdeStiff(y0,t0,ptrSys);
	}jacType;

	unsigned int	numStepFact,
		numStepJac;

	int	maxIterationsJacobian,
		iterConvergenceRate;

	virtual void FindCorrection(void);
	virtual void WhatToDoInCaseOfConvergenceFailure(void);
	virtual void ConvergenceRate(void) = 0;
	virtual void Jacobian(void) = 0;
	virtual void AnalyticalJacobian(void) = 0;
	virtual void BuildBzzMatrixG(void) = 0;
	virtual void SolveLinearSystem(void) = 0;
	void SetParameters(void);

public:

	void SetMaxOrder(int maxO);

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzOdeStiffBase(void)
		: BzzOdeMultiValue() {};

	BzzOdeStiffBase(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
		: BzzOdeMultiValue(y00, t00, ptrSys)
	{
		SetParameters();
	}
};

//	============================================================================
//	*********************< BzzOdeStiff for stiff problems >******************
//	============================================================================

class BzzOdeStiff : public BzzOdeStiffBase
{
private:
	//	BzzMatrix J,JJ;
	//	BzzVector d,dy,df;
	BzzVector dy, df;
	BzzFactorizedGauss G;

	char weMustChangeJ,
		changedJ;

	void (*sysJac)(BzzVector& yy, double tt, BzzMatrix& JJ);
	//	virtual void FindCorrection(void);
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void);
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);

public:
	BzzVector operator () (double tF);
	//		{
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		//ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00);
		numStepFact = numStepJac = 0;
		iterConvergenceRate = 0;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
		//	SetParameters();
	}

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzOdeStiff(void)
		:BzzOdeStiffBase() {};

	BzzOdeStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrix& JJ)
	{
		SetInitialConditions(y00, t00, ptrSys, JJ);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrix& JJ);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrix& JJ)
	{
		SetInitialConditions(y00, t00, ptrSys, JJ);
	}

	BzzOdeStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrJac)(BzzVector& yy, double tt, BzzMatrix& JJ))
	{
		SetInitialConditions(y00, t00, ptrSys, ptrJac);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrJac)(BzzVector& yy, double tt, BzzMatrix& JJ));
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrJac)(BzzVector& yy, double tt, BzzMatrix& JJ))
	{
		SetInitialConditions(y00, t00, ptrSys, ptrJac);
	}

	BzzOdeStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, ptrSys);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f));
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f))
	{
		SetInitialConditions(y00, t00, ptrSys);
	}
	void GetLastJacobian(BzzMatrix* JJ) { *JJ = J; }
};

//	============================================================================
//	***********< BzzOdeSparseStiff for stiff and sparse problems >********
//	============================================================================

class BzzOdeSparseStiff : public BzzOdeStiffBase
{
private:
	// generic new
// TODO
	/////////////////////////////////////////////////////////////////7
	BzzVectorInt rA, cA, rB, cB;
	BzzVector vA;
	BzzMatrixSparseLocked AL, AG;
	BzzMatrixCoefficientsLocked BL;
	BzzVectorInt newRowsOrder;
	BzzVectorInt newColumnsOrder;

	/////////////////////////////////////////////////////////////////7

	BzzMatrixCoefficientsExistence Je;
	//	BzzMatrixSparse JS,JJS;
	BzzFactorizedSparseUnspecified S;

	//	BzzMatrix J,JJ;
	BzzFactorizedGauss G;

	//	BzzMatrixBand JB,JJB;
	BzzFactorizedBandGauss B;

	BzzMatricesExistence Jm;
	//	BzzMatrixBlocks Bm;
	//	BzzMatrixBlocks BBm;
	BzzFactorizedMatrixBlocksGauss Fm;

	//	BzzMatrixTridiagonalBlocks Tm;
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
	//	BzzMatrixDiagonalBlocks A11,A12,A21,A22;
	//	BzzMatrixDiagonalBlocks M12,M21,M22;
	//	BzzMatrixTridiagonalBlocks T22;
	BzzMatrixBand B22;
	// Four TridiagonalGauss-Sparse-Sparse-Dense
	// Four BandGauss-Sparse-Sparse-Dense
	int nS12; // numColumns S12 and numRows S21
	int numBlocks;
	//	BzzMatrixTridiagonalBlocks T11;
	//	BzzMatrixBand B11;
	//	BzzMatrixSparse S12,S21,MS12,MS21;
	//	BzzMatrix D22,MD22;
	BzzMatrixCoefficientsExistence Je12, Je21;

	BzzFactorizedFourBlocksGauss Ffour;

	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked DaL;
	BzzMatrixSparseLockedByRows SL, SLh;

	void ProductFourBlocksDiagonalTridiagonal(void);
	void ProductFourBlocksDiagonalBand(void);
	void ProductFourBlocksDiagonalDiagonal(void);
	void ProductFourBlocksTridiagonalSparseSparseDense(void);
	void ProductFourBlocksBandSparseSparseDense(void);
	BzzVector ax1, ax2, bx1, bx2, cx1, cx2;

	void (*sysJacSparse)(BzzVector& yy, double tt, BzzMatrixSparse& JJ);
	void (*sysJacBand)(BzzVector& yy, double tt, BzzMatrixBand& JJB);
	void (*sysJacTriDiaBlo)(BzzVector& yy, double tt, BzzMatrixTridiagonalBlocks& JJT);
	void (*sysJacFourD)(BzzVector& yy, double tt,
		BzzMatrixDiagonalBlocks& AA11,
		BzzMatrixDiagonalBlocks& AA12,
		BzzMatrixDiagonalBlocks& AA21,
		BzzMatrixDiagonalBlocks& AA22);
	void (*sysJacFourT)(BzzVector& yy, double tt,
		BzzMatrixDiagonalBlocks& AA11,
		BzzMatrixDiagonalBlocks& AA12,
		BzzMatrixDiagonalBlocks& AA21,
		BzzMatrixTridiagonalBlocks& TT22);
	void (*sysJacFourB)(BzzVector& yy, double tt,
		BzzMatrixDiagonalBlocks& AA11,
		BzzMatrixDiagonalBlocks& AA12,
		BzzMatrixDiagonalBlocks& AA21,
		BzzMatrixBand& BB22);
	void (*sysJacFourTSSD)(BzzVector& yy, double tt,
		BzzMatrixTridiagonalBlocks& TT11,
		BzzMatrixSparse& SS12,
		BzzMatrixSparse& SS21,
		BzzMatrix& DD22);

	void (*sysDiffDiaLock1)(BzzVector& y, double t, BzzVector& f);
	void (*sysJacDiaLock1)(BzzVector& yy, double tt, BzzMatrixDiagonalBlocks& JJD);
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void);
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);

public:
	BzzVector operator () (double tF);
	//		{
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		//return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00);
		numStepFact = numStepJac = 0;
		iterConvergenceRate = 0;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
		//		SetParameters();
	}

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

	BzzOdeSparseStiff(void)
		:BzzOdeStiffBase() {};

	// Jacobian const
	//	BzzOdeSparseStiff(BzzVector &y00,double t00,
	//		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		 BzzMatrixSparse &JJ);
	//void SetInitialConditions(BzzVector &y00,
	//		double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		BzzMatrixSparse &JJJ);
	//void operator()(BzzVector &y00,
	//		double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		BzzMatrixSparse &JJJ)
	//	{
	//	SetInitialConditions(y00,t00,ptrSys,JJJ);
	//	}

	//	BzzOdeSparseStiff(BzzVector &y00,double t00,
	//		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		 void (*ptrJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ));

	// Jacobian analytic
	//void SetInitialConditions(BzzVector &y00,double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		void (*ptrJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ));
	//void operator()(BzzVector &y00,double t00,
	//		void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
	//		void (*ptrJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ))
	//	{
	//	SetInitialConditions(y00,t00,ptrSys,ptrJac);
	//	}

		// sparse generic new
		// sparse generic
	//////////////////////////////////////////////////////////////////////
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, BzzVector& f),
		BzzVectorInt* rr, BzzVectorInt* cc,
		void (*ptrSysJac)(BzzVector& y, BzzVector& f) = 0);

	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, BzzVector& f),
		BzzVectorInt* rr, BzzVectorInt* cc, BzzVector* vv, BzzVectorInt* nl,
		void (*ptrSysJac)(BzzVector& y, BzzVector& f) = 0);

	////////////////////////////////////////////////////////////////////////7
		// sparse generic obsolete
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		//		 		 BzzMatrixCoefficientsExistence &JE);
		BzzMatrixCoefficientsExistence* JE,
		void (*sysJacS)(BzzVector& yy, double tt, BzzMatrixSparse& JJ) = 0)
	{
		SetInitialConditions(y00, t00, ptrSys, JE, sysJacS);
	}
	// Jacobian numerical
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE,
		void (*sysJacS)(BzzVector& yy, double tt, BzzMatrixSparse& JJ) = 0);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatrixCoefficientsExistence* JE,
		void (*sysJacS)(BzzVector& yy, double tt, BzzMatrixSparse& JJ) = 0)
	{
		SetInitialConditions(y00, t00, ptrSys, JE, sysJacS);
	}

	// Jacobian Block Matrices
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, ptrSys, JBm);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, ptrSys, JBm);
	}

	// block tridiagonal
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, ptrSys, dim);
	}
	// Jacobian Block Tridiagonal
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dim)
	{
		SetInitialConditions(y00, t00, ptrSys, dim);
	}

	// band
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up,
		void (*sysJacB)(BzzVector& yy, double tt, BzzMatrixBand& JJB) = 0)
	{
		SetInitialConditions(y00, t00, ptrSys, low, up, sysJacB);
	}
	// Jacobian banded
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up,
		void (*sysJacB)(BzzVector& yy, double tt, BzzMatrixBand& JJB) = 0);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int low, int up,
		void (*sysJacB)(BzzVector& yy, double tt, BzzMatrixBand& JJB) = 0)
	{
		SetInitialConditions(y00, t00, ptrSys, low, up, sysJacB);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal tridiagonal
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrSys, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrSys, nV, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal diagonal
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrSys, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrSys, n1, n2, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal band
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, ptrSys, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, ptrSys, nV, nb1, nb2, low, up);
	}

	// Four TridiagonalBlocksGauss,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrSys, dimBlock, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrSys, dimBlock, J12, J21);
	}

	// Four BandGauss,sparse,sparse,dense
	// (....,int lowerBand,int upperBand,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrSys, lB, uB, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrSys, lB, uB, J12, J21);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear First Typology
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*sysDiffD)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, ptrSys, sysDiffD, dia, QQ);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*sysDiffD)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*sysDiffD)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, ptrSys, sysDiffD, dia, QQ);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear Second Typology
	BzzOdeSparseStiff(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM)
	{
		SetInitialConditions(y00, t00, ptrSys, dia, QQ, fileM);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM);
	void operator()(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM)
	{
		SetInitialConditions(y00, t00, ptrSys, dia, QQ, fileM);
	}
	/*
	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear Third Typology
	BzzOdeSparseStiff(BzzMatrix &y00,double t00,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
				 BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ)
		{
		SetInitialConditions(y00,t00,ptrSys,dia,QQ);
		}
	void SetInitialConditions(BzzMatrix &y00,double t00,
			void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
			BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ);
	void operator()(BzzMatrix &y00,double t00,
			void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
			BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ)
		{
		SetInitialConditions(y00,t00,ptrSys,dia,QQ);
		}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear Fourth Typology
	BzzOdeSparseStiff(BzzMatrix &y00,double t00,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
				 BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ,char *fileM)
		{
		SetInitialConditions(y00,t00,ptrSys,dia,QQ,fileM);
		}
	void SetInitialConditions(BzzMatrix &y00,double t00,
			void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
			BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ,char *fileM);
	void operator()(BzzMatrix &y00,double t00,
			void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
			BzzVectorInt dia,BzzMatrixSparseLockedByRows *QQ,char *fileM)
		{
		SetInitialConditions(y00,t00,ptrSys,dia,QQ,fileM);
		}
	*/
	void PrintSelectedMethod(void);
};

//	============================================================================
//	********************< BzzOdeDiagonalBlocks >***************************
//	============================================================================

class BzzOdeDiagonalBlocks : public BzzOdeStiffBase
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

public:
	BzzVector operator () (double tF);
	//		{
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		//ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValue::SetInitialConditions(y00, t00);
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
	BzzOdeDiagonalBlocks(void)
		:BzzOdeStiffBase() {};

	BzzOdeDiagonalBlocks(BzzVector& y00, double t00,
		void (*ptrSys)(BzzVector& y, double t, BzzVector& f),
		void (*ptrSys2)(BzzVector& y, double t, BzzVector& f),
		int numB, BzzMatrixSparse* SSS);
};

#endif // ODEMULTI_HPP