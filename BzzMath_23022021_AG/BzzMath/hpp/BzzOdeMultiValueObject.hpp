// BZZMATH: Release 7.0

//	==========================< BzzOdeMultiValueObject.hpp >==========================
//	* BzzOdeMultiValueObject class: base class for ODE solution with				*
//	* multi - value algorithms in double precision												*
//	* BzzOdeNonStiffObject class: Adams - Moulton algorithm							*
//	* BzzOdeStiffObject class: Gear algorithm for stiff problems					*
//	* BzzOdeSparseObject class: Gear algorithm for stiff problems					*
//	*											and sparse Jacobian										*
//	==================================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2004	Date Written
//	12-2004	Added GetOdeCalculationState function;
//	12-2004	Added GetTimeInMeshPoint function.
//	12-2004	Added GetYInMeshPoint function.
//	12-2004	Added GetY1InMeshPoint function.
//	01-2005	StopIntegrationWhenSumAbsY1IsLessThan

////////////////// Release 6.0
// 07-2009	Added OpenMP parallelization.
// 02-2010	Added StepDebug and StopStepDebug functions.

//	============================================================================
//	* These classes allow to solve initial value systems of first order ODEs:	*
//	*					dy/dt = f(y,t)																*
//	****************************************************************************

#ifndef ODEMULTI_OBJECT_HPP
#define ODEMULTI_OBJECT_HPP

//	============================================================================
//	=======================< class BzzOdeSystemObject >===================
//	============================================================================
class BzzOdeSystemObject : public BzzBaseClass
{
public:
	// -2: initial
	// -1: normal
	// 0: Jacobian base
	// > 0 and <= numVariables: Jacobian computation for jacobianIndex variable
	// > numVariables Jacobian computation for jacobianVariables variables
	int jacobianIndex;
	BzzVectorInt jacobianVariables;
	BzzMatrixSparseLockedByRows SL;
	BzzMatrix A;
	BzzFactorizedGauss G;

	char fileDiagonal[20];
	char fileDiagonalFactorized[20];

	virtual void GetSystemFunctions(BzzVector& y, double t, BzzVector& f) {};
	virtual void GetSystemFunctionsDBALMemo(BzzVector& y, double t, BzzVector& f) {};
	virtual void ObjectBzzPrint(void) {};
	virtual void GetJacobian(BzzVector& y, double t, BzzMatrix& JJ) {};
	virtual void GetJacobian(BzzVector& yy, double tt, BzzMatrixSparse& JJ) {};
	virtual void GetJacobian(BzzVector& yy, double tt, BzzMatrixBand& JJB) {};
	virtual void GetJacobian(BzzVector& yy, double tt, BzzMatrixTridiagonalBlocks& JJT) {};
	virtual void GetJacobian(BzzVector& yy, double tt) {};
	virtual void GetBuildJacobian(double hr, BzzMatrixSparseLockedByRows* Sh) {};
};

//	============================================================================
//	=================< class BzzOdeStiffObjectJacobian >========================
//	============================================================================
class BzzOdeStiffObjectJacobian
{
private:
	double ZERO_DER, ETA2, BETA;
	double hJ, hJf, yh, hInv, hf;
	BzzVector aaa, bbb;
public:
	BzzOdeSystemObject* ptrObject;
	void GetJacobianNoConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
	void GetJacobianConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& yMax, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
};

class BzzOdeMultiValueObject : public BzzBaseClass
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
	BzzMatrixSparseLockedByRows SLh;

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
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_FOURTH
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

	enum BzzOdeJacobianAssigned
	{
		JACOBIAN_NON_ASSIGNED,
		JACOBIAN_ASSIGNED
	}odeJacobianAssigned;

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
		numThread;

	BzzOdeStiffObjectJacobian* getJac;

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

	//	void (*sysDiff)(BzzVector &y,double t,BzzVector &f);
	BzzOdeSystemObject* ptrObject;
	void (*stepPrintOut)(BzzVector& y, double t);

	void Initialize(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO);
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
		BzzOdeSystemObject* ptrO);
	void SetInitialConditions(BzzVector& y00, double t00);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default
	BzzOdeMultiValueObject(void);

	BzzOdeMultiValueObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	void Deinitialize(void);
	~BzzOdeMultiValueObject(void);

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
	int GetOdeCalculationState(void) { return odeCalculationState; }
	double GetTimeInMeshPoint(void) { return t; }
	BzzVector GetYInMeshPoint(void) { return z[0]; }
	BzzVector GetY1InMeshPoint(void) { Division(z[1], hInNextStep, &va);return va; }
	void GetNumericalJacobian(BzzMatrix& J);
	void BzzPrintErrorState(void);
	void StepPrint(const char* sv, const char* pr = "Labels");
	void StepPrint(const char* sv, BzzVectorInt& isv, const char* pr = "Labels");
	void StepPrint(void (*stepPrintOut)(BzzVector& y, double t));
	void StepDebug(const char* sd, char deb = 1, double st = 0., double et = 0.);
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
	//ELIMINATOBzzVector Operator(double tF){return (*this)(tF);}
	BzzVector Operator(double tF);
	//{BzzVector result; result = (*this)(tF); return result;}

	BzzVector operator() (double tF, double tC);
	// ELIMINATOBzzVector Operator(double tF,double tC){return (*this)(tF,tC);}
	BzzVector Operator(double tF, double tC);
	//{BzzVector result; result = (*this)(tF, tC); return result;}
};

//	============================================================================
//	***********< BzzOdeNonStiffObject class: Adams Moulton algorithm >*************
//	============================================================================

class BzzOdeNonStiffObject : public BzzOdeMultiValueObject
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
	//		// ELIMINATO return Operator(tF);
	//		BzzVector result; result = Operator(tF); return result;
	//		}
	BzzVector operator () (double tF, double tC);
	//		{
	//		// ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
	{
		BzzOdeMultiValueObject::SetInitialConditions(y00, t00, ptrO);
		SetParameters();
	}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValueObject::SetInitialConditions(y00, t00);
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
	BzzOdeNonStiffObject(void)
		: BzzOdeMultiValueObject() {}

	BzzOdeNonStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
		: BzzOdeMultiValueObject(y00, t00, ptrO)
	{
		SetParameters();
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
	{
		BzzOdeMultiValueObject::SetInitialConditions(y00, t00, ptrO);
		SetParameters();
	}
};

class BzzOdeStiffBaseObject : public BzzOdeMultiValueObject
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
		JAC_CONST,	 // BzzOdeStiffObject(y0,t0,ptrSys,J);
		JAC_ANALYTIC, // BzzOdeStiffObject(y0,t0,ptrSys,ptrJac);
		JAC_NUMERICAL // BzzOdeStiffObject(y0,t0,ptrSys);
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
	void SetAnalyticalJacobian(void) { jacType = JAC_ANALYTIC; }

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzOdeStiffBaseObject(void)
		: BzzOdeMultiValueObject() {};

	BzzOdeStiffBaseObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
		: BzzOdeMultiValueObject(y00, t00, ptrO)
	{
		SetParameters();
	}
};

//	============================================================================
//	*********************< BzzOdeStiffObject for stiff problems >******************
//	============================================================================

class BzzOdeStiffObject : public BzzOdeStiffBaseObject
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
	//		// ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValueObject::SetInitialConditions(y00, t00);
		numStepFact = numStepJac = 0;
		iterConvergenceRate = 0;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
		//	SetParameters();
	}

	void SetInitialConditions(BzzVector& y00, BzzOdeSystemObject* ptrO);
	void SetInitialConditions(BzzVector& y00, double t00, BzzOdeSystemObject* ptrO,
		BzzMatrix* JJ);

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzOdeStiffObject(void)
		:BzzOdeStiffBaseObject() {};

	BzzOdeStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, ptrO);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, ptrO);
	}

	void GetLastJacobian(BzzMatrix* JJ) { *JJ = J; }
};

//	============================================================================
//	***********< BzzOdeSparseStiffObject for stiff and sparse problems >********
//	============================================================================

class BzzOdeSparseStiffObject : public BzzOdeStiffBaseObject
{
private:
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
	//	BzzMatrixBand B22;
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

	//	BzzMatrixSparseLockedByRows SLh;
	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked DaL;
	//	BzzMatrixSparseLockedByRows SL;

	void ProductFourBlocksDiagonalTridiagonal(void);
	void ProductFourBlocksDiagonalBand(void);
	void ProductFourBlocksDiagonalDiagonal(void);
	void ProductFourBlocksTridiagonalSparseSparseDense(void);
	void ProductFourBlocksBandSparseSparseDense(void);

	BzzVector ax1, ax2, bx1, bx2, cx1, cx2;

	//	BzzOdeSystemSparseObject *ptrObject;

	//	void (*sysJacSparse)(BzzVector &yy,double tt,BzzMatrixSparse &JJ);
	//	void (*sysJacBand)(BzzVector &yy,double tt,BzzMatrixBand &JJB);
	//	void (*sysJacTriDiaBlo)(BzzVector &yy,double tt,BzzMatrixTridiagonalBlocks &JJT);
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
	//		// ELIMINATO return Operator(tF,tC);
	//		BzzVector result; result = Operator(tF, tC); return result;
	//		}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzOdeMultiValueObject::SetInitialConditions(y00, t00);
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

	BzzOdeSparseStiffObject(void)
		:BzzOdeStiffBaseObject() {};

	// sparse generic
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, ptrO, JE);
	}
	// Jacobian numerical
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, ptrO, JE);
	}

	// Jacobian Block Matrices
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, ptrO, JBm);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, ptrO, JBm);
	}

	// block tridiagonal
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, ptrO, dim);
	}
	// Jacobian Block Tridiagonal
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dim);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, ptrO, dim);
	}

	// band
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, ptrO, low, up);
	}
	// Jacobian banded
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int low, int up);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, ptrO, low, up);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal tridiagonal
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrO, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrO, nV, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal diagonal
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrO, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, ptrO, n1, n2, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal band
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, ptrO, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, ptrO, nV, nb1, nb2, low, up);
	}

	// Four TridiagonalBlocks,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrO, dimBlock, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrO, dimBlock, J12, J21);
	}

	// Four BandGauss,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrO, lB, uB, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, ptrO, lB, uB, J12, J21);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear First and third Typologies
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		//		void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, ptrO, dia, QQ);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		//	   void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		//	   void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, ptrO, dia, QQ);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear Second and fourth Typologies
	BzzOdeSparseStiffObject(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, const char* fileM)
	{
		SetInitialConditions(y00, t00, ptrO, dia, QQ, fileM);
	}
	void SetInitialConditions(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, const char* fileM);
	void operator()(BzzVector& y00, double t00,
		BzzOdeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, const char* fileM)
	{
		SetInitialConditions(y00, t00, ptrO, dia, QQ, fileM);
	}

	void PrintSelectedMethod(void);
};

/*
//	============================================================================
//	********************< BzzOdeDiagonalBlocks >***************************
//	============================================================================

class BzzOdeDiagonalBlocks : public BzzOdeStiffBaseObject
	{
private:
	int	numMat,
			numVarBlock;
//	BzzVector d;
	BzzMatrix *E,EE;
	BzzFactorizedDiagonalBlocksAndSparse EF;
	BzzMatrixSparse *ES,SS;
	BzzVector en;
	void (*sysDiff2)(BzzVector &y,double t,BzzVector &f);

//	int lowerBand,upperBand;

//	void (*sysJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ);
//	virtual void ConvergenceRate(void){}
	virtual void ConvergenceRate(void);
	virtual void Jacobian(void);
	virtual void AnalyticalJacobian(void){};
	virtual void BuildBzzMatrixG(void);
	virtual void SolveLinearSystem(void);

public:
	BzzVector operator () (double tF)
		{
		return Operator(tF);
		}
	BzzVector operator () (double tF,double tC)
		{
		return Operator(tF,tC);
		}

	void SetInitialConditions(BzzVector &y00,double t00)
		{
		BzzOdeMultiValue::SetInitialConditions(y00,t00);
		}

	void SetInitialConditions(BzzVector &y00,double t00,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
		 void (*ptrSys2)(BzzVector &y,double t,BzzVector &f),
		 int numB,BzzMatrixSparse *SSS);
	void operator()(BzzVector &y00,double t00,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
		 void (*ptrSys2)(BzzVector &y,double t,BzzVector &f),
		 int numB,BzzMatrixSparse *SSS)
		{
		SetInitialConditions(y00,t00,ptrSys,ptrSys2,numB,SSS);
		}

//	============================================================================
//	******************************< constructors >******************************
//	============================================================================
	BzzOdeDiagonalBlocks(void)
		:BzzOdeStiffBaseObject(){};

	BzzOdeDiagonalBlocks(BzzVector &y00,double t00,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
		 void (*ptrSys2)(BzzVector &y,double t,BzzVector &f),
		 int numB,BzzMatrixSparse *SSS);
	};
*/

#endif // ODEMULTIOBJECT__HPP