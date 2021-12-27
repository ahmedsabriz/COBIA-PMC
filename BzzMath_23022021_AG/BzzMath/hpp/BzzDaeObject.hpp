// BZZMATH: Release 7.0

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-2006
//	05-2008	Added GetTimeInMeshPoint function.
//	05-2008	Added GetYInMeshPoint function.
//	05-2008	Added GetY1InMeshPoint function.

////////////////// Release 6.0
// 07-2009	Added OpenMP parallelization.
// 02-2010	Added StepDebug and StopStepDebug functions.
//	05-2010	Bug fixed in Linear System Solution.

#ifndef DAEMULTI_OBJECT_HPP
#define DAEMULTI_OBJECT_HPP

//	============================================================================
//	=======================< class BzzDaeSystemObject >===================
//	============================================================================
class BzzDaeSystemObject : public BzzBaseClass
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
//	=================< class BzzDaeObjectJacobian >========================
//	============================================================================
class BzzDaeObjectJacobian
{
private:
	double ZERO_DER, ETA2, BETA;
	double hJ, hJf, yh, hInv, hf;
	BzzVector aaa, bbb;
public:
	BzzDaeSystemObject* ptrObject;
	void GetJacobianNoConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
	void GetJacobianConstraint(int jStart, int jEnd, double t, double hff,
		BzzVector& vb, BzzVector& yMax, BzzVector& errorWeightBzzVector,
		BzzVector& f, BzzMatrix& J);
};

class BzzDaeMultiValueObject : public BzzBaseClass
{
	friend void bzzNlsDaeObject(BzzVector& x, BzzVector& g);
	friend class BzzMyNonLinearSystemObjectDaeObject;
	friend class BzzMyNonLinearSystemObjectSparseDaeObject;
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

	enum BzzDaeSparseJacobyanType
	{
		FULL,
		SPARSE,
		SPARSE_BAND, // ***
		SPARSE_FULL, // ***
		BAND,
		BLOCK,
		TRIDIAGONAL_BLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_DIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_TRIDIAGONALBLOCK,
		FOUR_A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_BAND,
		FOUR_A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,
		FOUR_A11_FACTORIZED_BAND_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_FIRST, // ***
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_SECOND, // ***
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_THIRD, // ***
		DIAGONAL_BLOCK_AND_LOCKED_LINEAR_FOURTH // ***
	}daeJacobianType;

	enum BzzDaeCalculationState
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
		numThread;

	BzzDaeObjectJacobian* getJac;
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
		printSubTasks,
		stopIntegrationBeforeRecalcuatingJacobian,
		stopIntegrationWhenSumAbsY1IsLessThan;

	BzzVectorInt varSaved;

	double tInMinH;
	BzzVector yInMinH;

	BzzVector	y0,
		y,
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

	//	void (*sysDiff)(BzzVector &y,double t,BzzVector &f);
	BzzDaeSystemObject* ptrObject;
	void (*stepPrintOut)(BzzVector& y, double t);

	void Initialize(BzzVector& y00, double t00,
		BzzDaeSystemObject* ptrO);
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
		BzzDaeSystemObject* ptrO);
	void SetInitialConditions(BzzVector& y00, double t00);

public:
	double tt;
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default
	BzzDaeMultiValueObject(void);

	BzzDaeMultiValueObject(BzzVector& y00, double t00,
		//		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f));
		BzzDaeSystemObject* ptrO);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	void Deinitialize(void);
	~BzzDaeMultiValueObject(void);

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
	void StopIntegrationBeforeRecalcuatingJacobian(int nJ = 1)
	{
		stopIntegrationBeforeRecalcuatingJacobian = nJ;
	}
	void StopIntegrationWhenSumAbsY1IsLessThan(double y1v)
	{
		stopIntegrationWhenSumAbsY1IsLessThan = 1;sumAbsY1 = y1v;
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

class BzzDaeBaseObject : public BzzDaeMultiValueObject
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
		JAC_CONST,	 // BzzDaeObject(y0,t0,ptrSys,J);
		JAC_ANALYTIC, // BzzDaeObject(y0,t0,ptrSys,ptrJac);
		JAC_NUMERICAL // BzzDaeObject(y0,t0,ptrSys);
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
	BzzDaeBaseObject(void)
		: BzzDaeMultiValueObject() {};

	BzzDaeBaseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO);

	BzzDaeBaseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO);
};

//	============================================================================
//	*********************< BzzDaeObject for stiff problems >******************
//	============================================================================

class BzzDaeObject : public BzzDaeBaseObject
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
		BzzDaeMultiValueObject::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
	}

	//	void SetInitialConditions(BzzVector &y00,	BzzDaeSystemObject *ptrO);
	//	void SetInitialConditions(BzzVector &y00, double t00,BzzDaeSystemObject *ptrO,
	//		BzzMatrix *JJ);

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzDaeObject(void)
		:BzzDaeBaseObject() {};

	BzzDaeObject(BzzVector& y00, double t00, BzzMatrix& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatrix& JJ);

	BzzDaeObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, AA, ptrO);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, AA, ptrO);
	}

	BzzDaeObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, iDer, ptrO);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO)
	{
		SetInitialConditions(y00, t00, iDer, ptrO);
	}

	void GetLastJacobian(BzzMatrix* JJ) { *JJ = J; }
};

//	============================================================================
//	***********< BzzDaeSparseObject for stiff and sparse problems >********
//	============================================================================

class BzzDaeSparseObject : public BzzDaeBaseObject
{
private:
	BzzMatrixCoefficientsExistence Je;
	BzzFactorizedGauss G;
	BzzFactorizedBandGauss B;
	BzzFactorizedMatrixBlocksGauss Fm;
	BzzFactorizedTridiagonalBlocksGauss Tf;
	BzzMatricesExistence Jm; // ***
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
	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked DaL; // ***

	void ProductFourBlocksDiagonalTridiagonal(void);
	void ProductFourBlocksDiagonalBand(void);
	void ProductFourBlocksDiagonalDiagonal(void);
	void ProductFourBlocksTridiagonalSparseSparseDense(void);
	void ProductFourBlocksBandSparseSparseDense(void);
	BzzVector ax1, ax2, bx1, bx2, cx1, cx2;

	//	void (*sysJac)(BzzVector &yy,double tt,BzzMatrixSparse &JJ);
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
	virtual int SolveNonLinearSystem(void);
	virtual int SolveStabilizeNonLinearSystem(void)
	{
		BzzError("Stabilize Non Linear System sparse DAE");
		return 1;
	}

public:
	BzzVector operator () (double tF);
	//{
	//// ELIMINATO return Operator(tF);
	//BzzVector result; result = Operator(tF); return result;
	//}
	BzzVector operator () (double tF, double tC);
	//{
	//// ELIMINATO return Operator(tF,tC);
	//BzzVector result; result = Operator(tF, tC); return result;
	//}

	void SetInitialConditions(BzzVector& y00, double t00)
	{
		BzzDaeMultiValueObject::SetInitialConditions(y00, t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for (int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
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

		// ??????
	void operator()(BzzVector& y00, double t00)
	{
		SetInitialConditions(y00, t00);
	}
	BzzDaeSparseObject(void)
		:BzzDaeBaseObject() {};

	// Jacobian numerical
		// sparse generic
	///////////////////////////////////////////////////////////////////////////////
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, AA, ptrO, JE);
	}
	// TODO
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, AA, ptrO, JE);
	}
	///////////////////////////////////////////////////////////////////////////////
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, JE);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatrixCoefficientsExistence* JE)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, JE);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Jacobian Block Matrices
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, AA, ptrO, JBm);
	}
	// TODO
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, AA, ptrO, JBm);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, JBm);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzMatricesExistence* JBm)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, JBm);
	}

	///////////////////////////////////////////////////////////////////////////////
	// block tridiagonal
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, AA, ptrO, dim);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int dim);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, AA, ptrO, dim);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dim);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dim);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dim)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dim);
	}

	///////////////////////////////////////////////////////////////////////////////
	// band
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrO, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int low, int up);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrO, low, up);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int low, int up);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, low, up);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal tridiagonal
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrO, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrO, nV, nb1, nb2);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, nV, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, nV, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal diagonal
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrO, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, AA, ptrO, n1, n2, nb1, nb2);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, n1, n2, nb1, nb2);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int n1, int n2, int nb1, int nb2)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, n1, n2, nb1, nb2);
	}

	///////////////////////////////////////////////////////////////////////////////
	// four blocks diagonal diagonal diagonal band
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrO, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00, BzzMatrixSparse& AA,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, AA, ptrO, nV, nb1, nb2, low, up);
	}

	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, nV, nb1, nb2, low, up);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int nV, int nb1, int nb2, int low, int up)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, nV, nb1, nb2, low, up);
	}

	// Four TridiagonalBlocksGauss,sparse,sparse,dense
	// (....,int dimBlock,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dimBlock, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int dimBlock, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dimBlock, J12, J21);
	}

	// Four BandGauss,sparse,sparse,dense
	// (....,int lB,int uB,BzzMatrixCoefficientsExistence &J12,BzzMatrixCoefficientsExistence &J21);
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, lB, uB, J12, J21);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		int lB, int uB, BzzMatrixCoefficientsExistence& J12, BzzMatrixCoefficientsExistence& J21)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, lB, uB, J12, J21);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear First and third Typologies
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		//		void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dia, QQ);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		//	   void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		//	   void (*sysDiffD)(BzzVector &y,double t,BzzVector &f),
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dia, QQ);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Diagonal Block and Sparse Locked Linear Second and fourth Typologies
	BzzDaeSparseObject(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dia, QQ, fileM);
	}
	void SetInitialConditions(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM);
	void operator()(BzzVector& y00, double t00, BzzVectorInt& iDer,
		BzzDaeSystemObject* ptrO,
		BzzVectorInt dia, BzzMatrixSparseLockedByRows* QQ, char* fileM)
	{
		SetInitialConditions(y00, t00, iDer, ptrO, dia, QQ, fileM);
	}

	void PrintSelectedMethod(void);
};

/*
//	============================================================================
//	********************< BzzDaeDiagonalBlocksObject >***************************
//	============================================================================
class BzzDaeDiagonalBlocksObject : public BzzDaeBaseObject
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
	void SetInitialConditions(BzzVector &y00,double t00)
		{
		BzzDaeMultiValueObject::SetInitialConditions(y00,t00);
		maxOrder = MAX_ORDER;
		maxOrderUsed = minOrder = 1;
		parOrderPM1 = PAR_ORDER_PM1;
		parOrderPP1 = PAR_ORDER_PP1;
		for(int i = 0;i <= MAX_ORDER;i++)
			parOrderP[i] = PAR_ORDER_P;
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
	BzzDaeDiagonalBlocksObject(void)
		:BzzDaeBaseObject(){};

	BzzDaeDiagonalBlocksObject(BzzVector &y00,double t00,BzzMatrixSparse &AA,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
		 void (*ptrSys2)(BzzVector &y,double t,BzzVector &f),
		 int numB,BzzMatrixSparse *SSS);
	BzzDaeDiagonalBlocksObject(BzzVector &y00,double t00,BzzVectorInt &iDer,
		 void (*ptrSys)(BzzVector &y,double t,BzzVector &f),
		 void (*ptrSys2)(BzzVector &y,double t,BzzVector &f),
		 int numB,BzzMatrixSparse *SSS);
	};
*/
//=============================================================================================
class BzzMyNonLinearSystemObjectDaeObject : public BzzMyNonLinearSystemObject
{
private:
	BzzDaeObject* ptrNlsDae;
	BzzVector vb;
public:
	BzzMyNonLinearSystemObjectDaeObject(void) {};
	BzzMyNonLinearSystemObjectDaeObject(BzzDaeObject* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	void operator()(BzzDaeObject* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

//=============================================================================================
class BzzMyNonLinearSystemObjectSparseDaeObject : public BzzMyNonLinearSystemSparseObject
{
private:
	BzzDaeSparseObject* ptrNlsDae;
	BzzVector vb;
public:
	BzzMyNonLinearSystemObjectSparseDaeObject(void) {}
	BzzMyNonLinearSystemObjectSparseDaeObject(BzzDaeSparseObject* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	void operator()(BzzDaeSparseObject* ptrNlD)
	{
		ptrNlsDae = ptrNlD;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

#endif // DAEMULTI_OBJECT_HPP