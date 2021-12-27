// TODO: JACOBIAN_BLOCK

// BZZMATH: Release 7.0

//	=================< BzzNonLinearSystemSparseObject.hpp >=====================
//	* BzzNonLinearSystemSparseObject class													*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 25, 26)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\BzzAdavanced\										*
// *					NonLinearSystems\NonLinearSystemSparse.cpp				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1997	Date Written
//	06-1998	Bug fixed in numFunctionsForGradient initialization.
//	09-1998	Introduction of new functions Restart.
//	09-1998	Modified JACOBIAN_BAND implementation
//	09-1998	Modified JACOBIAN_TRIDIAGONAL_BLOCK implementation
//	05-2001	SetWeights.
//	07-2001	Added SetTasksPrint, SetSubTasksPrint.

////////////////// Release 6.0
// 12-2011	Added GetCalculationState function.
// 12-2011	Added GetMaximumResidual function.

#ifndef NON_LINEAR_SYSTEM_DOUBLE_SPARSE_OBJECT_HPP
#define NON_LINEAR_SYSTEM_DOUBLE_SPARSE_OBJECT_HPP

//	============================================================================
//	=======================< class BzzNonLinearSystemObject >===================
//	============================================================================
class BzzMyNonLinearSystemSparseObject : public BzzBaseClass
{
public:
	virtual void GetResiduals(BzzVector& x, BzzVector& f) {};
	virtual void ObjectBzzPrint(void) {};
};

//	============================================================================
//	===================< class BzzNonLinearSystemSparseObject >=================
//	============================================================================
//class BzzCSTRNetwork;
class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked;
class BzzNonLinearSystemSparseObject : public BzzBaseClass
{
	friend class BzzBVP;
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
		JACOBIAN_FULL,
		JACOBIAN_BAND,
		JACOBIAN_SPARSE,
		JACOBIAN_BLOCK,
		JACOBIAN_TRIDIAGONAL_BLOCK,
		JACOBIAN_FOR_CSTR_NETWORK,
		JACOBIAN_LOCKED_GENERIC,
		JACOBIAN_LOCKED_BAND
	}nlsJacobianType;

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

	enum BzzNonLinearSystemSparseObjectWeightsType
	{
		AUTOMATIC_WEIGHTS,
		FIXED_WEIGHTS,
		UNITARY_WEIGHTS
	}nlsWeightsType;

	static const char* const BZZ_ERROR;

	static const double	DEFAULT_ABSOLUTE_TOLERANCE,
		DEFAULT_RELATIVE_TOLERANCE,
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

	int	numVariables, numVariablesNew, lowNew, upNew,
		maxNumFunctionsInTheCall,
		maxNumFunctionsTotal,
		numFunctionsInTheCall,
		numFunctionsTotal,
		numQuasiNewtons,
		numNewtons,
		numGradients,
		numFunctionsForGradient,/////////////
		numMonodimensional,/////////////////
		numNumericalJacobians,
		numFactorizations,
		numLinearSolutions,
		numTrials,
		numGroups,
		numDiagonalMatrices,
		lowerBand, upperBand,
		blockDimensions, // = 0 if different
		printTasks,
		printSubTasks,
		maxNewtonCallsNumber;

	char	//control,
		jacobianUpatating,
		singular,
		stepReduced,
		settingTollerance,
		stop,
		stopRestart,
		noJacobian,
		contraintsViolated,
		first;

	double	phiW,
		phiNew,
		phi1W,
		phi1New,
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
		u, v, z, fz, fu,
		timeSystem,
		timeJacobian,
		timeShubert,
		timeLinearSystem;

	BzzVector	x0v, x0, f0, f0v,
		xi, fi,
		xi1, fi1,
		xMin, xMax, xDimensions,
		pi,
		pi1,
		dxi, dfi,
		gi,
		weights,
		errorWeightBzzVector,
		tollAbs,
		tollRel,
		aux, aux1;

	BzzVectorInt* eqVar,
		numVar;

	//	BzzCSTRNetwork *cstr;
	BzzMatrix	Xi,
		Xi1,
		F0,
		Fi1,
		Pi,
		Pi1;
	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked Pcstr;

	BzzMatrixCoefficientsExistence Je;
	BzzMatrixSparse S;
	BzzMatrixSparse SS;
	BzzFactorizedSparseUnspecified U;

	BzzMatricesExistence Jm;
	BzzMatrixBlocks B;
	BzzMatrixBlocks BB;
	BzzFactorizedMatrixBlocksGauss F;

	BzzMatrixTridiagonalBlocks T;
	BzzMatrixTridiagonalBlocks TT;
	BzzFactorizedTridiagonalBlocksGauss JTrid;

	BzzMatrix Ji;
	BzzFactorizedGauss JiGauss;
	BzzMatrixBand JB;
	BzzFactorizedBandGauss JBand;

	//	BzzFactorizedLQ JiLQ;

	//	void (*ptrSystem)(BzzVector &xx,BzzVector &ff);
	BzzMyNonLinearSystemSparseObject* ptrObject;
	//	void (*ptrJacobian)(BzzVector &xx, BzzMatrix &JJ);

	void Initialize(void);
	void ReInitialize(void);
	void BzzPrintErrorState(void);
	void ControlObjectiveFunctions(void);
	char ControlPrevision(double gamma);
	double ControlNewton(void);
	void Controlxi(void);
	void Controlxi1(void);
	void Controlxi1Newton(void);
	void PhiWInxi1(void);
	void PhiNewInxi1(void);
	void Weights(void);
	void ErrorWeightBzzVector(void);
	void GetJacobian(void);
	void SparseJacobian(void);
	void BlockJacobian(void);
	void TridiagonalBlocksJacobian(void);
	void NewtonPrevision(void);
	void QuasiNewtonPrevision(void);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================

		// default constructor
	BzzNonLinearSystemSparseObject(void);

	// copy constructor
	BzzNonLinearSystemSparseObject
	(const BzzNonLinearSystemSparseObject& rval);

	// Sparse Jacobian
	BzzNonLinearSystemSparseObject
	(BzzVector& x0,
		//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrixCoefficientsExistence* S);
	void operator()
		(BzzVector& x0,
			//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
			BzzMyNonLinearSystemSparseObject* ptrO,
			BzzMatrixCoefficientsExistence* S);

	// Block Jacobian
	BzzNonLinearSystemSparseObject
	(BzzVector& x0,
		//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatricesExistence* B);

	// Tridiagonal Jacobian
	BzzNonLinearSystemSparseObject
	(BzzVector& x0,
		//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		int dim);
	void operator()
		(BzzVector& x0,
			//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
			BzzMyNonLinearSystemSparseObject* ptrO,
			int dim);

	// Banded Jacobian
	BzzNonLinearSystemSparseObject
	(BzzVector& x0,
		//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		int low, int up);
	void operator()
		(BzzVector& x0,
			//		 void (*ptrS)(BzzVector &xx,BzzVector &ff),
			BzzMyNonLinearSystemSparseObject* ptrO,
			int low, int up);

	// CSTR_NETWORK
//	BzzNonLinearSystemSparseObject
//		 (BzzMatrix *X0,
//		 BzzCSTRNetwork *P);
//	void operator()
//		 (BzzMatrix *X0,
//		 BzzCSTRNetwork *P);

//	============================================================================
//	*********************************< destructor >*****************************
//	============================================================================
	~BzzNonLinearSystemSparseObject(void);

	//	============================================================================
	//	******************************< Functions >*********************************
	//	============================================================================
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
	int NumNumericalJacobians(void) { return numNumericalJacobians; }
	int NumFactorizations(void) { return numFactorizations; }
	int NumLinearSolutions(void) { return numLinearSolutions; }
	int NumNewtons(void) { return numNewtons; }
	int NumQuasiNewtons(void) { return numQuasiNewtons; }
	int NumGradients(void) { return numGradients; }
	int NumFunctionsForGradient(void) { return numFunctionsForGradient; }
	int NumFunctionsForNumericalJacobian(void) { return numNumericalJacobians * numGroups; }
	int NumMonodimensional(void) { return numMonodimensional; }
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
	void SetWeights(BzzVector& w);
	void Restart(BzzVector& x00);
	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff));
		BzzMyNonLinearSystemSparseObject* ptrO);

	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrix& jac, char stopR);
	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrixSparse& jac, char stopR);
	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrixBand& jac, char stopR);
	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrixBlocks& jac, char stopR);
	void Restart(BzzVector& x00,
		//		void (*ptrS)(BzzVector &xx,BzzVector &ff),
		BzzMyNonLinearSystemSparseObject* ptrO,
		BzzMatrixBlocks& jac, char stopR, int dim);

	void GetSolution(BzzVector* xx, BzzVector* ff,
		double* pW = 0, double* pN = 0);
	void GetWeights(BzzVector* w)
	{
		*w = weights;
	}

	char operator ()(void) { return (*this)(BzzNonLinearSystemSparseObject::DEFAULT_MAX_FUNCTIONS_IN_THE_CALL); }
	char operator ()(int ns);
	char HyperSolve(int ni);
	void operator ()(BzzVector& x00);
	void SetMaxNewtonCallsNumber(int s)
	{
		maxNewtonCallsNumber = s;
	}
	void SetStopForBadConvergence(char s)
	{
		if (s > 0)stopForBadConvergence = s;
	}
};

#endif // NON_LINEAR_SYSTEM_DOUBLE_SPARSE_OBJECT_HPP