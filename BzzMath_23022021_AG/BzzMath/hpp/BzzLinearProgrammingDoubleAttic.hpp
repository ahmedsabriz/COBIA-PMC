// BZZMATH: Release 7.0

//	==============================< LINPROGRD.HPP >=============================
//	* BzzLinearProgrammingAttic: Class for linear programming							*
//	* Description: Numerical Analysis and Software in C++ Vol. 2 (Chapter )		*
//	*					by G. Buzzi-Ferraris														*
// * Examples: c:\bzzmath\exampled\dxlinprog.cpp						 				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-2000	Date Written.

//	============================================================================
//	* This class solve the following linear programming								*
//	* Min F = sTx																					*
//	* Ex = e	(mE linear equations)															*
//	* Dx >= d (mD linear disequations)														*
//	* x >= L (lower bounds)																		*
//	* x <= U (upper bounds)																		*
//	****************************************************************************
//	============================================================================
//	******* Constructor																			*
//	* BzzLinearProgrammingAttic lp(&x0,&s,&E,&e,&D,&d,&iL,&xL,&iU,&xU)					*
//	* BzzLinearProgrammingAttic lp;																	*
//	* lp(&x0,&s,&E,&e,&D,&d,&iL,&xL,&iU,&xU)												*
//	* BzzLinearProgrammingAttic lp("FILELP.DAT");												*
//	****************************************************************************
//	******* Functions																				*
// * lp(n);																							*
//	****************************************************************************

#ifndef BZZ_LINEAR_PROGRAMMING_DOUBLE
#define BZZ_LINEAR_PROGRAMMING_DOUBLE

//	============================================================================
//	====================< class BzzLinearProgrammingAttic >========================
//	============================================================================
class BzzLinearProgrammingAttic : public BzzBaseClass
{
private:
	enum LinearProgrammingCalculationState
	{
		PROBLEM_SOLVED_STATE = 0,
		INITIALIZATION_STATE = 1,
		CONTINUATION_STATE = 2,
		EXCESSIVE_WORK_STATE = -1,
		NON_FEASIBLE_STATE = -2,
		SINGULAR_MATRIX = -3,
		INSOLUBLE_DEGENERATION = -4
	}linearProgrammingCalculationState;

	enum BzzLinearProgrammingMethodUsed
	{
		UNKNOWN = -1,
		START = 0,
		DEALLOCATION_ONE = 1,
		DEALLOCATION_TWO = 2,
		DEALLOCATION_THREE = 3,
		DETACHMENT = 4,
		ALLOCATION_ONE = 5,
		ALLOCATION_TWO = 6,
		ALLOCATION_THREE = 7,
		//		FEASIBLE_ONE = 8,
		//		FEASIBLE_TWO = 9,
		DEGENERATION = -2,
		TERMINATION = 10
	}methodUsed, methodToBeUsed;

	enum BzzLinearProgrammingMatrixJStatus
	{
		REBUILD_J,
		UPDATE_J
	}matrixJStatus;

	enum bzzLinearProgrammingMatrixJType
	{
		J_BULK_NO,
		J_BULK_ONE,
		J_BULK_DENSE,
		J_BULK_SPARSE,
		J_BULK_NO_BORDER,
		J_BULK_ONE_BORDER,
		J_BULK_DENSE_BORDER,
		J_BULK_SPARSE_BORDER
	}matrixJType;

	//	enum bzzLinearProgrammingDegenerationStatus
	//		{
	//		NON_DEGENERATING,
	//		DEGENERATING
	//		}degenerationStatus;

	static const char* const BZZ_ERROR;
	static int count; // per whoAmI
	static int countInScope;
	static const int MAX_ITERATION;
	static const double DEG_FOR_TMAX;
	static const double DEG_FOR_TMAX_D;

	int whoAmI;
	FILE* fileLinearProgramming;
	FILE* fileBzzTmp;

	int	printTasks,
		printSubTasks;

	char bDefault;
	int	numVariables,
		originalNumVariables,
		original1NumVariables,
		mE, // num Equations
		mD; // num Disequations
//			mL, // num lower
//			mU; // num upper

	double startTime, endTime, initialStartTime;

	int	nD, // num Inequality active
		nL, // num lower active
		nU, // num upper active
		nDOut, // num Violated
		nLOut, // num Violated
		nUOut, // num Violated
		nEOut,
		nDDeg, // num Violated
		nLDeg, // num Violated
		nUDeg, // num Violated
//			nX, // nL + nU
nSelected; // numVariables - nX

//	int	nDisMark1Pos,
//			nInfMark1Pos,
//			nSupMark1Neg,
//			nDisMark1Neg,
//			nInfMark1Neg,
//			nSupMark1Pos,
//			nDisMark2Pos,
//			nInfMark2Pos,
//			nSupMark2Neg,
//			nDisMark2Neg,
//			nInfMark2Neg,
//			nSupMark2Pos,
//			nArt,
//			nPosMark1,
//			nNonArt,
//			nNonZero;

//	int	numAllocationOne,
//			numAllocationTwo,
//			numAllocationThree,
//			numDeallocationOne,
//			numDeallocationTwo,
//			numDeallocationThree,
//			numDetachment;

//	int	numB,numN; // variables in Jb and Jn
//	int	numRowsW,numColumnsW;
	int	singularAndIncompatibleE, iterDegeneration;

	//	int forceTwo;

	//	int indexFindx,typeFindx;  //type 0 xn1, 1 L, 2 U,3 D

	double	epsA,
		epsR,
		epsADeg,
		epsRDeg,
		//			xn1,
		//			sXn1,
		//			gXn1,
		//			pXn1,
		//			yXn1,
		//			xn1Ott,
		//			xn1Min,
		//			lambdaXn1,
		//			lXn1,
		F,
		fConst, // F = fConst + sixi
		fOtt,
		fFeasible;

	// lMark = -3 vincolo di uguaglianza!!!!!!!
	BzzVectorInt	lMark, // 0 initially, + 1 every time, dimensioned numVariables
		uMark, // 0 initially, + 1 every time, dimensioned numVariables
		dMark, // 0 initially, + 1 every time, dimensioned mD
		ivmark,
		//						lMarkDegener, // 0 initially, 1 if modified, dimensioned numVariables
		//						uMarkDegener, // 0 initially, 1 if modified, dimensioned numVariables
		//						dMarkDegener, // 0 initially, 1 if modified, dimensioned mD
		ilConstrained, // 1 active 0 passive, dimensioned numVariables
		iuConstrained, // 1 active 0 passive, dimensioned numVariables
		idConstrained, // 1 active 0 passive dimensioned mD
		ilConstrainedOut, // 1 active 0 passive, dimensioned numVariables
		iuConstrainedOut, // 1 active 0 passive, dimensioned numVariables
		idConstrainedOut, // 1 active 0 passive dimensioned mD
		ilConstrainedOld, // 1 active 0 passive, dimensioned numVariables
		iuConstrainedOld, // 1 active 0 passive, dimensioned numVariables
		idConstrainedOld, // 1 active 0 passive dimensioned mD
		ieConstrainedOut, // 1 active 0 passive, dimensioned numVariables
		ilDegener, // 1 active 0 passive, dimensioned numVariables
		iuDegener, // 1 active 0 passive, dimensioned numVariables
		idDegener, // 1 active 0 passive dimensioned mD
		ieDegener, // 1 active 0 passive dimensioned mE
		dRows, // D rows active, dimensioned nD
//						dRowsToAdd,
//						lConstraintToAdd,
//						uConstraintToAdd,
originalRowsE,
originalRowsD,
unboundedVariables,
ixFixed,
ixUsed,
ixFixed1,
ixUsed1,

//						idModified, // D rows modified dimensioned nDModified
//						ilModified, // inf modified dimensioned nLModified
//						iuModified, // sup modified dimensioned nUModified
//						iX, // variables on constraint dimensioned nX
//						iXl, // variables on inferior constraint dimensioned nL
//						iXlMark, // mark on iXl
//						iXu, // variables on superior constraint dimensioned nU
//						iXuMark, // mark on iXu
ordRows,
ordColumns,
selectedColumns,
//fixVariables, //??
eqType, // dimensioned numVariables
eqRefer, // dimensioned numVariables
lNonZero,
lNonArt,
lArt,
lPosMark1,
iB, iN, // variables in Jb and Jn
eNumEquationsForEachVariable,
eNumVariablesForEachEquation,
dNumEquationsForEachVariable,
dNumVariablesForEachEquation,
wNumEquationsForEachVariable,
wNumVariablesForEachEquation,
linearCombinations,
variableToBeUnconstrained;

	BzzVector	x0,
		x,
		xOtt,
		xFeasible,
		xFixed,
		xFixed1,
		s,
		g,
		e,
		d,
		L,
		U,
		dOriginal,
		LOriginal,
		UOriginal,
		ek, // dimensioned mE
		dk, // dimensioned mD
		lk, // dimensioned mL
		uk, // dimensioned mL
		epsL,
		epsU,
		//					ekXn1, // coeff for xn1 dimensioned mE
		//					dkXn1, // coeff for xn1 dimensioned mD
		//					lkXn1, // coeff for xn1 dimensioned mL
		//					ukXn1, // coeff for xn1 dimensioned mL
		//					dModified, // original constraint dimensioned nDModified
		//					lModified, //  original constraint dimensioned nLModified
		//					uModified, //  original constraint dimensioned nUModified
		unorm2E,
		norm2EX,
		norm2D,
		unorm2D,
		norm2Dl,
		norm2DX,
		unorm2Variables,
		unorm2VariablesInE,
		//					feasibleWeightD,
		//					feasibleWeightL,
		//					feasibleWeightU,
		p,//pB,pN,
		pFeasible,
		l, lB, lN,
		lambda,//lambdaB,lambdaN,
		dDenForTmax, // for tMax dimensioned mD
		dNumForTmax, // for tMax dimensioned mD
		lDenForTmax, // for tMax dimensioned mL
		lNumForTmax, // for tMax dimensioned mL
		uDenForTmax, // for tMax dimensioned mU
		uNumForTmax, // for tMax dimensioned mU
		lBase, lOtt,
		y,
		sB, sN,
		h;

	BzzVector at;
	BzzVector bt;
	BzzVectorInt it;
	BzzVectorInt jt;
	int nt;

	BzzMatrixSparseLockedByRows E, D, W, Jn;
	BzzMatrixSparse Jb; // divenera' undefined
	BzzMatrix Jbm;
	BzzMatrixSparseStore store;
	BzzFactorizedGauss gDense;

	///////////////////////////////////////////////////////////////////////////////
		// for Feasible and Degenerate
	char isActiveZ;
	int rowsTotZ, columnsTotZ;
	BzzMatrixSparseLockedByRows Z;
	BzzFactorizedSparseLockedLQ B;
	//	BzzFactorizedLQ B;
	int FindFeasiblePoint(double precision = 1.e-7, int maxIter = 5000,
		int maxBreak = 30);
	char InitializeBLocked(void);
	void BzzPrintErrorState(void);
	///////////////////////////////////////////////////////////////////////////////
	void Initialize(BzzVector* x00, BzzVector* ss,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* dd,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU);
	//	void Copy(const BzzMatrixSparse &rval);
	void SetInsideBounds(BzzVector* xxx);
	int StartingPoint(void);
	int GlobalStartingDegeneration(void);
	int SemiGlobalStartingDegeneration(void);
	int StartingDegeneration(void);
	void GetEqualityDegenerations(BzzVectorInt* vNonInE);
	void GetInequalityDegenerations(BzzVectorInt* vNonInD);
	void GetInequalityDegenerations(void);
	//	void StartingWithE(void);
	//	void StartingWithoutE(void);
	void Residuals(void);
	void ContraintsCount(void);
	void ContraintsDegCount(void);
	//	void ContraintsCountXn1(void);
	//	void ContraintsCountNew(void);
	//	void ConstraintCheck(void);
	//	void CheckFeasibility(void);
	//	void CheckFeasibilityNonDegenerate(void);
	//	void CheckFeasibilityDegenerate(void);
	//	void BuildMatrixW(void);
	//	void FindX(void);
	//	void BuildOrUpdateJ(void);
	//	void BuildMatrixJ(void);
	void BuildStartingMatrixJ(void);
	//	void UpdateMatrixJ(void);
	//	void LambdaEvaluation(void);
	//	void GetVectorP(void);
	//	void SelectMethod(void);
	//	void AllocationOne(void);
	//	void AllocationTwo(void);
	//	void DeallocationOne(void);
	//	void DeallocationTwo(void);
	//	void DeallocationThree(void);
	//	void Detachment(void);
	//	void FeasibleOne(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzLinearProgrammingAttic(void);

	// copy constructor
//	BzzLinearProgrammingAttic(const BzzLinearProgrammingAttic &rval);

	// constructor
	BzzLinearProgrammingAttic(BzzVector* x00, BzzVector* ss,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* dd,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU);

	// constructor FILE .BZZ
	BzzLinearProgrammingAttic(char* fileLP);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzLinearProgrammingAttic(void) { countInScope--; };

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
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
	void GetSolution(BzzVector* xo) { *xo = xOtt; }

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	void operator () (BzzVector* x00, BzzVector* ss,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* dd,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU);
	int operator () (int nIt = 0);
};

#endif // BZZ_LINEAR_PROGRAMMING_DOUBLE