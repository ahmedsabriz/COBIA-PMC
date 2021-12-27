// BZZMATH: Release 7.0

//	==============================< LINPROGRD.HPP >=============================
//	* BzzLinearProgramming: Class for linear programming							*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-2000	Date Written.
//	03-2010	Revision.

//	============================================================================
//	* This class solve the following linear programming								*
//	* Min F = sTx																					*
//	* Ex = e	(mE linear equations)															*
//	* Dx <= c (mD linear disequations)														*
//	* x = xE (equailitylower bounds)	iE constraints										*
//	* x >= xL (lower bounds)	iL constraints												*
//	* x <= xU (upper bounds)	iU constraint												*
//	****************************************************************************
//	============================================================================
//	******* Constructors																			*
// * E,D BzzMatrix or BzzMatrixSparse														*
//	* BzzLinearProgramming lp(&x0,&s,&E,&e,&D,&c,&iL,&xL,&iU,&xU)					*
//	* BzzLinearProgramming lp;																	*
//	* lp(&x0,&s,&iE,&xE,&iL,&xL,&iU,&xU,&E,&e,&D,&c)												*
//	* BzzLinearProgramming lp("FILELP.DAT");												*
//	****************************************************************************
//	******* Functions																				*
// * lp(n);																							*
//	****************************************************************************

#ifndef BZZ_LINEAR_PROGRAMMING_ATTIC
#define BZZ_LINEAR_PROGRAMMING_ATTIC

//	============================================================================
//	====================< class BzzLinearProgramming >========================
//	============================================================================
class BzzLinearProgramming : public BzzBaseClass
{
private:
	/*
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
		}methodUsed,methodToBeUsed;

	enum BzzLinearProgrammingMatrixJStatus
		{
		REBUILD_J,
		UPDATE_J
		}matrixJStatus;
*/
	enum bzzLinearProgrammingMatrixJType
	{
		FULL,
		SPARSE
	}matrixJType;
	/*
	//	enum bzzLinearProgrammingDegenerationStatus
	//		{
	//		NON_DEGENERATING,
	//		DEGENERATING
	//		}degenerationStatus;
	*/
	static const char* const BZZ_ERROR;
	static int count; // per whoAmI
	static int countInScope;
	static const int MAX_ITERATION;
	//	static const double DEG_FOR_TMAX;
	//	static const double DEG_FOR_TMAX_D;

	int whoAmI;
	FILE* fileLinearProgramming;
	FILE* fileBzzTmp;

	int	printTasks,
		printSubTasks;

	int	numVariables,
		numConstraints,
		numOriginalVariables,
		numU,
		numL,
		numE,
		mE, // num Equations
		mD, // num Disequations
		numRows,
		numColumns;

	BzzMatrix D, E;
	BzzMatrixSparse DS, ES;
	//////////////
	BzzMatrix A, U;
	//////////////
	BzzVector x0, s, xE, xL, xU, e, c;
	BzzVectorInt iE, iL, iU;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzLinearProgramming(void);

	// copy constructor
//	BzzLinearProgramming(const BzzLinearProgramming &rval);

	// constructor
	BzzLinearProgramming(BzzVector* x00, BzzVector* ss,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzMatrix* DD, BzzVector* cc);

	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzMatrix* DD, BzzVector* cc);

	BzzLinearProgramming(BzzVector* x00, BzzVector* ss,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* cc);

	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* cc);

	// constructor FILE .BZZ
	BzzLinearProgramming(char* fileLP);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzLinearProgramming(void);

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
};

#endif // BZZ_LINEAR_PROGRAMMING_DOUBLE