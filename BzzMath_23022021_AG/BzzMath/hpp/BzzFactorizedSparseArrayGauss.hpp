// BZZMATH: Release 7.0

// ====================< BzzFactorizedSparseArrayGauss.HPP >===============
// * Class BzzFactorizedSparseArrayGauss for sparse non structured linear	*
// * systems solution																			*
// * Examples: c:\bzzmath\examples\.cpp*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2006	Date Written.
//	11-2010	Added BzzMatrixSparseLocked class.

// ============================================================================
// ****** Constructors for BzzFactorizedSparseArrayGauss:						*
// * BzzFactorizedSparseArrayGauss s; // default									*
// * BzzFactorizedSparseArrayGauss s("SparseMatrix.dat");	 					*
// ****************************************************************************
// ****************************************************************************
// ***** Implemented operations :															*
// ** Solve(&s,&bx);																				*
// ****************************************************************************

#ifndef BZZ_BZZ_FACTORIZED_DOUBLE_SPARSE_ARRAY_GAUSS_HPP
#define BZZ_BZZ_FACTORIZED_DOUBLE_SPARSE_ARRAY_GAUSS_HPP

// ============================================================================
// ===============< class BzzFactorizedSparseArrayGauss >==================
// ============================================================================
class BzzMatrixSparseLockedByRows;
class BzzFactorizedSparseArrayGauss : public BzzBaseClass
{
private:
	enum FactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;
	static const char* const BZZ_ERROR;
	int numRows, numColumns, numElements;
	int	ordering, // default 0 (without ordering)
		iConditionError,
		iNumberOfCorrectDigits,
		iVariableCondition,
		iFile,
		iSolutionNumber,
		iSolutionNumberNumberOfCorrectDigits,
		iVariableSolutionNumber,
		numTriangle,
		singular,
		minColumns, minRows, maxColumns, maxRows,
		iMinColumns, iMinRows, iMaxColumns, iMaxRows;
	//			pivoting;// default 1 (with pivoting)
	double	conditionError, solutionPrecision, variableConditionValue,
		solutionNumber, solutionNumberPrecision,
		variableSolutionNumberValue,
		meanColumns, meanRows;
	char* fileM;

	BzzVector rowsCoefficients, columnsCoefficients;
	BzzVector diag, sum;
	BzzVectorArray vdaR, vdaL;
	BzzMatrixSparseLockedByRows MLR;
	BzzVectorInt ip, numR, jp, numL, iCol, dimC;
	BzzVectorIntArray viaR, viaL;
	BzzVectorIntArray columnRows;
	BzzVectorInt ordRows, ordColumns;
	BzzVectorInt originalRowsOrder;
	void BzzBalance(BzzVector* b);
	void BzzOrdering(void);
	void GetRowsSum(void);
	void FactorizationWithPivoting(void);
	//	void FactorizationWithoutPivoting(void);
	void SolutionWithPivoting(BzzVector* bx);
	void TransposeSolutionWithPivoting(BzzVector* bx);
	//	void SolutionWithoutPivoting(BzzVector *bx);
	void CountElements(void) { numElements = numR.Sum(); }
public:
	int Rows(void) { return numRows; }
	int Columns(void) { return numColumns; }
	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzFactorizedSparseArrayGauss v;
	BzzFactorizedSparseArrayGauss(void);

	// from file
	BzzFactorizedSparseArrayGauss(char* file);
	void operator()(char* file);
	// from BzzIntArray and BzzVectorArray objects
	BzzFactorizedSparseArrayGauss(BzzVectorIntArray* ia,
		BzzVectorArray* id);
	void operator()(BzzVectorIntArray* ia, BzzVectorArray* id);
	// from BzzIntArray and ()

	void operator()(BzzVectorInt& ia, int nColumns);
	double& operator()(int i, int j);

	// from BzzMatrix
	BzzFactorizedSparseArrayGauss(BzzMatrix* A);
	void operator()(BzzMatrix* A);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzFactorizedSparseArrayGauss(void) {};

	// ============================================================================
	// =================< Functions for linear system solution >===================
	// ============================================================================
	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void)
	{
		if (factorizationStatus == FACTORIZED)
			MLR.BzzPrintStructure();
	}

	void NonOrdering(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			ordering = 0;
	}
	void Ordering(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			ordering = 1;
	}

	void Save_ftx_FileFormatted(char* file);
	void Save_btx_FileUnformatted(char* file);
	void Load_ftx_FileFormatted(char* file);
	void Load_btx_FileUnformatted(char* file);
	int GetInitialElementsNumber(void) { return numElements; }
	int GetElementsNumberInRightMatrix(void) { return viaR.GetTotalSize(); }
	int GetElementsNumberInLeftMatrix(void) { return viaL.GetTotalSize(); }
	double ConditionError(int* im = 0, double* maxE = 0);
	double SolutionPrecision(void);
	int NumberOfCorrectDigits(void);
	double SolutionNumber(BzzVector& x, int* im = 0, double* maxE = 0);
	friend void Solve(BzzFactorizedSparseArrayGauss* G, BzzVector* bx);
	friend void TransposeSolve(BzzFactorizedSparseArrayGauss* G, BzzVector* bx);
	/////////////////////////////
	void Analysis(void);
	void PrintAnalysis(void);
	/////////////////////////////
};

#endif // BZZ_BZZ_FACTORIZED_DOUBLE_SPARSE_ARRAY_GAUSS_HPP