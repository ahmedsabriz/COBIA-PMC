// BZZMATH: Release 7.0
// TODO: Condition Number

//	==================< BzzFactorizedSparseGauss.hpp >======================
//	* BzzFactorizedSparseGauss: Class for solution of linear sparse			*
//	* systems																						*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\			 				*
// *				FactorizedSparseGauss\FactorizedSparseGauss.cpp			*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1994	Date Written.
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.
//	04-1997	Added OrderingGauss.
//	05-1997	Added BzzPrintStructure.
//	07-1997	Added TransposeSolve.
//	07-1999	Modified Factorization.

/////////// Release 4.0
//	04-2000	Added ReleaseAllElementsInMatrix.
//	04-2000	Added InsertStoredElement.

/////////// Release 5.0
//	09-2003	Added Solve and TransposeSolve for matrices.
//	01-2007	Added ConditionError function;
//	01-2007	Added SolutionPrecision function;
//	01-2007	Added NumberOfCorrectDigits function;

//	============================================================================
//	******* Functions for solving linear sparse systems:								*
//	* Solve(&A,&bx);																				*
//	* Solve(&A,&BX);																				*
//	* SolveAndDestroy(&A,&bx);																	*
//	* TransposeSolve(&A,&bx);																	*
//	* TransposeSolve(&A,&BX);																	*
//	****************************************************************************
//	* The matrix A is DESTROYED (note that the argument is &A) when you use		*
// * both the functions Solve or SolveAndDestroy.										*
// * In the first case the factorization is overwritten on the original			*
// * matrix while in the second the matrix is destroyed.								*
// * Only with the function Solve it is possible to solve more systems.			*
//	****************************************************************************
//	* The coefficients of the matrix are accessible if									*
// * matrixStatus != DESTROYED																*
//	* float val = A(i,j);																		*
//	****************************************************************************
//	* The matrix can be modified only if matrixStatus != DESTROYED					*
//	* A(i,j) = val;																				*
//	****************************************************************************

#ifndef BZZ_FACTORIZEDSPARSE_DOUBLE_HPP
#define BZZ_FACTORIZEDSPARSE_DOUBLE_HPP

class BzzMatrixSparse;
struct ElementBzzMatrixSparse;

//	============================================================================
//	==================< class BzzFactorizedSparseGauss >====================
//	============================================================================

class BzzFactorizedSparseGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================

	friend void SolveAndDestroy
	(BzzFactorizedSparseGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSparseGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSparseGauss* A, BzzMatrix* BX);

	friend void TransposeSolve
	(BzzFactorizedSparseGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedSparseGauss* A, BzzMatrix* BX);

private:
	enum SparseFactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;
	enum SparseBzzMatrixStatus
	{
		AVAILABLE,
		DESTROYED,
		MODIFIED
	}matrixStatus;

	enum SparseFactorizationOrdering
	{
		ORDERING,
		NON_ORDERING,
		NON_REORDERING
	}orderingStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // per whoAmI
	int numRows, numColumns;
	ElementBzzMatrixSparse** elementRow;
	BzzMatrixSparse L;
	int* indx;
	int signd;
	int	iConditionError,
		iNumberOfCorrectDigits,
		iVariableCondition;
	double conditionError, solutionPrecision, variableConditionValue;

	int whoAmI;
	double norm;
	BzzVectorInt ordRows, ordColumns;
	BzzVector sum;

	double& InsertElement
	(ElementBzzMatrixSparse* elem, int row,
		int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSparse& rval);

	void BzzFactorizedInitialize(void);
	void Deinitialize(void);
	void BzzFactorizedDeinitialize(void);
	void PrepCopy(int rRows, int rColumns);
	BzzFactorizedSparseGauss(char, int rows, int columns);
	void PrepOnlySolve(void);
	void PrepSolve(void);
	void Norm(void);
	void Factorization(void);
	void Solution(BzzVector* bx);
	void TransposeSolution(BzzVector* bx);
	void OrderingGauss(void);
	void GetRowsSum(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseGauss(void);

	// copy constructor
	BzzFactorizedSparseGauss(const BzzFactorizedSparseGauss& rval);

	// constructor from BzzMatrixSparse
	BzzFactorizedSparseGauss(const BzzMatrixSparse& rval);

	// sizing constructor
	BzzFactorizedSparseGauss(int rows, int columns);

	// FILE ASCII
	BzzFactorizedSparseGauss(char* filematrix);
	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedSparseGauss(void);

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
		// number of rows
	int Rows(void) const
	{
		return numRows;
	}

	// number of columns
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }

	// assigns and receives vector values with control
	double& operator ()
		(int row, int column);

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	BzzFactorizedSparseGauss& operator =
		(const BzzMatrixSparse& rval);

	// transforms a BzzMatrixSparse in BzzFactorizedSparseGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseGauss* rval);

	//	============================================================================
	//	==============================< OPERATIONS >================================
	//	============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================

	friend BzzFactorizedSparseGauss operator +
		(const BzzFactorizedSparseGauss& lval,
			const BzzFactorizedSparseGauss& rval);

	BzzFactorizedSparseGauss& operator +=
		(const BzzFactorizedSparseGauss& rval);

	//	============================================================================
	//	****************************< Difference >**********************************
	//	============================================================================

	friend BzzFactorizedSparseGauss operator -
		(const BzzFactorizedSparseGauss& lval,
			const BzzFactorizedSparseGauss& rval);

	BzzFactorizedSparseGauss& operator -=
		(const BzzFactorizedSparseGauss& rval);

	//	============================================================================
	//	*******************************< Product >**********************************
	//	============================================================================
	friend void Product
	(const BzzFactorizedSparseGauss& lval, const BzzVector& rval,
		BzzVector* result);

	friend BzzVector operator *
		(const BzzFactorizedSparseGauss& lval, const BzzVector& rval);

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	void RemoveElement(int row, int column);
	void RemoveAllElementsInRow(int row);
	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzFactorizedSparseGauss* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	void CleanMatrix(double eps); // eliminates those <= eps
	void ReleaseAllElementsInMatrix(BzzMatrixSparseStore& store);
	void InsertStoredElement(int row, int column, double value,
		BzzMatrixSparseStore& store);
	void InsertStoredElement(int row, int column, double value,
		ElementBzzMatrixSparse* elem);
	double Determinant(void);
	double ConditionNumber(void);
	char Singular(void); // 1 if singular 0 if not
	int CountElements(void);
	void BzzPrintExistingElements(void);
	void BzzPrintStructure(void);
	void NonOrdering(void);
	void NonReordering(void);
	friend void ChangeDimensions(int rows, int columns,
		BzzFactorizedSparseGauss* A);

	//	============================================================================
	//	===================< System solution functions >============================
	//	============================================================================
	void SolveRight(BzzVector* bx);
	void SolveLeft(BzzVector* bx);
	void TransposeSolveRight(BzzVector* bx);
	void TransposeSolveLeft(BzzVector* bx);
	//	double ConditionError(void);
	double ConditionError(int* im = 0, double* maxE = 0);
	double SolutionPrecision(void);
	int NumberOfCorrectDigits(void);
};

#endif // BZZ_FACTORIZEDSPARSE_DOUBLE_HPP