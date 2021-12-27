// BZZMATH: Release 7.0

///////////////////////////////////////
// TODO:
// GetBzzMatrixN when singular or without NonOrdering!!!!!!!!!!!!!!!
// ConditionNumber
// NProduct
// NTProduct
// ProductN
// ProductNT
///////////////////////////////////////

//	====================< BzzFactorizedSparseLQ.HPP >===============================
//	* BzzFactorizedSparseLQ: Class for solution of square and					*
//	* underdimensioned linear sparse systems												*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 10, 11)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\							*
// *				FactorizedSparseLQ\FactorizedSparseLQ.cpp	 				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2000	Date Written.

////////////////// Release 5.0
// 01-2007	Added IntrinsicConditionNumber function.

//	============================================================================
//	******* Functions for solving square and underdimensioned linear				*
//	* sparse systems:																				*
//	* Solve(&A,&bx);																				*
//	*																									*
//	****************************************************************************
//	* The matrix A is DESTROYED (note that the argument is &A) when you use		*
// * both the functions Solve or SolveAndDestroy.										*
// * In the first case the factorization is overwritten on the original			*
// * matrix while in the second the matrix is destroyed.								*
// * Only with the function Solve it is possible to solve more systems.			*
//	****************************************************************************
//	* The coefficients of the matrix are accessible if									*
// * matrixStatus != DESTROYED																*
//	* double val = A(i,j);																		*
//	****************************************************************************
//	* The matrix can be modified only if matrixStatus != DESTROYED					*
//	* A(i,j) = val;																				*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_SPARSE_LQ_HPP
#define BZZ_FACTORIZED_DOUBLE_SPARSE_LQ_HPP

class BzzMatrixSparse;
struct ElementBzzMatrixSparse;

//	============================================================================
//	=====================< class BzzFactorizedSparseLQ >==========================
//	============================================================================
class BzzFactorizedSparseLQ : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	//friend void SolveAndDestroy
	//		 (BzzFactorizedSparseLQ *A,BzzVector *bx,double precision = 10.F*MachEpsFloat());
	friend void Solve
	(BzzFactorizedSparseLQ* A, BzzVector* bx);//,double precision = 10.*MachEpsFloat());
	friend void TransposeSolve
	(BzzFactorizedSparseLQ* A, BzzVector* bx);//,double precision = 10.*MachEpsFloat());

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
	enum SparseFactorizationOrderingLQ
	{
		ORDERING,
		NON_ORDERING,
		NON_REORDERING
	}orderingStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	int numRows, numColumns;
	ElementBzzMatrixSparse** elementRow;
	BzzVector d;
	char singular;	// singular = 1 when singular
	double sinTiny;
	int whoAmI;
	BzzVector norm;
	BzzVectorInt ordRows, ordColumns;

	double& InsertElement
	(ElementBzzMatrixSparse* elem, int row, int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSparse& rval);

	void BzzFactorizedInitialize(void);
	void Deinitialize(void);
	void BzzFactorizedDeinitialize(void);
	void PrepCopy(int rRows, int rColumns);
	BzzFactorizedSparseLQ(char, int rows, int columns);
	void PrepOnlySolve(void);
	void PrepSolve(void);
	void Norm(void);
	double Condition(void);
	void Factorization(void);
	void Solution(BzzVector* bx);
	void TransposeSolution(BzzVector* bx);
	void OrderingLQ(void);
	//	void SolveRight(BzzVector *bx,double mach);
	//	void SolveLeft(BzzVector *bx);
	int iIntrinsicConditionNumber;
	double intrinsicConditionNumber;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseLQ(void);

	// copy constructor
	BzzFactorizedSparseLQ(const BzzFactorizedSparseLQ& rval);

	// constructor from BzzMatrixSparse
	BzzFactorizedSparseLQ(const BzzMatrixSparse& rval);

	// sizing constructor
	BzzFactorizedSparseLQ(int rows, int columns);

	// FILE ASCII
	BzzFactorizedSparseLQ(char* filematrix);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedSparseLQ(void);

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
	BzzFactorizedSparseLQ& operator =
		(const BzzMatrixSparse& rval);

	// transforms a BzzMatrixSparse in BzzFactorizedSparseLQ
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseLQ* rval);

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	BzzVectorInt GetLinearCombinations(void);
	void PrintLinearCombinations(void);
	int GetTheNumberOfIndependentEquations(void);
	void GetBzzMatrixN(BzzMatrix* N);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	void SetPrecision(double precision);
	void RemoveElement(int row, int column);
	void RemoveAllElementsInRow(int row);
	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzFactorizedSparseLQ* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	void CleanMatrix(double eps); // eliminates those <= eps
	double Determinant(void);
	double ConditionNumber(void);
	double IntrinsicConditionNumber(void)
	{
		if (iIntrinsicConditionNumber != 0)
			return intrinsicConditionNumber;
		else
			return 0.;
	}
	char Singular(void); // 1 if singular 0 if not
	int CountElements(void);
	void BzzPrintExistingElements(void);
	void BzzPrintStructure(void);
	void NonOrdering(void);
	void NonReordering(void);
	void FindBands(int* low, int* up);
};

#endif // BZZ_FACTORIZED_DOUBLE_SPARSE_LQ_HPP