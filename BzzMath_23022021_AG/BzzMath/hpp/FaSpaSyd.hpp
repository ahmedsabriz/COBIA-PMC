// TODO: Condition Number
// TODO: BzzPrintStructure

// BZZMATH: Release 3.1 version alfa

//	=============================< FASPASYD.HPP >===============================
//	* BzzFactorizedSparseCholesky:														*
//	* Class for sparse, symmetric and definite positive matrices					*
//	* BzzFactorizedSparseSymmetric:														*
//	* Class for sparse, symmetric matrices													*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\ExFaSpSy.cpp							 				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	03-1997	Date Written.
//	03-1998	Added OrderingSymmetric function.
//	03-1998	Added BzzFactorizedSparseSymmetric class.

//	============================================================================
//	******* Functions for solving linear sparse systems:								*
//	* Solve(&A,&bx);																				*
//	* SolveAndDestroy(&A,&bx);																	*
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
//	============================================================================
//	******	BzzFactorizedSparseCholesky constructors:							*
//	* BzzFactorizedSparseCholesky A; // default										*
//	* BzzFactorizedSparseCholesky A = B; // copy-initializer					*
//	* BzzFactorizedSparseCholesky A(n,n); // sizes									*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A(i,j);																					*
//	* xf = A.GetValue(i,j);																		*
//	* A(i,j) = xf;																					*
//	* A.SetValue(i,j,xf);																		*
//	* char t = A.Existing(i,j);																*
//	* int count = BzzMatrixSparse::ObjectCount();										*
//	* int countInScope = BzzMatrixSparse::ObjectCountInScope();						*
//	****************************************************************************
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixSparse																*
//	****************************************************************************
//	***** BzzMessage and BzzPrint																*
//	* A.BzzPrint("Comment"); // always														*
//	* A.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	****************************************************************************
//	***** Remove, Delete, Clean, ChangeDimensions and Swap							*
//	* A.RemoveElement(3,151);																	*
//	* A.RemoveAllElementsInMatrix();		// eliminates all elements					*
//	* A.CleanMatrix(eps);	// eliminates those <= eps									*
//	* Delete(&A);																					*
//	* ChangeDimensions(newr,newc,&A);														*
//	* Swap(&A,&B);																					*
//	****************************************************************************
//	***** Function for solution of symmetrical sparse systems						*
//	* A.Solve(A,B,&x);																			*
//	****************************************************************************

#ifndef FACTORED_DOUBLE_SPARSE_CHOLESKY
#define FACTORED_DOUBLE_SPARSE_CHOLESKY

class BzzMatrixSparseSymmetric;
struct ElementBzzMatrixSparse;

//	============================================================================
//	=================< class BzzFactorizedSparseCholesky >==================
//	============================================================================
class BzzFactorizedSparseCholesky : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================

	friend void SolveAndDestroy
	(BzzFactorizedSparseCholesky* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSparseCholesky* A, BzzVector* bx);
protected:
	enum SparseCholeskyFactorizationStatus
	{
		UNFACTORED,
		FACTORED
	}factorizationStatus;
	enum SparseCholeskyMatrixStatus
	{
		AVAILABLE,
		DESTROYED,
		MODIFIED
	}matrixStatus;

	enum SparseSymmetricFactorizationOrdering
	{
		ORDERING,
		NON_ORDERING,
		NON_REORDERING
	}orderingStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;
	int numRows, numColumns;
	int whoAmI;

	ElementBzzMatrixSparse** elementRow;
	double norm;
	BzzVector diagonal;
	BzzVectorInt ordRows;

	double& InsertElement
	(ElementBzzMatrixSparse* elem, int row, int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSparseSymmetric& rval);
	void BzzFactorizedInitialize(void);
	void Deinitialize(void);
	void BzzFactorizedDeinitialize(void);
	void PrepCopy(int rRows, int rColumns);
	void PrepOnlySolve(void);
	void PrepSolve(void);
	void Norm(void);
	double Condition(void);
	virtual void Factorization(void);
	virtual void Solution(BzzVector* bx);
	void OrderingSymmetric(void);
	int eigenvalues;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseCholesky(void);

	// copy constructor
	BzzFactorizedSparseCholesky(const BzzFactorizedSparseCholesky& rval);

	// constructor from BzzMatrixSparseSymmetric
	BzzFactorizedSparseCholesky(const BzzMatrixSparseSymmetric& rval);

	// sizing constructor
	BzzFactorizedSparseCholesky(int rows, int columns);

	// FILE ASCII
	BzzFactorizedSparseCholesky(char* filematrix);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFactorizedSparseCholesky(void);

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
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// assigns and receives vector values with control
	double& operator ()
		(int row, int column);

	void SetValue(int i, int j, double val)
	{
		(*this)(i, j) = val;
	}
	double GetValue(int i, int j);
	char Existing(int i, int j);

	//	============================================================================
	//	****************************< assignment operators >************************
	//	============================================================================
	BzzFactorizedSparseCholesky& operator =
		(const BzzMatrixSparseSymmetric& rval);

	// transforms a BzzMatrixSparse in BzzFactorizedSparseGauss
// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparseSymmetric* lval, BzzFactorizedSparseCholesky* rval);

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================

	//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	============================< Modifying Functions >=========================
	//	============================================================================
	void RemoveElement(int row, int column);
	friend void Delete(BzzFactorizedSparseCholesky* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	void CleanMatrix(double eps); // eliminates those <= eps
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
};

//	============================================================================
//	====================< class BzzFactorizedSymmetric >====================
//	============================================================================

class BzzFactorizedSparseSymmetric : public BzzFactorizedSparseCholesky
{
private:
	BzzVectorInt kMemo, pivot;
	BzzVector lMemo, rMemo;
	char pos, neg, nul;
	void RemoveDiagonal(void);
	void RemoveRowColumn(int row);

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseSymmetric(void)
		: BzzFactorizedSparseCholesky() {};

	// copy constructor
	BzzFactorizedSparseSymmetric(const BzzFactorizedSparseCholesky& rval)
		: BzzFactorizedSparseCholesky(rval) {};

	// constructor from BzzMatrixSparseSymmetric
	BzzFactorizedSparseSymmetric(const BzzMatrixSparseSymmetric& rval)
		: BzzFactorizedSparseCholesky(rval) {};

	// sizing constructor
	BzzFactorizedSparseSymmetric(int rows, int columns)
		: BzzFactorizedSparseCholesky(rows, columns) {};

	// FILE ASCII
	BzzFactorizedSparseSymmetric(char* filematrix)
		: BzzFactorizedSparseCholesky(filematrix) {};

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	//	void operator = (const BzzFactorizedSymmetric &rval);
	//	BzzFactorizedSymmetric &operator = (const BzzMatrixSymmetric &rval);
	int Eigenvalues(void) { return eigenvalues; }

	//	============================================================================
	//	=======================< Functions for linear algebra >=====================
	//	============================================================================
	virtual void Factorization(void);
	virtual void Solution(BzzVector* bx);

	//	virtual void Factorization(void);
	//	virtual void Solution(int n,double **a,double *b);
	//	virtual void ChangeFactorization
	//		(int n,double **a,double alfa,double *z);
	//	double ConditionNumber(void);
	double Determinant(void);
};

#endif // FACTORED_DOUBLE_SPARSE_CHOLESKY