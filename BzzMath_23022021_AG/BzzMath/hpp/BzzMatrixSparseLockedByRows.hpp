// Serve per fare il prodotto fra una matrice sparsa ed un vettore
// nel caso in cui la matrice non debba essere ulteriormente modificata
// la matrice viene memorizzata in modo piu' efficiente

// BZZMATH: Release 7.0

//	==================< BzzMatrixSparseLockedByRows.hpp >=======================
//	* Class BzzMatrixSparseLockedByRows for operations									*
// *	between matrices and vectors															*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixSparseLocked\MatrixSparseLocked.cpp				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1996	Date Written.

/////////// Release 4.0
//	04-2000	Added operator = Sparse.
//	04-2000	Added TProduct.
//	04-2000	Added GetRowsNorm2.
//	04-2000	Added GetColumnsNorm2.
//	05-2000	Added ProductForSelectedRowsAndColumns.
//	05-2000	Added operator =.
//	05-2000	Added CopyMatrixForSelectedRowsAndColumns.
//	05-2000	Added ProductForSelectedColumns.
//	05-2000	Added TProductForSelectedColumns.
//	05-2000	Added CopyMatrixForOrderedRowsAndSelectedColumns.
//	05-2000	Added CopyMatrixForselectedColumns.
//	06-2000	Bug fixed in TProduct.
//	06-2000	Added GetNumEquationsForEachVariable.
//	06-2000	Added GetNumVariablesForEachEquation.
//	09-2000	Added OrderingLQ.
//	10-2000	Added Division.
//	11-2000	Added GetNumVariablesForEachEquationAndNumEquationsForEachVariable.
//	11-2000	Added GetEquationsContainingSingletons.
//	11-2000	Added GetAnalysis.
//	11-2000	Added BzzPrintStructure.
//	03-2002	Added SolveLeft and SolveRight.
//	03-2002	Added constructor from FILE.
//	07-2002	Added Product between two BzzMatrixSparseLockedByRows.

/////////// Release 5.0
//	11-2003	Added GetRow function.
//	11-2003	Added GetColumn function.
//	11-2003	Added GetRowForSelectedColumns function.
//	11-2003	Added GetColumnForSelectedRows function.
//	10-2004	Added another Product function.
//	11-2004	Added ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked function.
//	11-2004	Added Minus function.
//	12-2004	Added ProductForSelectedRow function.
//	12-2004	Added ResetLockedCoefficients function.
//	05-2006	Added ProductAndDifference function.
//	05-2006	Added ProductAndDifferenceWithControl function.
//	12-2006	Added GetRowsSum function.
//	12-2006	Added GetColumnsSum function.
//	10-2007	Added operator operator ()(BzzVectorArray &V,BzzVectorIntArray &K).
// 11-2007  Added ProductForSelectedRows function.
// 11-2007  Added ProductFromRowjToRowk function.

/////////// Release 6.0
//	11-2008	Added compact constructors.
//	11-2008	Added operator(i,j).
// 04-2012	Added BzzMatrixSparseLockedByRowsAndByColumns class.

//	============================================================================
//	******	BzzMatrixSparseLockedByRows constructors:									*
//	* BzzMatrixSparseLockedByRows A; // default											*
//	* BzzMatrixSparseLockedByRows A = B; // from Sparse								*
//	* BzzMatrixSparseLockedByRows A("SPARSE:DAT"); // from FILE						*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A.GetValue(i,j);																		*
//	* int count = BzzMatrixSparseLockedByRows::ObjectCount();						*
//	* int countInScope = BzzMatrixSparseLockedByRows::ObjectCountInScope();		*
//	* double *ptrVal = A.Scanning(&i,&j,&val);											*
//	* A.BeginScanning();																			*
//	****************************************************************************
//	***** Assignments:																			*
//	* ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked(&S,&A);						*
//	* ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked(&S,&A);					*
//	* BzzMatrixSparseLockedByRows A;													*
//	* A = B; // from Sparse																		*
//	* A.ResetLockedCoefficients(B) // from Sparse										*
//	****************************************************************************
//	***** Functions:																				*
//	A.GetRowsNorm2(&nrom2R);																	*
//	A.GetColumnsNorm2(&nrom2C);																*
//	****************************************************************************
//	***** Implemented operations :															*
//	* Product(A,x,&y);	// y = A*x;															*
//	* Product(A,B,&C);	// C = A*B;															*
//	* ProductAndDifference(d,D,x,r); // r = d - D*x										*
//	* int i = ProductAndDifferenceWithControl(d,D,x,r); // r = d - D*x			*
//	* i = 0 if r <= 0																				*
//	* i = j if r[j] > 0																			*
//	* TProduct(A,v,&y);	// y = AT*x;														*
//	* ProductForSelectedRow(A,B,row,&v);													*
//	* ProductForSelectedRowsAndColumns(A,selR,selC,v,&y);								*
//	* ProductForSelectedRows(A,selR,v,&y);													*
//	* ProductFromRowjToRowk(A,j,k,v,&y);
//	* ProductForSelectedColumns(A,selC,v,&y);												*
//	* TProductForSelectedColumns(A,selC,v,&y);											*
//	* Division(&A,2.);																			*
//	* Division(&A,x);																				*
//	****************************************************************************

#ifndef BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_BY_ROWS_HPP
#define BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_BY_ROWS_HPP

//	============================================================================
//	==================< class BzzMatrixSparseLockedByRows >=====================
//	============================================================================
class BzzVectorArray;
class BzzMatrixSparseLockedByRows : public BzzBaseClass
{
	friend class BzzMatrixSparse;
	friend class BzzLinearProgramming;
	friend class BzzFactorizedLQ;
	friend class BzzFactorizedSparseLockedLQ;
	friend class BzzLinearProgrammingAttic;
	friend void ProLockLock(int iStart, int iEnd, BzzMatrixSparseLockedByRows& A, BzzMatrixSparseLockedByRows& b,
		BzzMatrix* C);
	friend void ProLockFull(int iStart, int iEnd, BzzMatrixSparseLockedByRows& A, BzzMatrix& B,
		BzzMatrix* C);

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double BZZ_ERROR_PRECISION;
	static const int BZZ_ITER_TOTAL;

	double* matrix;
	int* columns;
	int* coeffInRow;
	int* ptrRow;
	int numRows, numColumns;
	int size;
	int whoAmI;
	//	char shadow;
	BzzVectorInt ordSelectedColumns;
	BzzMatrixInt IJ;

	int iterativeSolve;
	double* matrixLeft, * matrixRight;
	int* columnsLeft, * columnsRight;
	int* coeffInRowLeft, * coeffInRowRight;
	//	int *ptrRowLeft,*ptrRowRight;
	int sizeLeft, sizeRight;
	BzzVector diagonal;
	void InitializeIterative(void);

	// initialisation constructors
	void Initialize(int mrows, int ncolumns, int nsize);

	// deinitialisation
	void Deinitialize(void);
	// prepares assignments
//	void PrepCopy(int rRows,int rColumns);
	void BuildIJ(void);

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default // BzzMatrixSparseLockedByRows A;
	BzzMatrixSparseLockedByRows(void);

	// copy-initializer // BzzMatrixSparseLockedByRows A = B;
	BzzMatrixSparseLockedByRows(BzzMatrixSparseLockedByRows& E);

	// from BzzMatrixSparse
	BzzMatrixSparseLockedByRows(BzzMatrixSparse& rval);

	// from FILE
	BzzMatrixSparseLockedByRows(char* filematrix);

	void operator ()(BzzVectorArray& V, BzzVectorIntArray& K, int nC);

	BzzMatrixSparseLockedByRows(int nr, int nc, BzzVector* v, BzzVectorInt* c,
		BzzVectorInt* r);

	BzzMatrixSparseLockedByRows(int nr, int nc, BzzVectorInt* c,
		BzzVectorInt* r);

	BzzMatrixSparseLockedByRows(int nr, int nc, int nt, BzzVector* v,
		BzzVectorInt* r, BzzVectorInt* c);

	//	============================================================================
	//	****************************< destructor >**********************************
	//	============================================================================
	~BzzMatrixSparseLockedByRows(void)
	{
		countInScope--;Deinitialize();
	}

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
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
	double* Scanning(int* i, int* j, double* val);
	void BeginScanning(void);
	void BzzPrintStructure(void);

	//	============================================================================
	//	**********************************< Access functions >**********************
	//	============================================================================
		// assigns and receives vector values with control
	double& operator () (int row, int col);

	//	============================================================================
	//	**********************< assignment operators >******************************
	//	============================================================================

		// transforms a BzzMatrixSparse in BzzMatrixSparseLockedByRows
		// BzzMatrixSparse gets destroyed
	friend void ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked
	(BzzMatrixSparse* lval, BzzMatrixSparseLockedByRows* rval);
	friend void ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked
	(BzzMatrixSparse* lval, BzzMatrixSparseLockedByRows* rval);
	friend void ReplaceBzzMatrixWithBzzMatrixSparseLocked
	(BzzMatrix* lval, BzzMatrixSparseLockedByRows* rval);

	friend void Swap(BzzMatrixSparseLockedByRows* lval,
		BzzMatrixSparseLockedByRows* rval);

	BzzMatrixSparseLockedByRows& operator =
		(BzzMatrixSparse& S);
	BzzMatrixSparseLockedByRows& operator =
		(BzzMatrixSparseLockedByRows& E);

	void ResetLockedCoefficients(BzzMatrixSparse& S);
	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByRows* G);
	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByRows& E,
		BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByRows* G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByRows* G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& selectedColumns,
		BzzMatrix* G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& selectedColumns,
		BzzMatrixSparseStore& store,
		BzzMatrixSparse* G);
	friend void CopyMatrixForOrderedRowsAndColumns
	(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& orderedRows, BzzVectorInt& orderedColumns,
		BzzMatrix* G);
	friend void CopyMatrixForOrderedRowsAndColumns
	(BzzMatrixSparseLockedByRows& P,
		BzzVectorInt& orderedRows, BzzVectorInt& orderedColumns,
		BzzMatrixSparseStore& store,
		BzzMatrixSparse* G);
	friend void CopyMatrixForLinProMat(BzzMatrixSparse& E,
		BzzMatrixSparse& D, BzzMatrixSparseLockedByRows* W);
	friend void CopyMatrixForLinProMat(BzzMatrixSparseLockedByRows& E,
		BzzMatrixSparseLockedByRows& D, BzzMatrixSparseLockedByRows* W);
	friend void CopyMatrixForLinProMat(BzzVectorInt& ie,
		BzzMatrixSparseLockedByRows& E,
		BzzMatrixSparseLockedByRows& D, BzzMatrixSparseLockedByRows* W);

	//	============================================================================
	//	*****************************< Other functions >****************************
	//	============================================================================
	void GetRowsNorm2(BzzVector* norm2R);
	void GetColumnsNorm2(BzzVector* norm2R);
	void GetRowsSum(BzzVector* sumR);
	void GetColumnsSum(BzzVector* sumC);
	void GetRowsNorm2X(BzzVector& x, BzzVector* norm2X);
	void GetRowsWeieghtedX(BzzVector& x, BzzVector* norm2X);
	void GetRow(int i, BzzVector* x);
	void GetColumn(int i, BzzVector* x);
	void GetRowForSelectedColumns(int i, BzzVectorInt& iv, BzzVector* x);
	void GetColumnForSelectedRows(int i, BzzVectorInt& iv, BzzVector* x);
	void GetNumEquationsForEachVariable
	(BzzVectorInt* numEquationsForEachVariable);
	void GetNumVariablesForEachEquation
	(BzzVectorInt* numVariablesForEachEquation);
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable);
	void GetEquationsContainingSingletons
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable,
		BzzVectorInt* equationsContainingSingletons);
	void GetAnalysis
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable,
		BzzVectorInt* equationsContainingSingletons,
		BzzVectorInt* equationsContainingMoreSingletons,
		BzzVectorInt* twoVariablesEquationsContainingSingletons,
		BzzVectorInt* threeVariablesEquationsContainingSingletons);

	int OrderingLinearProgramming(BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);

	int ChoiceLinearProgrammingVariables(BzzVector& unorm2Columns,
		BzzVectorInt* varForFactorization);
	void OrderingLQ(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);

	//	============================================================================
	//	=====================================	OPERATIONS	=======================
	//	============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================

	//	============================================================================
	//	*****************************< Difference >*********************************
	//	============================================================================

	//	============================================================================
	//	*********************************< Minus >**********************************
	//	============================================================================
	void Minus(void);

	//	============================================================================
	//	*******************************< Product >**********************************
	//	============================================================================
	friend void Product(double d, BzzMatrixSparseLockedByRows& P,
		BzzMatrixSparseLockedByRows* R);
	friend void Product(BzzMatrixSparseLockedByRows& P,
		BzzVector& x, BzzVector* y);
	friend void Product(BzzMatrixSparseLockedByRows& A,
		BzzMatrixSparseLockedByRows& B, BzzMatrix* C);
	friend void Product(BzzMatrixSparseLockedByRows& A,
		BzzMatrix& B, BzzMatrix* C);

	friend void ProductAndDifference(BzzVector& d,
		BzzMatrixSparseLockedByRows& P, BzzVector& x,
		BzzVector* y);
	friend int ProductAndDifferenceWithControl(BzzVector& d,
		BzzMatrixSparseLockedByRows& P, BzzVector& x,
		BzzVector* y);

	friend void TProduct(BzzMatrixSparseLockedByRows& P,
		BzzVector& x, BzzVector* y);

	friend void ProductForSelectedRow(BzzMatrixSparseLockedByRows& A,
		BzzMatrix& B, int row, BzzVector* c);
	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByRows& D,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzVector& x, BzzVector* y);
	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByRows& E,
		BzzMatrixSparseLockedByRows& D,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzVector& x, BzzVector* y);
	friend void ProductForSelectedRows(BzzMatrixSparseLockedByRows& D,
		BzzVectorInt& selectedRows, BzzVector& x, BzzVector* y);
	friend void ProductFromRowjToRowk(BzzMatrixSparseLockedByRows& D,
		int j, int k, BzzVector& x, BzzVector* y);
	friend void ProductForSelectedColumns(BzzMatrixSparseLockedByRows& D,
		BzzVectorInt& selectedColumns, BzzVector& x, BzzVector* y);
	friend void TProductForSelectedColumns(BzzMatrixSparseLockedByRows& D,
		BzzVectorInt& selectedColumns, BzzVector& x, BzzVector* y);

	//	============================================================================
	//	******************************< TProduct >**********************************
	//	============================================================================

	//	============================================================================
	//	******************************< ProductT >**********************************
	//	============================================================================

	//	============================================================================
	//	*******************************< Division >*********************************
	//	============================================================================
	friend void Division(BzzMatrixSparseLockedByRows* P, double x);
	friend void Division(BzzMatrixSparseLockedByRows* P, BzzVector& x);

	//	============================================================================
	//	=========================	Triangular Systems Solution	===================
	//	============================================================================
	friend void SolveLeft(BzzMatrixSparseLockedByRows& P,
		BzzVector& x, BzzVector* y);
	friend void SolveLeft(BzzMatrixSparseLockedByRows& P, BzzVector* y);
	friend void SolveRight(BzzMatrixSparseLockedByRows& P,
		BzzVector& x, BzzVector* y);
	friend void SolveRight(BzzMatrixSparseLockedByRows& P, BzzVector* y);

	//	============================================================================
	//	===========================	Iterative Systems Solution	===================
	//	============================================================================
	int JacobiSolve(BzzVector& b, BzzVector* x);
	int GaussSeidelForwardSolve(BzzVector& b, BzzVector* x);
	int GaussSeidelBackwardSolve(BzzVector& b, BzzVector* x);
	int GaussSeidelSolve(BzzVector& b, BzzVector* x);
	int GaussSeidelForwardSolve(double omega, BzzVector& b, BzzVector* x);
	int GaussSeidelBackwardSolve(double omega, BzzVector& b, BzzVector* x);
	int GaussSeidelSolve(double omega, BzzVector& b, BzzVector* x);

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};
/*
//	============================================================================
//	===========< class BzzMatrixSparseLockedByRowsAndByColumns >================
//	============================================================================
class BzzVectorArray;
class BzzMatrixSparseLockedByRowsAndByColumns : public BzzBaseClass
	{
friend class BzzMatrixSparse;
friend class BzzLinearProgramming;
friend class BzzFactorizedLQ;
friend class BzzFactorizedSparseLockedLQ;
friend class BzzLinearProgrammingAttic;
friend void ProLockLock(int iStart,int iEnd,BzzMatrixSparseLockedByRowsAndByColumns &A,BzzMatrixSparseLockedByRowsAndByColumns &b,
			BzzMatrix *C);
friend void ProLockFull(int iStart,int iEnd,BzzMatrixSparseLockedByRowsAndByColumns &A,BzzMatrix &B,
			BzzMatrix *C);
	private:
	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double BZZ_ERROR_PRECISION;
	static const int BZZ_ITER_TOTAL;

	double *matrix;
	int *columns;
	int *coeffInRow;
	int *ptrRow;
	int numRows,numColumns;
	int size;
	int whoAmI;
//	char shadow;
	BzzVectorInt ordSelectedColumns;
	BzzMatrixInt IJ;

	int iterativeSolve;
	double *matrixLeft,*matrixRight;
	int *columnsLeft,*columnsRight;
	int *coeffInRowLeft,*coeffInRowRight;
//	int *ptrRowLeft,*ptrRowRight;
	int sizeLeft,sizeRight;
	BzzVector diagonal;
	void InitializeIterative(void);

	// initialisation constructors
	void Initialize(int mrows,int ncolumns,int nsize);

	// deinitialisation
	void Deinitialize(void);
	// prepares assignments
//	void PrepCopy(int rRows,int rColumns);
	void BuildIJ(void);

public:
//	============================================================================
//	***************************< constructors >*********************************
//	============================================================================
	// default // BzzMatrixSparseLockedByRowsAndByColumns A;
	BzzMatrixSparseLockedByRowsAndByColumns(void);

	// copy-initializer // BzzMatrixSparseLockedByRowsAndByColumns A = B;
	BzzMatrixSparseLockedByRowsAndByColumns(BzzMatrixSparseLockedByRowsAndByColumns &E);

	// from BzzMatrixSparse
	BzzMatrixSparseLockedByRowsAndByColumns(BzzMatrixSparse &rval);

		// from FILE
	BzzMatrixSparseLockedByRowsAndByColumns(char *filematrix);

	void operator ()(BzzVectorArray &V,BzzVectorIntArray &K,int nC);

	BzzMatrixSparseLockedByRowsAndByColumns(int nr,int nc,BzzVector *v,BzzVectorInt *c,
		BzzVectorInt *r);

	BzzMatrixSparseLockedByRowsAndByColumns(int nr,int nc,BzzVectorInt *c,
		BzzVectorInt *r);

	BzzMatrixSparseLockedByRowsAndByColumns(int nr,int nc,int nt,BzzVector *v,
		BzzVectorInt *r,BzzVectorInt *c);

//	============================================================================
//	****************************< destructor >**********************************
//	============================================================================
	~BzzMatrixSparseLockedByRowsAndByColumns(void)
		{countInScope--;Deinitialize();}

//	============================================================================
//	********************< Non-modifying access functions >**********************
//	============================================================================
	// number of rows
	int Rows(void) const
		{return numRows;}

	// number of columns
	int Columns(void) const
		{return numColumns;}

	int WhoAmI(void) const {return whoAmI;}
	static int ObjectCount(void){return count;}
	static int ObjectCountInScope(void){return countInScope;}
	double *Scanning(int *i,int *j,double *val);
	void BeginScanning(void);
	void BzzPrintStructure(void);

//	============================================================================
//	**********************************< Access functions >**********************
//	============================================================================
	// assigns and receives vector values with control
	double &operator () (int row,int col);

//	============================================================================
//	**********************< assignment operators >******************************
//	============================================================================

	// transforms a BzzMatrixSparse in BzzMatrixSparseLockedByRowsAndByColumns
	// BzzMatrixSparse gets destroyed
	friend void ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked
	(BzzMatrixSparse *lval,BzzMatrixSparseLockedByRowsAndByColumns *rval);
	friend void ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked
	(BzzMatrixSparse *lval,BzzMatrixSparseLockedByRowsAndByColumns *rval);
	friend void ReplaceBzzMatrixWithBzzMatrixSparseLocked
	(BzzMatrix *lval,BzzMatrixSparseLockedByRowsAndByColumns *rval);

	friend void Swap(BzzMatrixSparseLockedByRowsAndByColumns *lval,
					  BzzMatrixSparseLockedByRowsAndByColumns *rval);

	BzzMatrixSparseLockedByRowsAndByColumns &operator =
		 (BzzMatrixSparse &S);
	BzzMatrixSparseLockedByRowsAndByColumns &operator =
		 (BzzMatrixSparseLockedByRowsAndByColumns &E);

	void ResetLockedCoefficients(BzzMatrixSparse &S);
	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
		BzzMatrixSparseLockedByRowsAndByColumns *G);
	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByRowsAndByColumns &E,
		BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
		BzzMatrixSparseLockedByRowsAndByColumns *G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &selectedColumns,
		BzzMatrixSparseLockedByRowsAndByColumns *G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &selectedColumns,
		BzzMatrix *G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &selectedColumns,
		BzzMatrixSparseStore &store,
		BzzMatrixSparse *G);
	friend void CopyMatrixForOrderedRowsAndColumns
		(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &orderedRows,BzzVectorInt &orderedColumns,
		BzzMatrix *G);
	friend void CopyMatrixForOrderedRowsAndColumns
		(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVectorInt &orderedRows,BzzVectorInt &orderedColumns,
		BzzMatrixSparseStore &store,
		BzzMatrixSparse *G);
	friend void CopyMatrixForLinProMat(BzzMatrixSparse &E,
		BzzMatrixSparse &D,BzzMatrixSparseLockedByRowsAndByColumns *W);
	friend void CopyMatrixForLinProMat(BzzMatrixSparseLockedByRowsAndByColumns &E,
		BzzMatrixSparseLockedByRowsAndByColumns &D,BzzMatrixSparseLockedByRowsAndByColumns *W);
	friend void CopyMatrixForLinProMat(BzzVectorInt &ie,
		BzzMatrixSparseLockedByRowsAndByColumns &E,
		BzzMatrixSparseLockedByRowsAndByColumns &D,BzzMatrixSparseLockedByRowsAndByColumns *W);

//	============================================================================
//	*****************************< Other functions >****************************
//	============================================================================
	void GetRowsNorm2(BzzVector *norm2R);
	void GetColumnsNorm2(BzzVector *norm2R);
	void GetRowsSum(BzzVector *sumR);
	void GetColumnsSum(BzzVector *sumC);
	void GetRowsNorm2X(BzzVector &x,BzzVector *norm2X);
	void GetRowsWeieghtedX(BzzVector &x,BzzVector *norm2X);
	void GetRow(int i,BzzVector *x);
	void GetColumn(int i,BzzVector *x);
	void GetRowForSelectedColumns(int i,BzzVectorInt &iv,BzzVector *x);
	void GetColumnForSelectedRows(int i,BzzVectorInt &iv,BzzVector *x);
	void GetNumEquationsForEachVariable
		(BzzVectorInt *numEquationsForEachVariable);
	void GetNumVariablesForEachEquation
		(BzzVectorInt *numVariablesForEachEquation);
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable
		(BzzVectorInt *numVariablesForEachEquation,
		BzzVectorInt *numEquationsForEachVariable);
	void GetEquationsContainingSingletons
		(BzzVectorInt *numVariablesForEachEquation,
		BzzVectorInt *numEquationsForEachVariable,
		BzzVectorInt *equationsContainingSingletons);
	void GetAnalysis
		(BzzVectorInt *numVariablesForEachEquation,
		BzzVectorInt *numEquationsForEachVariable,
		BzzVectorInt *equationsContainingSingletons,
		BzzVectorInt *equationsContainingMoreSingletons,
		BzzVectorInt *twoVariablesEquationsContainingSingletons,
		BzzVectorInt *threeVariablesEquationsContainingSingletons);

	int OrderingLinearProgramming(BzzVectorInt *ordRows,
		BzzVectorInt *ordColumns);

	int ChoiceLinearProgrammingVariables(BzzVector &unorm2Columns,
		BzzVectorInt *varForFactorization);
	void OrderingLQ(BzzVectorInt *ordRows,BzzVectorInt *ordColumns);

//	============================================================================
//	=====================================	OPERATIONS	=======================
//	============================================================================

//	============================================================================
//	********************************< Sum >*************************************
//	============================================================================

//	============================================================================
//	*****************************< Difference >*********************************
//	============================================================================

//	============================================================================
//	*********************************< Minus >**********************************
//	============================================================================
	void Minus(void);

//	============================================================================
//	*******************************< Product >**********************************
//	============================================================================
	friend void Product(double d,BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzMatrixSparseLockedByRowsAndByColumns *R);
	friend void Product(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVector &x,BzzVector *y);
	friend void Product(BzzMatrixSparseLockedByRowsAndByColumns &A,
		BzzMatrixSparseLockedByRowsAndByColumns &B,BzzMatrix *C);
	friend void Product(BzzMatrixSparseLockedByRowsAndByColumns &A,
		BzzMatrix &B,BzzMatrix *C);

	friend void ProductAndDifference(BzzVector &d,
		BzzMatrixSparseLockedByRowsAndByColumns &P,BzzVector &x,
		BzzVector *y);
	friend int ProductAndDifferenceWithControl(BzzVector &d,
		BzzMatrixSparseLockedByRowsAndByColumns &P,BzzVector &x,
		BzzVector *y);

	friend void TProduct(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVector &x,BzzVector *y);

	friend void ProductForSelectedRow(BzzMatrixSparseLockedByRowsAndByColumns &A,
		BzzMatrix &B,int row,BzzVector *c);
	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByRowsAndByColumns &D,
		BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
		BzzVector &x,BzzVector *y);
	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByRowsAndByColumns &E,
		BzzMatrixSparseLockedByRowsAndByColumns &D,
		BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
		BzzVector &x,BzzVector *y);
	friend void ProductForSelectedRows(BzzMatrixSparseLockedByRowsAndByColumns &D,
		BzzVectorInt &selectedRows,BzzVector &x,BzzVector *y);
	friend void ProductFromRowjToRowk(BzzMatrixSparseLockedByRowsAndByColumns &D,
		int j,int k,BzzVector &x,BzzVector *y);
	friend void ProductForSelectedColumns(BzzMatrixSparseLockedByRowsAndByColumns &D,
		BzzVectorInt &selectedColumns,BzzVector &x,BzzVector *y);
	friend void TProductForSelectedColumns(BzzMatrixSparseLockedByRowsAndByColumns &D,
		BzzVectorInt &selectedColumns,BzzVector &x,BzzVector *y);

//	============================================================================
//	******************************< TProduct >**********************************
//	============================================================================

//	============================================================================
//	******************************< ProductT >**********************************
//	============================================================================

//	============================================================================
//	*******************************< Division >*********************************
//	============================================================================
friend void Division(BzzMatrixSparseLockedByRowsAndByColumns *P,double x);
friend void Division(BzzMatrixSparseLockedByRowsAndByColumns *P,BzzVector &x);

//	============================================================================
//	=========================	Triangular Systems Solution	===================
//	============================================================================
	friend void SolveLeft(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVector &x,BzzVector *y);
	friend void SolveLeft(BzzMatrixSparseLockedByRowsAndByColumns &P,BzzVector *y);
	friend void SolveRight(BzzMatrixSparseLockedByRowsAndByColumns &P,
		BzzVector &x,BzzVector *y);
	friend void SolveRight(BzzMatrixSparseLockedByRowsAndByColumns &P,BzzVector *y);

//	============================================================================
//	===========================	Iterative Systems Solution	===================
//	============================================================================
	int JacobiSolve(BzzVector &b,BzzVector *x);
	int GaussSeidelForwardSolve(BzzVector &b,BzzVector *x);
	int GaussSeidelBackwardSolve(BzzVector &b,BzzVector *x);
	int GaussSeidelSolve(BzzVector &b,BzzVector *x);
	int GaussSeidelForwardSolve(double omega,BzzVector &b,BzzVector *x);
	int GaussSeidelBackwardSolve(double omega,BzzVector &b,BzzVector *x);
	int GaussSeidelSolve(double omega,BzzVector &b,BzzVector *x);
//	============================================================================
//	=======================< Non-modifying functions >==========================
//	============================================================================
	virtual void ObjectBzzPrint(void);
};

*/
#endif // BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_BY_ROWS_HPP