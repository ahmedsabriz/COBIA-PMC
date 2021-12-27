// Serve per fare il prodotto fra una matrice sparsa ed un vettore
// nel caso in cui la matrice non debba essere ulteriormente modificata
// la matrice viene memorizzata in modo piu' efficiente

// BZZMATH: Release 7.0

//	=====================< BzzMatrixSparseLockedByColumns.hpp >==========================
//	* Class BzzMatrixSparseLockedByColumns for operations between matrices and vectors	*
// *																									*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixSparseLockedByColumns\MatrixSparseLockedByColumns.cpp		*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-2004	Date Written.

/////////// Release 6.0
//	11-2008	Added compact constructor.

//	============================================================================
//	******	BzzMatrixSparseLockedByColumns constructors:											*
//	* BzzMatrixSparseLockedByColumns A; // default													*
//	* BzzMatrixSparseLockedByColumns A = B; // from Sparse										*
//	* BzzMatrixSparseLockedByColumns A("SPARSE:DAT"); // from FILE								*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A.GetValue(i,j);																		*
//	* int count = BzzMatrixSparseLockedByColumns::ObjectCount();					*
//	* int countInScope = BzzMatrixSparseLockedByColumns::ObjectCountInScope();	*
//	* double *ptrVal = A.Scanning(&i,&j,&val);												*
//	* A.BeginScanning();																			*
//	****************************************************************************
//	***** Assignments:																			*
//	* ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked(&S,&A);						*
//	* BzzMatrixSparseLockedByColumns A;														*
//	* A = B; // from Sparse																		*
//	****************************************************************************
//	***** Functions:																				*
//	A.GetRowsNorm2(&nrom2R);																	*
//	A.GetColumnsNorm2(&nrom2C);																*
//	* ProductAndDifference(d,D,x,r); // r = d - D*x										*
//	****************************************************************************
//	***** Implemented operations :															*
//	* Product(A,x,&y);	// y = A*x;															*
//	* Product(A,B,&C);	// C = A*B;															*
//	* TProduct(A,v,&y);	// y = AT*x;														*
//	* ProductForSelectedRowsAndColumns(A,selR,selC,v,&y);								*
//	* ProductForSelectedColumns(A,selC,v,&y);												*
//	* TProductForSelectedColumns(A,selC,v,&y);											*
//	* Division(&A,2.);																			*
//	* Division(&A,v);																				*
//	****************************************************************************

#ifndef BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_COUMNS_HPP
#define BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_COUMNS_HPP

//	============================================================================
//	=========================< class BzzMatrixSparseLockedByColumns >====================
//	============================================================================
class BzzMatrixSparseLockedByColumns : public BzzBaseClass
{
	friend class BzzMatrixSparse;
	friend class BzzLinearProgramming;
	friend void ProductAndDifference(BzzVector& d,
		BzzMatrixSparseLockedByColumns& P, BzzVector& x,
		BzzVector* y);
	friend void TProduct(BzzMatrixSparseLockedByColumns& P, BzzVector& x, BzzVector* y);

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	double* matrix;
	int* rows;
	int* coeffInColumn;
	int* ptrColumn;
	int numRows, numColumns;
	int size;
	int whoAmI;
	//	BzzVectorInt ordSelectedColumns;
	//	char shadow;

		// initialisation constructors
	void Initialize(int mrows, int ncolumns, int nsize);

	// deinitialisation
	void Deinitialize(void);
	// prepares assignments
//	void PrepCopy(int rRows,int rColumns);

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default // BzzMatrixSparseLockedByColumns A;
	BzzMatrixSparseLockedByColumns(void);

	// copy-initializer // BzzMatrixSparseLockedByColumns A = B;
	BzzMatrixSparseLockedByColumns(BzzMatrixSparseLockedByColumns& E);

	// from BzzMatrixSparse
	BzzMatrixSparseLockedByColumns(BzzMatrixSparse& rval);

	// from FILE
	BzzMatrixSparseLockedByColumns(char* filematrix);

	// compact
	BzzMatrixSparseLockedByColumns(int nr, int nc, BzzVector* v,
		BzzVectorInt* r, BzzVectorInt* c);

	//	============================================================================
	//	****************************< destructor >**********************************
	//	============================================================================
	~BzzMatrixSparseLockedByColumns(void)
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
	void BzzPrintTransposeStructure(void);

	//	============================================================================
	//	**********************< assignment operators >******************************
	//	============================================================================

		// transforms a BzzMatrixSparse in BzzMatrixSparseLockedByColumns
		// BzzMatrixSparse gets destroyed
	friend void ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked
	(BzzMatrixSparse* lval, BzzMatrixSparseLockedByColumns* rval);
	friend void ReplaceBzzMatrixWithBzzMatrixSparseLocked
	(BzzMatrix* lval, BzzMatrixSparseLockedByColumns* rval);

	friend void Swap(BzzMatrixSparseLockedByColumns* lval,
		BzzMatrixSparseLockedByColumns* rval);

	BzzMatrixSparseLockedByColumns& operator =
		(BzzMatrixSparse& S);
	BzzMatrixSparseLockedByColumns& operator =
		(BzzMatrixSparseLockedByColumns& E);

	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByColumns& P,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByColumns* G);
	friend void CopyMatrixForSelectedRowsAndColumns(BzzMatrixSparseLockedByColumns& E,
		BzzMatrixSparseLockedByColumns& P,
		BzzVectorInt& selectedRows, BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByColumns* G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByColumns& P,
		BzzVectorInt& selectedColumns,
		BzzMatrixSparseLockedByColumns* G);
	friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByColumns& P,
		BzzVectorInt& selectedColumns,
		BzzMatrix* G);
	/*
		friend void CopyMatrixForSelectedColumns(BzzMatrixSparseLockedByColumns &P,
			BzzVectorInt &selectedColumns,
			BzzMatrixSparseStore &store,
			BzzMatrixSparse *G);
		friend void CopyMatrixForOrderedRowsAndColumns(BzzMatrixSparseLockedByColumns &P,
			BzzVectorInt &orderedRows,BzzVectorInt &orderedColumns,
			BzzMatrix *G);
		friend void CopyMatrixForOrderedRowsAndColumns(BzzMatrixSparseLockedByColumns &P,
			BzzVectorInt &orderedRows,BzzVectorInt &orderedColumns,
			BzzMatrixSparseStore &store,
			BzzMatrixSparse *G);
	*/

	//	============================================================================
	//	*****************************< Other functions >****************************
	//	============================================================================
	void GetRowsNorm2(BzzVector* norm2R);
	void GetColumnsNorm2(BzzVector* norm2C);

	void NormalizeColumns(BzzVector* norm2C);
	double TProductForSelectedColumn(int j, BzzVector& a);
	// a = a + cj*w
	void ProductForSelectedColumnAndSum(int j, double w, BzzVector* a);

	//	void GetRowsNorm2X(BzzVector &x,BzzVector *norm2X);
	//	void GetRow(int i,BzzVector *x);

	void GetColumn(int i, BzzVector* x);

	//	void GetRowForSelectedColumns(int i,BzzVectorInt &iv,BzzVector *x);
	//	void GetColumnForSelectedRows(int i,BzzVectorInt &iv,BzzVector *x);

	void GetNumEquationsForEachVariable
	(BzzVectorInt* numEquationsForEachVariable);
	void GetNumVariablesForEachEquation
	(BzzVectorInt* numVariablesForEachEquation);
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable);
	/*
		void GetEquationsContainingSingletons
			(BzzVectorInt *numVariablesForEachEquation,
			BzzVectorInt *numEquationsForEachVariable,
			BzzVectorInt *equationsContainingSingletons);

	//	void OrderingLinearProgramming(BzzVectorInt &fixVariables,
	//		BzzMatrixSparse *A,BzzMatrixSparse *B,
	//		BzzVectorInt *ordRows,BzzVectorInt *ordColumns);
		int OrderingLinearProgramming(BzzVectorInt *ordRows,
			BzzVectorInt *ordColumns);
		// questo credo si possa eliminare
		void OrderingLinearProgramming(BzzVectorInt &fixVariables,
			BzzVectorInt *ordRows,BzzVectorInt *ordColumns);
		// questo credo si possa eliminare
		int ChoiceLinearProgrammingVariables(BzzVectorInt &fixVariables,
			BzzVector &unorm2Columns,BzzVectorInt *varForFactorization);

		int ChoiceLinearProgrammingVariables(BzzVector &unorm2Columns,
			BzzVectorInt *varForFactorization);

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
	*/
	//	============================================================================
	//	*******************************< Product >**********************************
	//	============================================================================
	friend void Product(BzzMatrixSparseLockedByColumns& P,
		BzzVector& x, BzzVector* y);
	friend void Product(BzzMatrixSparseLockedByColumns& A, BzzMatrixSparseLockedByColumns& B,
		BzzMatrix* C);
	/*
		friend void TProduct(BzzMatrixSparseLockedByColumns &P,BzzVector &x,BzzVector *y);

		friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByColumns &D,
			BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
			BzzVector &x,BzzVector *y);
		friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLockedByColumns &E,
			BzzMatrixSparseLockedByColumns &D,
			BzzVectorInt &selectedRows,BzzVectorInt &selectedColumns,
			BzzVector &x,BzzVector *y);
	*/
	friend void ProductForSelectedColumn(BzzMatrixSparseLockedByColumns& D,
		int j, BzzVector& x, BzzVector* y);
	/*
	  friend void ProductForSelectedColumns(BzzMatrixSparseLockedByColumns &D,
			BzzVectorInt &selectedColumns,BzzVector &x,BzzVector *y);
		friend void TProductForSelectedColumns(BzzMatrixSparseLockedByColumns &D,
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
		friend void Division(BzzMatrixSparseLockedByColumns *P,double c);
		friend void Division(BzzMatrixSparseLockedByColumns *P,BzzVector &x);

	//	============================================================================
	//	=========================	Triangular Systems Solution	===================
	//	============================================================================
		friend void SolveLeft(BzzMatrixSparseLockedByColumns &P,BzzVector &x,BzzVector *y);
		friend void SolveLeft(BzzMatrixSparseLockedByColumns &P,BzzVector *y);
		friend void SolveRight(BzzMatrixSparseLockedByColumns &P,BzzVector &x,BzzVector *y);
		friend void SolveRight(BzzMatrixSparseLockedByColumns &P,BzzVector *y);
	*/
	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_MATRIX_DOUBLE_SPARSE_LOCKED_COUMNS_HPP