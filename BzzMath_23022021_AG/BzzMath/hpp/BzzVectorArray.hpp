// BZZMATH: Release 7.0

// ==========================< VectorArray.hpp >===============================
// * Class BzzVectorArray for *
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-2001	Date Written.

////////////////// Release 5.0
//	04-2005	Added Save and Load.
//	12-2006	Added Swap function.
//	12-2006	Added SwapRows function.
//	12-2006	Added Reorder function.
//	01-2007	Added MoveRowFromkToj function.
//	01-2007	Added DeleteElement function.
// 10-2007  Added GetTotalSize function.
// 10-2007  Added SetSubVector function.
// 10-2007  Added SetSize function.
// 07-2008  Added GetGrid function.
// 12-2009  Added InsertAndSumInSortedVectors function.
// 12-2009  Added InsertAndSubractInSortedVectors function.
// 12-2010  Added NormalizeRows function.

// ============================================================================
// ****** Constructors for BzzVectorArray:														*
// * BzzVectorArray v; // default																	*
// * BzzVectorArray v = x; // copy-initializer													*
// * BzzVectorArray v(n); // 			*
// ****************************************************************************
// ***** Access functions :																	*
// **
// ****************************************************************************
// ***** Assignment:																				*
// **
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// * v.BzzPrint();																				*
// * v.BzzMessage();																				*
// **
// ****************************************************************************
// ***** Implemented operations :															*
// **
// ****************************************************************************
// ***** Operators for tests:																	*
// **
// * if(v == w)																					*
// * if(v != w)																					*
// ****************************************************************************
// ***** Other functions:																		*
// **
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_VECTORARRAY_DOUBLE_HPP
#define BZZ_VECTORARRAY_DOUBLE_HPP

// ============================================================================
// ============================< class BzzVectorArray >=============================
// ============================================================================

class BzzVectorArray : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzMatrixSparseLockedByRows;

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int numVectors;
	int numColumns;

	// initialise constructors
	int	whoAmI;

public:
	BzzVector* v;

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzVectorArray v;
	BzzVectorArray(void);

	// copy-initializer
	BzzVectorArray(BzzVectorArray& rval)
	{
		BzzError("Non Implemented");
	}

	// other constructor
	BzzVectorArray(int nv);
	void operator()(int nv);

	BzzVectorArray(int nv, int nc);
	void SetDimensions(int nv, int nc);

	BzzVectorArray(int nc, BzzVectorInt& nr);
	void operator()(int nc, BzzVectorInt& nr);

	BzzVectorArray(BzzVectorInt& nr);
	void operator()(BzzVectorInt& nr);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzVectorArray(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	int GetNumVectors(void) { return numVectors; }
	int GetSize(int elem);
	void GetSize(BzzVectorInt* elements);
	int GetTotalSize(void);
	void SetSize(int elem, int dim);

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzVectorArray& operator =
		(BzzVectorArray& rval);
	void operator()(int i, BzzVector& b);
	void operator()(int i, BzzVector* b);
	void operator()(int i, int n, BzzVector& b);
	void SetVector(int i, BzzVector& b);
	void SetVector(int i, BzzVector* b);
	void SetSubVector(int i, BzzVector& b, int dim);
	BzzVector GetVector(int i);
	void GetVector(int i, BzzVector* b);
	friend void Swap(BzzVectorArray* V, BzzVectorArray* W);
	friend void InsertAndSumInSortedVectors(BzzVectorInt* iv, BzzVector* dv,
		BzzVectorInt& iw, BzzVector& dw);
	friend void InsertAndSubractInSortedVectors(BzzVectorInt* iv, BzzVector* dv,
		BzzVectorInt& iw, BzzVector& dw);

	// ============================================================================
	// *************************< access functions >*******************************
	// ============================================================================
		// assigns and receives values with control
	double& operator () (int row, int col);

	// assigns and receives values without control
	double* operator [] (int r)
	{
		return v[r].GetHandle() - 1;
	}

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// ======================< Modifying Functions >===============================
	// ============================================================================
	friend void Delete(BzzVectorArray* result); // eliminates BzzVectorArray
	void DeleteVector(int i);
	void DeleteElement(int i, int j);
	void SwapVector(int i, BzzVector* v);
	void SwapRows(int i, int j);
	void Reorder(BzzVectorInt& iS);
	void MoveRowFromkToj(int k, int j);
	void NormalizeRows(BzzVector* nr);
	friend void FindVariableSingletonsInEquations(int nr, int nc, BzzVectorInt& r, BzzVectorInt& c,
		BzzVector& v, BzzVectorIntArray* rowV, BzzVectorArray* valRowV,
		BzzVectorIntArray* colV, BzzVectorArray* valColV,
		BzzVector* e, BzzVector* x);
	// recovery from Save
//	friend void Load
//			(BzzVectorArray *result,char *file***); // formatted
//	friend void Load
//			(BzzVectorArray *result,char,char *file***);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
};

// ============================================================================
// ======================< class BzzMatrixLeftSparse >=========================
// ============================================================================
class BzzMatrixLeftSparse : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int whoAmI;

	int	numRows,
		numColumns,
		numLowerElements;
	int diagonalUnitary;
	BzzVector d;
	BzzVectorInt r, c;
	BzzVector l;
	BzzVectorIntArray rI;
	BzzVectorArray rV;
	BzzVectorIntArray cI;
	BzzVectorArray cV;
	BzzVectorInt auxi; // serve per i Build
	void BuildRowsFromRCV(void);
	void BuildColumnsFromRCV(void);
	void BuildRCVFromRows(void);
	void BuildRCVFromColumns(void);

public:
	void SetDiagonalUnitary(void)
	{
		diagonalUnitary = 1;
	}
	BzzVectorInt	numVariablesInEachEquation,
		numEquationsForEachVariable;
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft L
	BzzMatrixLeftSparse(void);
	void operator()(BzzVector* dd, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* ll);
	void operator()(int nr, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* ll);

	// copy-initializer // BzzMatrixLeft L = left;
//	BzzMatrixLeftSparse(BzzMatrixLeftSparse &rval);

	// sizes and initialises at 0
//	BzzMatrixLeftSparse(int rows,int columns);

		// sizes and initialises at 0
//	BzzMatrixLeftSparse(int rows,int columns,BzzVector *dd,
//		BzzVectorInt *rr,BzzVectorInt *cc,BzzVector *ll);

	// from formatted File // BzzMatrixLeft L("LEFT.DAT");
//	BzzMatrixLeft(char *filematrix);

	// from binary File // BzzMatrixLeft L('*',"LEFT.BIN");
//	BzzMatrixLeft(char,char *filematrix);

// ============================================================================
// ************************< Modifying access functions >**********************
// ============================================================================
	void RemoveLeftRow(int row);
	void SetRow(int row, BzzVectorInt* cc, BzzVector* vv, double dd);
	void AppendRow(BzzVectorInt* cc, BzzVector* vv, double dd);
	void RemoveLeftColumn(int col);
	void SetColumn(int col, double dd, BzzVectorInt* rr, BzzVector* vv);
	void SetDiagonal(BzzVector* dd);
	void SetDiagonalElement(int id, double dd);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixLeftSparse(void) {};

	// ============================================================================
	// ********************< Non-modifying access functions >**********************
	// ============================================================================
		// row number
	int Rows(void) const
	{
		return numRows;
	}

	// column number
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================

	// ============================================================================
	// ****************************< Product >*************************************
	// ============================================================================
	friend void Product(BzzMatrixLeftSparse& L, BzzVector& x, BzzVector* y);

	// ============================================================================
	// =======================< Functions for linear algebra >=====================
	// ============================================================================
	friend void Solve(BzzMatrixLeftSparse& L, BzzVector& b, BzzVector* x);
};

// ============================================================================
// ======================< class BzzMatrixLeftSparseForAttic >=========================
// ============================================================================
class BzzMatrixLeftSparseForAttic : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int whoAmI;
	int	numRows,
		numColumns;
public:
	BzzVectorIntArray rI;
	BzzVectorArray rV;
	BzzVectorInt	numVariablesInEachEquation;

	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft L
	BzzMatrixLeftSparseForAttic(void);
	BzzMatrixLeftSparseForAttic(int numr, int numc);
	void operator()(int numr, int numc);

	// ============================================================================
	// ************************< Modifying access functions >**********************
	// ============================================================================
	void RemoveRow(int row);
	void SetRow(int row, BzzVectorInt* cc, BzzVector* vv);
	void SetColumn(int column, BzzVectorInt& rr, BzzVector& vv);
	void AppendRow(BzzVectorInt* cc, BzzVector* vv);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixLeftSparseForAttic(void) {};

	// ============================================================================
	// ********************< Non-modifying access functions >**********************
	// ============================================================================
		// row number
	int Rows(void) const
	{
		return numRows;
	}

	// column number
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================

	// ============================================================================
	// ****************************< Product >*************************************
	// ============================================================================
	//friend void Product(BzzMatrixLeftSparseForAttic &L,BzzVector &x,BzzVector *y);

	// ============================================================================
	// =======================< Functions for linear algebra >=====================
	// ============================================================================
	//friend void Solve(BzzMatrixLeftSparseForAttic &L,BzzVector &b,BzzVector *x);
};

//==============================================================================
//==============================================================================
//==============================================================================

// ============================================================================
// ======================< class BzzMatrixRightSparse >=========================
// ============================================================================
class BzzMatrixRightSparse : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int whoAmI;

	int	numRows,
		numColumns,
		numUpperElements;
	int diagonalUnitary;
	BzzVector d;
	BzzVectorInt r, c;
	BzzVector u;
	BzzVectorIntArray rI;
	BzzVectorArray rV;
	BzzVectorIntArray cI;
	BzzVectorArray cV;
	BzzVectorInt auxi; // serve per i Build
	void BuildRowsFromRCV(void);
	void BuildColumnsFromRCV(void);
	void BuildRCVFromRows(void);
	void BuildRCVFromColumns(void);

public:
	void SetDiagonalUnitary(void)
	{
		diagonalUnitary = 1;
	}
	BzzVectorInt	numVariablesInEachEquation,
		numEquationsForEachVariable;
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft R
	BzzMatrixRightSparse(void);
	void operator()(BzzVector* dd, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* uu);
	void operator()(int nr, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* uu);

	// ============================================================================
	// ************************< Modifying access functions >**********************
	// ============================================================================
	void RemoveRightRow(int row);
	void SetRow(int row, BzzVectorInt* cc, BzzVector* vv, double dd);
	void AppendColumn(BzzVectorInt* rr, BzzVector* vv, double dd);
	void RemoveRightColumn(int col);
	void SetColumn(int col, double dd, BzzVectorInt* rr, BzzVector* vv);
	void SetDiagonal(BzzVector* dd);
	void SetDiagonalElement(int id, double dd);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixRightSparse(void) {};

	// ============================================================================
	// ********************< Non-modifying access functions >**********************
	// ============================================================================
		// row number
	int Rows(void) const
	{
		return numRows;
	}

	// column number
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================

	// ============================================================================
	// ****************************< Product >*************************************
	// ============================================================================
	friend void Product(BzzMatrixRightSparse& R, BzzVector& x, BzzVector* y);

	// ============================================================================
	// =======================< Functions for linear algebra >=====================
	// ============================================================================
	friend void Solve(BzzMatrixRightSparse& R, BzzVector& b, BzzVector* x);
};

void GetGrid(BzzVectorArray& x, BzzMatrix* X);

#endif // BZZ_VECTORARRAY_DOUBLE_HPP