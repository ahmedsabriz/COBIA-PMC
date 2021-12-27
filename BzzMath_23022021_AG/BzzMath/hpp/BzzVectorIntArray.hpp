// BZZMATH: Release 7.0

// =======================< VectorIntArray.hpp >===============================

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

////////////////// Release 6.0
// 06-2009  Added constructors from FILES.
// 08-2010  Added BzzCuthillMcKleeOrdering functions.
// 08-2010  Added BzzReverseCuthillMcKleeOrdering functions.
// 08-2010  Added GetLowerEnvelope function.
// 08-2010  Added GetUpperEnvelope function.
// 11-2010  Added BzzPrintStructure function.
//	03-2012	Added GetNonDuplicateVector function.

#ifndef BZZ_VECTOR_INT_ARRAY_HPP
#define BZZ_VECTOR_INT_ARRAY_HPP

// ============================================================================
// ====================< class BzzVectorIntArray >=============================
// ============================================================================

class BzzVectorIntArray : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzMatrixSparseLocked;
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int numVectors;
	int numColumns;
	//	BzzVectorInt *v;
	int	whoAmI;
public:
	BzzVectorInt* v;

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzVectorIntArray v;
	BzzVectorIntArray(void);

	// copy-initializer
	BzzVectorIntArray(BzzVectorIntArray& rval)
	{
		BzzError("Non Implemented");
	}

	// other constructor
	BzzVectorIntArray(int nv);
	void operator()(int nv);

	BzzVectorIntArray(int nr, int nc);
	void SetDimensions(int nr, int nc);

	BzzVectorIntArray(int nc, BzzVectorInt& nr);
	void operator()(int nc, BzzVectorInt& nr);

	BzzVectorIntArray(BzzVectorInt& nr);
	void operator()(BzzVectorInt& nr);

	BzzVectorIntArray(char* filematrix);
	void operator()(char* filematrix);

	BzzVectorIntArray(char* filematrix, double val);
	void operator()(char* filematrix, double val);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzVectorIntArray(void);
	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	int GetNumVectors(void) { return numVectors; }
	int Rows(void) { return numVectors; }
	int Columns(void) { return numColumns; }
	int GetSize(int elem);
	void GetSize(BzzVectorInt* elements);
	int GetTotalSize(void);
	void SetSize(int elem, int dim);

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzVectorIntArray& operator =
		(BzzVectorIntArray& rval);

	//	void operator()(int i,BzzVectorInt &b);
	void operator()(int i, BzzVectorInt* b);
	void operator()(int i, int n, BzzVectorInt& b);

	void SetVector(int i, BzzVectorInt& b);
	void SetVector(int i, BzzVectorInt* b);
	void SetSubVector(int i, BzzVectorInt& b, int dim);

	BzzVectorInt GetVector(int i);
	void GetVector(int i, BzzVectorInt* b);

	friend void Swap(BzzVectorIntArray* V, BzzVectorIntArray* W);
	friend void BzzCuthillMcKleeOrdering(int start, BzzVectorIntArray& V, BzzVectorInt* v);
	friend void BzzCuthillMcKleeOrdering(BzzVectorIntArray& V, BzzVectorInt* v);
	friend void BzzCuthillMcKleeOrdering(int start, BzzVectorIntArray* V,
		BzzVectorArray* D, BzzVectorInt* v);
	friend void BzzCuthillMcKleeOrdering(BzzVectorIntArray* V,
		BzzVectorArray* D, BzzVectorInt* v);
	friend void BzzReverseCuthillMcKleeOrdering(int start, BzzVectorIntArray& V, BzzVectorInt* v);
	friend void BzzReverseCuthillMcKleeOrdering(int start, BzzVectorIntArray* V,
		BzzVectorArray* D, BzzVectorInt* v);

	// ============================================================================
	// *************************< access functions >*******************************
	// ============================================================================
		// assigns and receives vector values with control
	int& operator () (int row, int col);

	// assigns and receives vector values without control
	int* operator [] (int r)
	{
		return v[r].GetHandle() - 1;
	}

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void);
	void	Save(char* filematrix); // formatted
	void Transpose(BzzVectorIntArray* T);

	// ============================================================================
	// ======================< Modifying Functions >===============================
	// ============================================================================
	friend void Delete(BzzVectorIntArray* result); // eliminates BzzVectorIntArray
	void DeleteVector(int i);
	void DeleteElement(int i, int j);
	void SwapVector(int i, BzzVectorInt* v);
	void SwapRows(int i, int j);
	void Reorder(BzzVectorInt& iS);
	void MoveRowFromkToj(int k, int j);

	// recovery from Save
//	friend void Load
//			(BzzVectorIntArray *result,char *file***); // formatted
//	friend void Load
//			(BzzVectorIntArray *result,char,char *file***);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
	int Sum(void);
	void GetMatrixForOdebDerivatives(BzzMatrixInt* A);
	void GetLowerEnvelope(BzzVectorInt* lowerEnvelpe);
	void GetUpperEnvelope(BzzVectorInt* upperEnvelope);
	void GetNonDuplicateVector(BzzVectorInt& iv, BzzVectorInt* nm);
};

#endif // BZZ_VECTOR_INT_ARRAY_HPP

/*

// BZZMATH: Release 7.0

// =========================< VECTORINTARRAY.HPP >=============================
// * Class BzzVectorIntArray 																	*
// * Examples: c:\bzzmath\examples\exvectorintarray.cpp							 	*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2003	Date Written.

#ifndef BZZ_VECTOR_INT_ARRAY_HPP
#define BZZ_VECTOR_INT_ARRAY_HPP

// ============================================================================
// =========================< class BzzVectorIntArray >========================
// ============================================================================

class BzzVectorIntArray : public BzzBaseClass
	{
	friend class BzzMatrixIntArray;
private:
	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int **vector;
	int numRows,size;
	int whoAmI;
//	char shadow;

	BzzVectorInt numElementsVector;
	// initialise constructors
//	void Initialize(int nc);

public:

// ============================================================================
// ***************************< constructors >*********************************
// ============================================================================
	// default constructor BzzVectorIntArray v;
	BzzVectorIntArray(void);

	// copy-initializer
	BzzVectorIntArray(BzzVectorIntArray &rval)
		{BzzError("Copy-initializer not yet available");}

	// sizes
	BzzVectorIntArray(int nR,BzzVectorInt &iv);
	void operator()(int nR,BzzVectorInt &iv);
	void SetVector(int nR,BzzVectorInt &iv);
	void SetVector(int nR,BzzVectorInt *iv);

// ============================================================================
// *****************************< destructor >*********************************
// ============================================================================
	~BzzVectorIntArray(void);
	friend void Delete(BzzVectorIntArray *v);

// ============================================================================
// *******************< Non-modifying access functions >***********************
// ============================================================================
	int Size(void) const {return size - 1;} // dimensions

	int WhoAmI(void) const {return whoAmI;}
	static int ObjectCount(void){return count;}
	static int ObjectCountInScope(void){return countInScope;}
	int GetNumVectors(void){return numRows;}

// ============================================================================
// **********************< Modifying access functions >************************
// ============================================================================
	int *operator [] (int r)
		{return vector[r];}
	int &operator()(int row,int column);

// ============================================================================
// *************************< Assignment operators >***************************
// ============================================================================

	void SetVectorRow(int i,BzzVectorInt &iv);
	void operator()(BzzVectorInt &iTot); // initialize all elements
	BzzVectorIntArray &operator =
		 (BzzVectorIntArray &rval);

// ============================================================================
// ====================< Non-modifying functions >=============================
// ============================================================================

// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	int Sum(void);
// ============================================================================
// ======================< Modifying Functions >===============================
// ============================================================================

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
	int *GetHandle(void){return vector[0] + 1;}
};

#endif // BZZ_VECTOR_INT_ARRAY_HPP
*/