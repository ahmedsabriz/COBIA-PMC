// BZZMATH: Release 7.0

// =============================< VECINT.HPP >=================================
// * Class BzzVectorInt																			*
//	*					Metodi numerici e Software in C++ (Capitolo 3)					*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: c:\bzzmath\examples\BzzMathBasic\LinearAlgebra\						*
// *           VectorInt\VectorInt.cpp													 	*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-1995	Date Written.
//	12-1996	Added DeleteOneElementEveryN.
//	01-1997	Added Stretch.
//	04-1997	Added InsertSortedVectorInSortedVector.
//	04-1997	Added InsertNonMatchingSortedVectorInSortedVector.
//	04-1997	Added InsertElementInSortedVector.
//	04-1997	Added InsertNonMatchingElementInSortedVector.
//	04-1997	Added ExistingElementInVector.
//	04-1997	Added ExistingElementInSortedVector.
//	04-1997	Added OperatorAndForSortedVectors.
//	04-1997	Added DeleteMatchingElementInSortedVector.
//	10-1998	Added IsVectorSorted.
//	12-1999	Added AutoSwap.
//////// Release 4.0
//	03-2000	Added DeleteLastNElements.
// 01-2001	Added GetHandle function for mixed language: C++ / FORTRAN.
//	05-2001	Added new version for GetBzzVector.
//	09-2001	Added new version for Sort.
//	04-2002	Added Reorder.
//	04-2002	LocateInFirstNSortedElements.
//	04-2002	ExistingElementInFirstNSortedElements.
//	04-2002	Added ExistingOrLocateElementInSortedVector.
//	04-2002	Added MoveElementFromkToj.
//	12-2002	Added Sum.
//	04-2003	Added InsertElementInFirstNSortedElements.

////////////////// Release 5.0
//	10-2003	Added DeleteFirstNElements.
//	10-2003	Added DeleteFirstAndLastNElements.
//	10-2003	Added SetToZero function.
//	11-2004	Added AppendNonMatchingElement function.
//	11-2004	Added AppendNonMatchingElement function.
//	11-2004	GetSumElements function.
//	02-2006	Added constructor from file.
//	02-2006	Added ExistingDuplicationsInVector function.
//	02-2006	Added ExistingNonZeroDuplicationsInVector function.
//	02-2006	Added CopyDataFromVector function.
//	02-2006	Added UseSubVectorAsVector function.
//	12-2006	Added MinMax function.
//	06-2007	Added CompleteTheVector function.
//	06-2008	Added SwapMaxMinOrder function.

////////////////// Release 6.0
//	11-2008	Added DeleteMin, DeleteMax and DelteMinMax functions.
//	01-2009	Added ExistingDuplicationsInSortedVector function.
//	01-2009	Added ExistingNonZeroDuplicationsInSortedVector function.
//	01-2009	Added DeleteDuplicationsInVector function.
//	01-2009	Added DeleteNonZeroDuplicationsInVector function.
//	01-2009	Added DeleteDuplicationsInSortedVector function.
//	01-2009	Added DeleteNonZeroDuplicationsInSortedVector function.
//	06-2009	Added CrescentInitialization and DecrescentInitialization functions.
//	12-2009	Modified DeleteElement function.
//	12-2009	Added InsertNonMatchingElementInFirstNSortedElements.
//	12-2009	Added InsertAndSumNonMatchingElementInFirstNSortedElements function.
//	12-2009	Added InsertAndSubtractNonMatchingElementInFirstNSortedElements function.
//	01-2010	Added MinMaxMean function.
//	01-2010	Added ReplaceNElements function.
// 04-2011	Added MatchingElementsNumberInTwoSortedVectors function.
// 04-2011	Added DeleteMatchingVectorElementsInSortedVector function.
// 04-2011	Added DeleteMatchingFirstNVectorElementsInSortedVector function.
// 04-2012	Added InsertNonMatchingAndSumSortedVectorInSortedVector function.
// 04-2012	Added InsertNonMatchingAndSubractSortedVectorInSortedVector function.
//	04-2012	Added DeleteMatchingElementInFirstNSortedElements function.
// 12-2012	Added GetDuplicationsInVector function

////////////////// Release 7.0
//	06-2013	Added SwappedReordering function.
//	06-2013	Added ReorderSubVector function.
//	06-2013	Added SwappedSubVectorReordering function.

// ============================================================================
// ****** Constructors for BzzVectorInt:													*
// * BzzVectorInt iv; // default																*
// * BzzVectorInt iv = jx; // copy-initializer											*
// * BzzVectorInt iv(n);// sizes BzzVectorInt iv and set its coifficients = 0	*
// * int i[5] = {1,2,3,4,5};																	*
// * BzzVectorInt iv(5,i); // from array													*
// * BzzVector v(5,1,2,3,4,5); // vector 1,2,3,4,5										*
// ****************************************************************************
// ***** Access functions :																	*
// * int i = iv.Size(); // BzzVectorInt v size											*
// * int who = iv.WhoAmI(); // BzzVectorInt count										*
// * int j = iv(i);	// with control														*
// * iv(i) = 4;			// with control													*
// * int j = iv[i];	// inline without control											*
// * iv[i] = 7;			// inline without control										*
// * int j = v.GetValue(i); // with control												*
// * int count = BzzVectorInt::ObjectCount();											*
// * int countInScope = BzzVectorInt::ObjectCountInScope();							*
// ****************************************************************************
// ***** Assignment:																				*
// * iv = jw; // Assignment																	*
// * iv = 3; // Assignment from int															*
// * iv = jw(ix); // Assignment from a selection of jw								*
// * w(5,v); // w subvector of v	dimensioned 5											*
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// * iv.BzzPrint();																				*
// * iv.BzzMessage();																			*
// * iv.BzzPrint("The value of iv(3) is %d",iv[3]);									*
// * iv.BzzMessage("The value of iv(3) is %e",iv[3]);									*
// ****************************************************************************
// ***** Operators for tests:																	*
// * if(iv == jw)																					*
// * if(iv != jw)																					*
// * if(v == c)																					*
// * if(v > c)																						*
// * if(v >= c)																					*
// * if(v < c)																						*
// * if(v <= c)																					*
// ****************************************************************************
// ***** Other functions:																		*
// * int j;																							*
// * int *vector = iv.GetHandle();															*
// * js = v.Sum();	// Sum of elements													*
// * j = v.Max();																					*
// * j = v.Max(&imax);																			*
// * j = v.MaxAbs();																				*
// * j = v.MaxAbs(&imax);																		*
// * j = v.Min();																					*
// * j = v.Min(&imin);																			*
// * v.MinMax(&iMin,&min,&iMax,&max);														*
// *	v.MinMaxMean(&iMin,&min,&iMax,&max,&mean);										*
// * j = v.MinAbs();																				*
// * j = v.MinAbs(&imin);																		*
// * j = v.MinExludedK(k,&imin);	//k BzzVectorInt										*
// * BzzVectorInt w = GetBzzVector(nc,i,v); // w subvector of v					*
// * w.GetBzzVector(nc,i,v); // w subvector of v										*
// * w.UseSubVectorAsVector(nc,i,v); // w subvector of v with the same data!!	*
// * Delete(&iv); // eliminate a BzzVectorInt											*
// * v.SetToZeroDelete(); // set all elements to zero									*
// * ChangeDimensions(newdim,&iv);															*
// * Swap(&ix,&iy);																				*
// * AutoSwap(&ix);																				*
// * Reverse(&iv);																				*
// * Sort(&iv);																					*
// * SwapMaxMinOrder(&iv);																		*
// * Reorder(&v,iS);																				*
// * j = LocateInSortedVector(v,f);															*
// * j = v.LocateInSortedVector(f);															*
// * j = v.LocateInFirstNSortedElements(n,f);											*
// * j = v.IsVectorSorted();// return 0 non sorted; 1 ascending; 2 descending	*
// * v.FirstInLastOut(j);																		*
// * v.LastInFirstOut(j);																		*
// * v.Append(j);																					*
// * v.Append(vj);																				*
// * v.AppendNonMatchingElement(j);															*
// * v.Insert(3,j);																				*
// * v.Insert(3,vj);																				*
// * v.InsertSortedVectorInSortedVector(w);												*
// * v.InsertNonMatchingSortedVectorInSortedVector(w);								*
// * v.InsertElementInSortedVector(4);														*
// * v.InsertElementInFirstNSortedElements(10,5);										*
// * v.MoveElementFromkToj(k,j);																*
// * char c = v.InsertNonMatchingElementInSortedVector(4);							*
// * int i = v.ExistingElementInVector(5);												*
// * int i = v.ExistingElementInSortedVector(5);										*
// * int i = v.ExistingElementInFirstNSortedElements(10,5);							*
// * int i = v.ExistingOrLocateElementInSortedVector(5);								*
// * char c = v.OperatorAndForSortedVectors(v,w);										*
// * char c = v.OperatorAndForSortedVectors(v,w,5);									*
// * char c = v.OperatorAndForSortedVectors(v,w,&x);									*
// * char c = v.OperatorAndForSortedVectors(v,w,&x,5);								*
// * j = v.DeleteElement(j);																	*
// * v.DeleteLastNElements(n);																*
// * v.DeleteFirstNElements(n);																*
// * v.DeleteFirstAndLastNElements(n);														*
// * v.DeleteElements(iv);																		*
// * v.DeleteOneElementEveryN(n);															*
// *	v.DeleteMin(); v.DeleteMax(); v.DelteMinMax();									*
// * char c = v.DeleteMatchingElementInSortedVector(f);								*
// * v.CrescentInitialization();																*
// * v.DecrescentInitialization();															*
// * v.Stretch(n);																				*
// * v.SwapElements(i,j);																		*
// * v.Save("Vecint.DAT");																		*
// * v.Save('*',"Vecint.BIN");																*
// * Load(&v,"Vecint.DAT");																	*
// * Load(&v,'*',"Vecint.BIN");																*
// ****************************************************************************

#ifndef BZZ_VECTOR_INT_HPP
#define BZZ_VECTOR_INT_HPP

// ============================================================================
// ============================< class BzzVectorInt >==========================
// ============================================================================

class BzzVectorInt : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrixSparse;
	friend class BzzMatrixSparseLokedByRows;
	friend class BzzMatrixSparseLokedByColumns;
	friend class BzzMatrixInt;
	friend class BzzSave;
	friend class BzzLoad;
	friend void Load(BzzMatrixSparse* A, char, char* filematrix);

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int* vector;
	int dimensions;
	int whoAmI;
	char shadow;
	char matrixAsVector;
	char subVectorAsVector;

	// initialise constructors
	void Initialize(int nc);

	// private constructor BzzVectorInt('*',nc)
	BzzVectorInt(char, int nc);

public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzVectorInt v;
	BzzVectorInt(void);

	// copy-initializer

#if BZZ_COMPILER == 101
	BzzVectorInt(const BzzVectorInt& rval);
#else
	BzzVectorInt(BzzVectorInt& rval);
#endif
	// sizes and initialises at 0
	BzzVectorInt(int nc);

	// sizes and initialises
	BzzVectorInt(int nc, int v1, ...);

	// data pointer constructor
	// int i[] = {1,2,3,4};
	// BzzVectorInt iv(4,i);
	BzzVectorInt(int nc, int* initvalues);

	// file v(FILE) constructor;
	BzzVectorInt(char* filevector);

	// unformatted v('*',FILE) constructor saved with Save
	BzzVectorInt(char, char* filevector);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzVectorInt(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int Size(void) const { return dimensions; } // dimensions

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// receives the value of the vector with control
	int GetValue(int i) const;

	// subvector
	friend BzzVectorInt GetBzzVector
	(int nc, int ielem, const BzzVectorInt& rval);

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// subvector
	void GetBzzVector
	(int nc, int ielem, const BzzVectorInt& rval);

	void UseSubVectorAsVector(int nc, int ielem, const BzzVectorInt& rval);

	// assigns and receives vector values with control
	int& operator () (int i);

	// assigns and receives values without control
	int& operator []
	(int i) {
		return vector[i];
	}

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzVectorInt& operator =
		(const BzzVectorInt& rval);

	void operator =
		(int c);

	//w.(n,v);
	void operator() (int n, BzzVectorInt& v);
	void operator() (int n, int start, BzzVectorInt& v);
	void CrescentInitialization(void);
	void DecrescentInitialization(void);

	// ============================================================================
	// ****************************< operator (iv) >*******************************
	// ============================================================================

	BzzVectorInt operator () (const BzzVectorInt& iv);

	// ============================================================================
	// *************************< operators for tests >****************************
	// ============================================================================
	friend char operator ==
		(const BzzVectorInt& lval, const BzzVectorInt& rval);

	friend char operator !=
		(const BzzVectorInt& lval, const BzzVectorInt& rval);

	char operator == (int c);
	char operator > (int c);
	char operator >= (int c);
	char operator < (int c);
	char operator <= (int c);

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);

	// ******************************< Max and Min >*******************************
	int Max(int* imax = 0); //to have the Max position (im)
	int MaxAbs(int* imax = 0);

	int Min(int* imin = 0);
	int MinAbs(int* imin = 0);

	void MinMax(int* iMin, int* min, int* iMax, int* max);
	void MinMaxMean(int* iMin, int* min, int* iMax, int* max, double* mean);
	int GetCleverMeanDeviationMaxOutliers(int* mean, int* deviation,
		BzzVectorInt* iout, BzzVectorInt* out);

	int MinExludedK(BzzVectorInt& k, int* imin = 0);

	void Save(char* filevector); // formatted
	void Save(char, char* filevector);// binary

// ============================================================================
// ======================< Modifying Functions >===============================
// ============================================================================
	friend void Delete(BzzVectorInt* result); // eliminates BzzVectorInt
	void SetToZero(void)
	{
		memset(vector, 0, (dimensions + 1) * sizeof(int));
	}
	void SetToZero(int start, int end);
	friend void ChangeDimensions(int dim, BzzVectorInt* result,
		char);

	friend void Reverse(BzzVectorInt* result);	//inverts the vector

	friend void Sort(BzzVectorInt* result); // orders the vector
	friend void Sort(BzzVector* result, BzzVectorInt* iS); // orders the vector
	friend void Sort(BzzVectorInt* result, BzzVectorInt* iS); // orders the vector

	// return vector v ordered with iS
	friend void Reorder(BzzVectorInt* v, BzzVectorInt& iS);
	// return vector v ordered with iS
	friend void SwappedReording(BzzVectorInt* v, BzzVectorInt& iS);
	// return a portion of vector v starting from k ordered with iS
	friend void ReorderSubVector(int k, BzzVectorInt* v, BzzVectorInt& iS);
	// return a portion of vector v starting from k ordered with iS
	friend void SwappedSubVectorReordering(int k, BzzVectorInt* v, BzzVectorInt& iS);

	// swaps the contents of two vectors
	friend void Swap(BzzVectorInt* lval, BzzVectorInt* rval);

	// swaps the contents of a vector with its order
	friend void AutoSwap(BzzVectorInt* ix);

	//reverse the order Max min
	friend void SwapMaxMinOrder(BzzVectorInt* ix);

	// recovery from Save
	friend void Load
	(BzzVectorInt* result, char* filevector); // formatted
	friend void Load
	(BzzVectorInt* result, char, char* filevector);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================

	int* GetHandle(void) { return vector + 1; }
	int Sum(void);
	friend int LocateInSortedVector(BzzVectorInt& v, int f);
	int LocateInSortedVector(int f);
	int LocateInFirstNSortedElements(int n, int f);
	void FirstInLastOut(int f);
	void LastInFirstOut(int f);
	void Append(int f);
	void Append(BzzVectorInt& v);
	int AppendNonMatchingElement(int f);
	void Insert(int i, int f);
	void Insert(int i, BzzVectorInt& v);
	void InsertSortedVectorInSortedVector(BzzVectorInt& vj);
	void InsertNonMatchingSortedVectorInSortedVector(BzzVectorInt& vj);
	friend void InsertNonMatchingAndSumSortedVectorInSortedVector(BzzVectorInt& j1, BzzVector& v1,
		BzzVectorInt& j2, BzzVector& v2, BzzVectorInt* jf, BzzVector* vf);
	friend void InsertNonMatchingAndSubtractSortedVectorInSortedVector(BzzVectorInt& j1, BzzVector& v1,
		BzzVectorInt& j2, BzzVector& v2, BzzVectorInt* jf, BzzVector* vf);
	friend int NonMatchingElementsInTwoSortedVectors(BzzVectorInt& v1, BzzVectorInt& v2);
	int InsertElementInSortedVector(int f);
	int InsertElementInFirstNSortedElements(int nE, int f);
	char InsertNonMatchingElementInSortedVector(int f);
	int InsertNonMatchingElementInFirstNSortedElements(int nE, int f);
	int InsertAndSumNonMatchingElementInFirstNSortedElements
	(BzzVector* dv, int n, int iw, double dw);
	int InsertAndSubtractNonMatchingElementInFirstNSortedElements
	(BzzVector* dv, int n, int iw, double dw);
	void ReplaceNElements(BzzVectorInt& b);
	void ReplaceNElements(int j, BzzVectorInt& b);
	void CompleteTheVector(int m);
	void MoveElementFromkToj(int k, int j);
	int ExistingElementInVector(int f);
	int ExistingElementInSortedVector(int f);
	int ExistingElementInFirstNSortedElements(int n, int f);
	int ExistingOrLocateElementInSortedVector(int f);
	friend char OperatorAndForSortedVectors(BzzVectorInt& v, BzzVectorInt& w,
		int);
	friend char OperatorAndForSortedVectors(BzzVectorInt& v, BzzVectorInt& w,
		BzzVectorInt* x, int);
	int DeleteElement(int i);
	void DeleteLastNElements(int n);
	void DeleteFirstNElements(int n);
	void DeleteFirstAndLastNElements(int n);
	void DeleteElements(BzzVectorInt& v);
	void DeleteOneElementEveryN(int n);
	int DeleteMatchingElementInSortedVector(int f);
	int DeleteMatchingVectorElementsInSortedVector(BzzVectorInt& v);
	int DeleteMatchingVectorElementsInSortedVector(BzzVectorInt& v,
		BzzVector* dv);
	int DeleteMatchingElementInFirstNSortedElements(int n, int f);
	int DeleteMatchingElementInFirstNSortedElements(int n, int f, BzzVector* dv);
	int DeleteMatchingVectorElementsInFirstNSortedElements(int n, BzzVectorInt& v);
	int DeleteMatchingVectorElementsInFirstNSortedElements(int n, BzzVectorInt& v,
		BzzVector* dv);

	int DeleteMatchingFirstNVectorElementsInSortedVector(int N, BzzVectorInt& v);
	void DeleteMax(void);
	void DeleteMin(void);
	void DeleteMinMax(void);
	void Stretch(int i);
	void SwapElements(int i, int j);
	int IsVectorSorted(void); // return 0 non sorted; 1 ascending; 2 descending
	int GetSumElements(void);
	int ExistingDuplicationsInVector(void);
	int ExistingNonZeroDuplicationsInVector(void);
	int ExistingDuplicationsInSortedVector(void);
	int ExistingNonZeroDuplicationsInSortedVector(void);
	void DeleteDuplicationsInVector(void);
	void GetDuplicationsInVector(BzzVectorInt* dup);
	void DeleteNonZeroDuplicationsInVector(void);
	void DeleteDuplicationsInSortedVector(void);
	void DeleteNonZeroDuplicationsInSortedVector(void);
	friend void CopyDataFromVector(BzzMatrixInt* lval, BzzVectorInt& rval);
};

// Friend functions with default arguments
void ChangeDimensions(int dim, BzzVectorInt* result,
	char zero = 0);
char OperatorAndForSortedVectors(BzzVectorInt& v, BzzVectorInt& w,
	int f = BZZ_BIG_INT);
char OperatorAndForSortedVectors(BzzVectorInt& v, BzzVectorInt& w,
	BzzVectorInt* x, int f = BZZ_BIG_INT);


#endif // BZZ_VECTOR_INT_HPP