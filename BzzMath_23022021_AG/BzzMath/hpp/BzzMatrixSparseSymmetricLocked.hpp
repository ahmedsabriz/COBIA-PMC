// BZZMATH: Release 7.0

//===============================================================================
//======================================================================================
//	==============< BzzMatrixSparseSymmetricLocked.hpp >========================
//	* BzzMatrixSparseSymmetricLocked: Class for sparse and symmetrical matrices*
// ============================================================================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2010	Date Written.

#ifndef BZZ_MATRIX_SPARSE_SYMMETRICLLOCKED
#define BZZ_MATRIX_SPARSE_SYMMETRICLLOCKED

//	============================================================================
//	=================< class BzzBuildGradientAndHessian >=======================
//	============================================================================

class BzzBuildGradientAndHessian
{
private:
	//	double zerDer,eta2,hj,xh,xdh,hInv;
	double hi, fpi;
	BzzVector xaaa, xbbb;
public:
	double (*ptrFunction)(BzzVector& xx);
	//	void GetBuildGradientAndHessian(int jStart,int jEnd,BzzVector &x00,double f0,
	//		double (*F)(BzzVector &x));
	//	void GetBuildGradientAndPositiveHessian(int jStart,int jEnd,BzzVector &x00,double f0,
	//		double (*F)(BzzVector &x));
	//
	void GetBuildGradientAndDiagonalHessian(int iStart, int iEnd, BzzVector& x,
		double f0, BzzVector& h, BzzVector& fp, BzzVector& fa,
		BzzVector& gg, BzzVector& dd);
	void GetBuildLowerDiagonalHessian(int iStart, int iEnd, BzzVector& x,
		double f0, BzzVector& h, BzzVector& fp, BzzVector& fa,
		BzzVector& fm, BzzVectorInt& r, BzzVectorInt& c, BzzVector& l);
	void GetBuildForwardGradient(int iStart, int iEnd, BzzVector& x,
		double f0, BzzVector& gg);
	void GetBuildForwardGradientAndDiagonalHessian(int iStart,
		int iEnd, BzzVector& x, double f0, BzzVector& h,
		BzzVector& fp, BzzVector& fa, BzzVector& gg, BzzVector& dd);
	//	void GetBuildGradientAndUpdateHessian(int jStart,int jEnd,BzzVector &x00,double f0,
	//		double (*F)(BzzVector &x));
	//	void GetBuildGradientAndUpdatePositiveHessian(int jStart,int jEnd,BzzVector &x00,double f0,
	//		double (*F)(BzzVector &x));
	//
	//	void GetJacobian(int jStart,int jEnd,BzzVector &x,BzzVector &xDimensions,
	//		BzzVector &fi,BzzMatrix &Ji);
};

class BzzMatrixSparseSymmetricLocked : public BzzBaseClass
{
	friend class BzzVectorArray;
	friend class BzzQuadraticProgramming;
	friend class BzzMatrixSparseLocked;
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double ETA;
	int numRows, numColumns, numThread;
	int whoAmI;
	int maxIter;
	double eps;
	//	BzzVectorInt r,c,iv;
	BzzVectorInt low;
	/////////////////////////////////////////////
	//======< for gradient an Hessian >==========
	BzzBuildGradientAndHessian* getGradAndHess;
	BzzVector x0, x, h, fp, fa, fm, aux, aug;
	double f0;
	double (*ptrFun)(BzzVector& x);
	BzzVectorInt	numVariablesInEachEquation;
	BzzVectorIntArray rI;
	BzzVectorArray rV;
	BzzVectorInt auxi; // serve per i Build
/////////////////////////////////////////////

public:
	int numElements;
	BzzVector d, l;
	BzzVectorInt r, c;
	BzzVector g;

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixSparseSymmetricLocked(void);

	// copy constructor
	BzzMatrixSparseSymmetricLocked(const BzzMatrixSparseSymmetricLocked& rval);

	// constructor
	BzzMatrixSparseSymmetricLocked(BzzVector* dd,
		BzzVectorInt* rr, BzzVectorInt* cc, BzzVector* vv);
	void operator()(BzzVector* dd,
		BzzVectorInt* rr, BzzVectorInt* cc, BzzVector* vv);

	//	============================================================================
	//	*********************************< destructor >*****************************
	//	============================================================================
	~BzzMatrixSparseSymmetricLocked(void)
	{
		if (numThread != 0)
#if BZZ_COMPILER != 0
			delete[] getGradAndHess;
#else
			//delete[numThread] getGradAndHess;
			delete[] getGradAndHess;
#endif
	}

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
	friend void Product
	(BzzMatrixSparseSymmetricLocked& lval,
		BzzVector& rval, BzzVector* result);
	friend BzzVector operator *(BzzMatrixSparseSymmetricLocked& lval,
		BzzVector& rval);

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
	void SetMaxIterations(int max)
	{
		maxIter = max;
	}
	void SetPrecision(double e)
	{
		eps = e;
	}

	//	============================================================================
	//	=============================< Modifying functions >========================
	//	============================================================================
	BzzMatrixSparseSymmetricLocked& operator =
		(BzzMatrixSparseSymmetricLocked& rval);

	void BuildGradientAndHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));
	void BuildGradientAndPositiveHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));

	void BuildGradientAndDiagonalHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));
	void BuildLowerDiagonalHessian(void);
	void BuildForwardGradient(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));
	void BuildForwardGradientAndDiagonalHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));
	void BuildGradientAndUpdateHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));
	void BuildGradientAndUpdatePositiveHessian(BzzVector& x00, double f0,
		double (*F)(BzzVector& x));

	void GetGradient(BzzVector* g);
	void MatrixReordering(BzzVectorInt& ordering);
	void SetLowerArray(void); // TODO
	friend void Swap(BzzMatrixSparseSymmetricLocked* A, BzzMatrixSparseSymmetricLocked* B);

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================
	void GetLowerEnvelope(BzzVectorInt* lowerEnvelope);
	void GetNumberOfVariablesInEachEquation(BzzVectorInt* iiv);
	void BzzCuthillMcKleeOrdering(BzzVectorInt* ordering);
	void BzzReverseCuthillMcKleeOrdering(BzzVectorInt* ordering);
	int CountElementsInsideLowerBand(int lowerB);
	int CountElementsOutsideLowerBand(int lowerB);
	void GetDiagonal(int i, BzzVector* diag);
	int GetNumElements(void) { return numElements; }

	//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);
	//	friend int Solve(BzzMatrixSparseSymmetricLocked &A,BzzVector &b,BzzVector *x,
	//			 double eps = 1.e-15,int maxIter = 0,BzzVector *gg = 0);
	friend int Solve(BzzMatrixSparseSymmetricLocked& A, BzzVector& b, BzzVector* x,
		BzzVector*);
	friend int Solve(BzzMatrixSparseSymmetricLocked& A, int dim, BzzVector& b, BzzVector* x,
		BzzVector*);

	//			 double eps = 1.e-15,int maxIter = 0,BzzVector *gg = 0);
};

// Friend functions with default arguments
int Solve(BzzMatrixSparseSymmetricLocked& A, BzzVector& b, BzzVector* x,
	BzzVector* gg = 0);
int Solve(BzzMatrixSparseSymmetricLocked& A, int dim, BzzVector& b, BzzVector* x,
	BzzVector* gg = 0);


#endif // BZZ_MATRIX_SPARSE_SYMMETRICLLOCKED