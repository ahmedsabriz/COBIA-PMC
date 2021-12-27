// TODO: ChangeDimensions
// BZZMATH: Release 7.0

//	====================< BzzMatrixSparseSymmetric.hpp >==================
//	* BzzMatrixSparseSymmetric: Class for sparse and symmetrical matrices*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixSparseSymmetric\MatrixSparseSymmetric.cpp		*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	11-1994	Date Written.
//	11-1994	Conversion to double precision done.
//	02-1994	Added the Bzz prefix to the names of the classes.
//	01-1998	Added FindBand.
//	03-1998	Added BzzPrintStructure.
//	03-1998	Added BzzPrintExistingElements;
//	03-1998	Added OrderingSymmetric.

//	============================================================================
//	******	BzzMatrixSparseSymmetric constructors:								*
//	* BzzMatrixSparseSymmetric A; // default										*
//	* BzzMatrixSparseSymmetric A = B; // copy-initializer						*
//	* BzzMatrixSparseSymmetric A(n,n); // sizes									*
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
//	* A.FindBand(&low);																			*
//	****************************************************************************
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixSparse																*
//	****************************************************************************
//	***** Implemented operations :															*
//	* C = A + B;	// C = A + B																*
//	* A += B;		// A = A + B																*
//	* C = A - B;	// C = A - B																*
//	* A -= B;		// A = A - B																*
//	* Product(A,B,&C);	// C = AB															*
//	* Product(A,x,&y);	// y = Ax															*
//	* y = A*x;;		// y = Ax																	*
//	* Product(3.,&A);	// y = Ax																*
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

#ifndef SPARSYMD_HPP
#define SPARSYMD_HPP

//	============================================================================
//	=====================< class BzzMatrixSparseSymmetric >===============
//	============================================================================
class BzzMatrixSparseSymmetric : public BzzBaseClass
{
	friend class BzzFactorizedSparseCholesky;
	friend class BzzMatrixBandSymmetric;
	friend class BzzFactorizedBandSymmetricPositive;
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparseSymmetric* lval, BzzFactorizedSparseCholesky* rval);

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int numRows, numColumns;
	int whoAmI;
	int lowerBand;

	ElementBzzMatrixSparse** elementRow;

	double& InsertElement
	(ElementBzzMatrixSparse* elem, int row,
		int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSparseSymmetric& rval);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixSparseSymmetric(void);

	// copy constructor
	BzzMatrixSparseSymmetric(const BzzMatrixSparseSymmetric& rval);

	// sizing constructor
	BzzMatrixSparseSymmetric(int rows, int columns);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMatrixSparseSymmetric(void);

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
	void FindBand(int* low);

	//	============================================================================
	//	****************************< assignment operators >************************
	//	============================================================================
	BzzMatrixSparseSymmetric& operator =
		(const BzzMatrixSparseSymmetric& rval);

	//	============================================================================
	//	================================< OPERATIONS >==============================
	//	============================================================================

	//	============================================================================
	//	***********************************< Sum >**********************************
	//	============================================================================

	friend BzzMatrixSparseSymmetric operator +
		(const BzzMatrixSparseSymmetric& lval,
			const BzzMatrixSparseSymmetric& rval);

	BzzMatrixSparseSymmetric& operator +=
		(const BzzMatrixSparseSymmetric& rval);

	//	============================================================================
	//	**********************************< Difference >****************************
	//	============================================================================

	friend BzzMatrixSparseSymmetric operator -
		(const BzzMatrixSparseSymmetric& lval,
			const BzzMatrixSparseSymmetric& rval);

	BzzMatrixSparseSymmetric& operator -=
		(const BzzMatrixSparseSymmetric& rval);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
	friend void Product
	(BzzMatrixSparseSymmetric& lval,
		BzzVector& rval, BzzVector* result);

	friend BzzVector operator *(BzzMatrixSparseSymmetric& lval,
		BzzVector& rval);

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================

	//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	============================< Modifying Functions >=========================
	//	============================================================================
	void RemoveElement(int row, int column);
	friend void Delete(BzzMatrixSparseSymmetric* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	void CleanMatrix(double eps); // eliminates those <= eps
	void BzzPrintStructure(void);
	void BzzPrintExistingElements(void);
	void OrderingSymmetric(BzzMatrixSparseSymmetric* A,
		BzzVectorInt* ordRows);
};

void Solve(BzzMatrixSparseSymmetric& A, BzzVector& b,
	BzzVector* x);
/*
//	==============< BzzMatrixSparseSymmetricLocked.hpp >========================
//	* BzzMatrixSparseSymmetricLocked: Class for sparse and symmetrical matrices*
// ============================================================================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2010	Date Written.

class BzzMatrixSparseSymmetricLocked : public BzzBaseClass
	{
friend class BzzVectorArray;
friend class BzzQuadraticProgramming;
friend class BzzMatrixSparseLocked;
private:
	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double ETA;
	int numRows,numColumns,numElements;
	int whoAmI;
	int maxIter;
	double eps;
	BzzVector d,l;
//	BzzVectorInt r,c,iv;
	BzzVectorInt r,c;
	BzzVectorInt low;
/////////////////////////////////////////////
//======< for gradient an Hessian >==========
	BzzVector x0,x,h,fp,fa,fm,aux,aug;
	double f0;
	double (*ptrFun)(BzzVector &x);
	BzzVectorInt	numVariablesInEachEquation;
	BzzVectorIntArray rI;
	BzzVectorArray rV;
	BzzVectorInt auxi; // serve per i Build
/////////////////////////////////////////////

public:
	BzzVector g;
//	============================================================================
//	*****************************< constructors >*******************************
//	============================================================================
	// default constructor
	BzzMatrixSparseSymmetricLocked(void);

	// copy constructor
	BzzMatrixSparseSymmetricLocked(const BzzMatrixSparseSymmetricLocked &rval);

	// constructor
	BzzMatrixSparseSymmetricLocked(BzzVector *dd,
		BzzVectorInt *rr,BzzVectorInt *cc,BzzVector *vv);
	void operator()(BzzVector *dd,
		BzzVectorInt *rr,BzzVectorInt *cc,BzzVector *vv);
//	============================================================================
//	********************************< Product >*********************************
//	============================================================================
	friend void Product
		(BzzMatrixSparseSymmetricLocked &lval,
			BzzVector &rval,BzzVector *result);
	friend BzzVector operator *(BzzMatrixSparseSymmetricLocked &lval,
			BzzVector &rval);

//	============================================================================
//	*****************************< Access functions >***************************
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
	void SetMaxIterations(int max)
		{maxIter = max;}
	void SetPrecision(double e)
		{eps = e;}

//	============================================================================
//	=============================< Modifying functions >========================
//	============================================================================
	BzzMatrixSparseSymmetricLocked &operator =
	(BzzMatrixSparseSymmetricLocked &rval);
	void BuildGradientAndHessian(BzzVector &x00,double f0,
		double (*F)(BzzVector &x));
	void BuildGradientAndDiagonalHessian(BzzVector &x00,double f0,
		double (*F)(BzzVector &x));
	void BuildLowerDiagonalHessian(void);
	void BuildForwardGradient(BzzVector &x00,double f0,
		double (*F)(BzzVector &x));
	void BuildForwardGradientAndDiagonalHessian(BzzVector &x00,double f0,
		double (*F)(BzzVector &x));
	void BuildGradientAndUpdateHessian(BzzVector &x00,double f0,
		double (*F)(BzzVector &x));

	void GetGradient(BzzVector *g);
	void MatrixReordering(BzzVectorInt &ordering);
	void SetLowerArray(void); // TODO
	friend void Swap(BzzMatrixSparseSymmetricLocked *A,BzzMatrixSparseSymmetricLocked *B);

//	============================================================================
//	=========================< Non-modifying functions >========================
//	============================================================================
	void GetLowerEnvelope(BzzVectorInt *lowerEnvelope);
	void GetNumberOfVariablesInEachEquation(BzzVectorInt *iiv);
	void BzzCuthillMcKleeOrdering(BzzVectorInt *ordering);
	void BzzReverseCuthillMcKleeOrdering(BzzVectorInt *ordering);
	int CountElementsInsideLowerBand(int lowerB);
	int CountElementsOutsideLowerBand(int lowerB);
	void GetDiagonal(int i,BzzVector *diag);
	int GetNumElements(void){return numElements;}

//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);
//	friend int Solve(BzzMatrixSparseSymmetricLocked &A,BzzVector &b,BzzVector *x,
//			 double eps = 1.e-15,int maxIter = 0,BzzVector *gg = 0);
	friend int Solve(BzzMatrixSparseSymmetricLocked &A,BzzVector &b,BzzVector *x,
		BzzVector *gg = 0);

//			 double eps = 1.e-15,int maxIter = 0,BzzVector *gg = 0);
	};

*/
#endif // SPARSYMD_HPP