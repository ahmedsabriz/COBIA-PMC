// BZZMATH: Release 7.0

// =========================< BzzMatrixLeft.hpp >========================
// * Class BzzMatrixLeft for Left (Lower) square matrices					*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 12)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 12)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixLeft\MatrixLeft.cpp									*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
// 04-1991	Date Written
// 11-1992	English version
// 01-1994	Modified Max, Min, MaxAbs, MinAbs.
// 03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
// 09-1994	Added shadow variable for returning object.
// 09-1994	Added functions ObjectCount, ObjectCountInScope.
// 10-1994	Minor changes throughout.
// 11-1994	Conversion to double precision done.
// 02-1995	Added the Bzz prefix to the names of the classes.
// 09-1995	Added TProduct, PtT,IPt,PtI, ITPt, PtIT, IPtT, TPtI

// Version 7.0
// 04-2013	Added BzzMatrixLeftSparse class
// 04-2013	Added Solve(lowerBound,L,b,&x) function

// ============================================================================
// ****** BzzMatrixLeft constructors:												*
// * BzzMatrixLeft L; // default														*
// * BzzMatrixLeft L = left; // copy-initializer								*
// * BzzMatrixLeft L(3,3); // sizes and places at 0							*
// * BzzMatrixLeft L(3,3,																*
// *	 1.,																							*
// *	 2.,3.,																						*
// *	 4.,5.,6.);// sizes and initialises													*
// * float x[]=																					*
// *		{																							*
// *		1.,																						*
// *		2.,3.,																					*
// *		4.,5.,6.																					*
// *		};																							*
// * BzzMatrixLeft L(3,3,x); // from array										*
// * BzzMatrixLeft L("LEFT.DAT"); // Formatted File							*
// * BzzMatrixLeft L('*',LEFT.BIN"); // Binary File							*
// ****************************************************************************
// ***** Access functions :																	*
// *	i = L.Rows(); // numRows																*
// *	i = L.Columns(); // numColumns														*
// *	xf = L.GetValue(i);																		*
// *	xf = L(i,j);																				*
// *	L(i,j) = xf;																				*
// *	xf = L[i][j];																				*
// *	L[i][j] = xf;																				*
// *	L.SetValue(i,j,xf);																		*
// *	int who = L.WhoAmI();																	*
// *	int count = BzzMatrixLeft::ObjectCount();									*
// *	int countInScope = BzzMatrixLeft::ObjectCountInScope();				*
// ****************************************************************************
// ***** Assignments:																			*
// *	L = left; // left BzzMatrixLeft												*
// ****************************************************************************
// *****	Operators for tests:																	*
// * if(L1 == L2)																					*
// * if(L1 != L2)																					*
// ****************************************************************************
// ***** BzzMessage and BzzPrint																*
// * L.BzzPrint("Comment"); // always														*
// * L.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
// ****************************************************************************
// ***** Norms																						*
// * xf = L.NormT();																				*
// * xf = L.NormR();																				*
// * xf = L.NormC();																				*
// * xf = L.NormF();																				*
// * xf = L.NormI();																				*
// * xf = L.Norm1();																				*
// * xf = L.Norm2();																				*
// ****************************************************************************
// ***** Max and Min																				*
// * xf = L.Max(&imax,&jmax);																	*
// * xf = L.MaxAbs(&imax,&jmax);																*
// * xf = L.Min(&imin,&jmin);																	*
// * xf = L.MinAbs(&imin,&jmin);																*
// ****************************************************************************
// ***** Save and Load																			*
// * L.Save("LEFT.DAT");	// formatted													*
// * L.Save('*',"LEFT.DAT");	// unformatted												*
// * Load(&L,"LEFT.DAT");																		*
// * Load(&L,'*',"LEFT.DAT");																	*
// ****************************************************************************
// ***** Delete, ChangeDimensions and Swap												*
// * Delete(&L);																					*
// * ChangeDimensions(newr,newc,&L);														*
// * Swap(&L1,&L2);																				*
// ****************************************************************************
// ***** Implemented operations :															*
//	* Sum(A,L,&B);				// B = A + L													*
//	* Sum(&A,L);				// A = A + L													*
//	* Sum(L,A,&B);				// B = L + A													*
//	* Sum(L,&A);				// A = L + A													*

//	* Difference(A,L,&B);	// B = A - L													*
//	* Difference(&A,L);		// A = A - L													*
//	* Difference(L,A,&B);	// B = L - A													*
//	* Difference(L,&A);		// A = L - A													*

//	* Product(A,L,&B);		// B = AL														*
//	* Product(&A,L);			// A = AL														*
//	* Product(L,A,&B);		// B = LA														*
//	* Product(L,&A);			// A = LA														*

//	* TProduct(A,L,&B);		// B = ATL														*
//	* TProduct(&A,L);			// A = ATL														*
//	* TProduct(L,A,&B);		// B = LTA														*
//	* TProduct(L,&A);			// A = LTA														*

//	* ProductT(A,L,&B);		// B = ALT														*
//	* ProductT(&A,L);			// A = ALT														*
//	* ProductT(L,A,&B);		// B = LAT														*
//	* ProductT(L,&A);			// A = LAT														*

//	* IProduct(A,L,&B);		// B = A-1L														*
//	* IProduct(&A,L);			// A = A-1L														*
//	* IProduct(L,A,&B);		// B = L-1A														*
//	* IProduct(L,&A);			// A = L-1A														*

//	* ProductI(A,L,&B);		// B = AL-1														*
//	* ProductI(&A,L);			// A = AL-1														*
//	* ProductI(L,A,&B);		// B = LA-1														*
//	* ProductI(L,&A);			// A = LA-1														*

//	* ITProduct(A,L,&B);		// B = A-TL														*
//	* ITProduct(&A,L);		// A = A-TL														*
//	* ITProduct(L,A,&B);		// B = L-TA														*
//	* ITProduct(L,&A);		// A = L-TA														*

//	* ProductIT(A,L,&B);		// B = AL-T														*
//	* ProductIT(&A,L);		// A = AL-T														*
//	* ProductIT(L,A,&B);		// B = LA-T														*
//	* ProductIT(L,&A);		// A = LA-T														*

//	* IProductT(A,L,&B);		// B = A-1LT													*
//	* IProductT(&A,L);		// A = A-1LT													*
//	* IProductT(L,A,&B);		// B = L-1AT													*
//	* IProductT(L,&A);		// A = L-1AT													*

//	* TProductI(A,L,&B);		// B = ATL-1													*
//	* TProductI(&A,L);		// A = ATL-1													*
//	* TProductI(L,A,&B);		// B = LTA-1													*
//	* TProductI(L,&A);		// A = LTA-1													*

// * Sum(L1,L2,&L3);			// L3 = L1 + L2;												*
// * L3 = L1 + L2;			 // L3 = L1 + L2;												*
// * L1 += L2;					// L1 = L1 + L2;												*
// * L1 = L1 + L2;			 // L1 = L1 + L2;												*
// * Sum(L1,L2,&L1);			// L1 = L1 + L2;												*
// * Sum(&L1,L2);				// L1 = L1 + L2;												*
// * Sum(L1,L2,&L2);			// L2 = L1 + L2;												*
// * Sum(L1,&L2);				// L2 = L1 + L2;												*
// * Sum(L1,L1,&L1);			// L1 = L1 + L1;												*
// * Sum(&L1);					// L1 = L1 + L1;												*
// * Difference(L1,L2,&L3); // L3 = L1 - L2;												*
// * L3 = L1 - L2;			 // L3 = L1 - L2;												*
// * Difference(L1,L2,&L1); // L1 = L1 - L2;												*
// * Difference(&L1,L2);	 // L1 = L1 - L2;												*
// * L1 -= L2;					// L1 = L1 - L2;												*
// * Difference(L1,L2,&L2); // L2 = L1 - L2;												*
// * Difference(L1,&L2);	 // L2 = L1 - L2;												*
// * Difference(L1,L1,&L1); // L1 = L1 - L1;												*
// * Difference(&L1);		 // L1 = L1 - L1;												*
// * Minus(L1,&L2);			// L2 = -L1;													*
// * L2 = -L1;					// L2 = -L1;													*
// * Minus(&L1);				// L1 = -L1;													*
// * L3 = L1*L2;																					*
// * y = L1*x;																						*
// * Product(L,x,&y);		// y = Lx;														*
// * Product(L,&x);			// x = Lx;														*
// * L2 = 3.*L1;																					*
// * L1 *= 3.;					// L1 = 3.*L1;													*
// * y = L%x;					// y = LTx;														*
// * TProduct(L,x,&y);		// y = LTx;														*
// * TProduct(L,&x);			// x = LTx;														*
// * L.SelfTProduct(&S);	// S = LTL														*
// * L.SelfProductT(&S);	// S = LLT														*
// * IProduct(L,x,&y);		// y = L-1x;													*
// * IProduct(L,&x);			// x = L-1x;													*
// * ITProduct(L,x,&y);		// y = L-1x;													*
// * ITProduct(L,&x);			// x = L-1x;												*
// * L2 = L1/3.;				// L2 = L1/3.													*
// * L1 /= 3.;					// L1 = L1/3.													*
// * Sum(L,R,&A);				// A = L + R;													*
// * Sum(R,L,&A);				// A = R + L;													*
// ****************************************************************************
// ***** Functions for linear systems:														*
// * Solve(L,b,&x);																				*
// * Solve(L,&bx);	else Solve(L,b,&b);													*
// * Solve(L,B,&X);																				*
// * Solve(l,&BX);	else Solve(L,B,&B);													*
// * Solve(lowerBound,L,b,&x);																*
// * TransposeSolve(L,b,&x);																	*
// * TransposeSolve(L,&bx); else TransposeSolve(L,b,&b);								*
// * TransposeSolve(L,B,&X);																	*
// * TransposeSolve(L,&BX); else TransposeSolve(L,B,&B);								*
// * L.Determinant();																			*
// * L.ConditionNumber();																		*
// * Inverse(&L);																					*
// ****************************************************************************

#ifndef BZZ_LEFT_DOUBLE_HPP
#define BZZ_LEFT_DOUBLE_HPP

// preventive declarations
class BzzVector;
class BzzMatrix;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;

// ============================================================================
// ======================< class BzzMatrixLeft >=========================
// ============================================================================
class BzzMatrixLeft : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrix;
	friend class BzzMatrixRight;
	friend class BzzMatrixSymmetric;
	friend class BzzSave;
	friend class BzzLoad;

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	double** matrix;
	int numRows, numColumns;
	int size;
	int whoAmI;
	char shadow;

	// initialisation constructors
	void Initialize(int rows, int columns);

	// deinitialisation
	void Deinitialize(void);

	// BzzMatrixLeft A('*',3,3);
	BzzMatrixLeft(char, int rows, int columns);

public:
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft L
	BzzMatrixLeft(void);

	// copy-initializer // BzzMatrixLeft L = left;
	BzzMatrixLeft(BzzMatrixLeft& rval);

	// sizes and initialises at 0
	// BzzMatrixLeft L(3,3);
	BzzMatrixLeft(int rows, int columns);

	// sizes and initialises
	// BzzMatrixLeft L(2,2,1.,2.,3.);
	BzzMatrixLeft(int rows, int columnsn, double a11, ...);

	// from array // BzzMatrixLeft L(3,3,w)
	BzzMatrixLeft(int rows, int columns, double* initvalues);

	// from formatted File // BzzMatrixLeft L("LEFT.DAT");
	BzzMatrixLeft(char* filematrix);

	// from binary File // BzzMatrixLeft L('*',"LEFT.BIN");
	BzzMatrixLeft(char, char* filematrix);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixLeft(void);

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

	// receives the values with control
	double GetValue(int row, int col) const;

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// assigns values with control
	void SetValue(int row, int col, double val);

	// assigns and receives vector values with control
	double& operator () (int row, int col);

	// assigns and receives vector values without control
	double* operator [](int i)
	{
		return matrix[i];
	}

	// ============================================================================
	// **************************< assignment operators >**************************
	// ============================================================================
	BzzMatrixLeft& operator =
		(const BzzMatrixLeft& rval);

	// ============================================================================
	// **************************< operators for tests >***************************
	// ============================================================================
	friend char operator ==
		(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval);

	friend char operator !=
		(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	//	============================================================================
	//	********************************< Sum >*********************************
	//	============================================================================
		// Sum(A,D,&B);
	friend void Sum
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// Sum(&A,D);
	friend void Sum
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// Sum(D,A,&B);
	friend void Sum
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// Sum(D,&A);
	friend void Sum
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< Difference >*********************************
	//	============================================================================
		// Difference(A,D,&B);
	friend void Difference
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// Difference(&A,D);
	friend void Difference
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// Difference(D,A,&B);
	friend void Difference
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// Difference(D,&A);
	friend void Difference
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
		// Product(A,D,&B);
	friend void Product
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// Product(&A,D);
	friend void Product
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// Product(D,A,&B);
	friend void Product
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// Product(D,&A);
	friend void Product
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< TProduct >*********************************
	//	============================================================================
		// TProduct(A,D,&B);
	friend void TProduct
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// TProduct(&A,D);
	friend void TProduct
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// TProduct(D,A,&B);
	friend void TProduct
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// TProduct(D,&A);
	friend void TProduct
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,D,&B);
	friend void ProductT
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// ProductT(&A,D);
	friend void ProductT
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// ProductT(D,A,&B);
	friend void ProductT
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductT(D,&A);
	friend void ProductT
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< IProduct >*********************************
	//	============================================================================
		// IProduct(A,D,&B);
	friend void IProduct
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// IProduct(&A,D);
	friend void IProduct
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// IProduct(D,A,&B);
	friend void IProduct
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// IProduct(D,&A);
	friend void IProduct
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< ProductI >*********************************
	//	============================================================================
		// ProductI(A,D,&B);
	friend void ProductI
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// ProductI(&A,D);
	friend void ProductI
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// ProductI(D,A,&B);
	friend void ProductI
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductI(D,&A);
	friend void ProductI
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< ITProduct >********************************
	//	============================================================================
		// ITProduct(A,D,&B);
	friend void ITProduct
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// ITProduct(&A,D);
	friend void ITProduct
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// ITProduct(D,A,&B);
	friend void ITProduct
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// ITProduct(D,&A);
	friend void ITProduct
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< ProductIT >*******************************
	//	============================================================================
		// ProductIT(A,D,&B);
	friend void ProductIT
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// ProductIT(&A,D);
	friend void ProductIT
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// ProductIT(D,A,&B);
	friend void ProductIT
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductIT(D,&A);
	friend void ProductIT
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< IProductT >*******************************
	//	============================================================================
		// IProductT(A,D,&B);
	friend void IProductT
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// IProductT(&A,D);
	friend void IProductT
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// IProductT(D,A,&B);
	friend void IProductT
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// IProductT(D,&A);
	friend void IProductT
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< TProductI >*******************************
	//	============================================================================
		// TProductI(A,D,&B);
	friend void TProductI
	(BzzMatrix& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// TProductI(&A,D);
	friend void TProductI
	(BzzMatrix* lvalAndResult, BzzMatrixLeft& rval);

	// TProductI(D,A,&B);
	friend void TProductI
	(BzzMatrixLeft& lval, BzzMatrix& rval, BzzMatrix* result);

	// TProductI(D,&A);
	friend void TProductI
	(BzzMatrixLeft& rval, BzzMatrix* rvalAndResult);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================
		 // Sum(L1,L2,&L3);
	friend void Sum
	(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval,
		BzzMatrixLeft* result);

	// L3 = L1 + L2;
	friend BzzMatrixLeft operator +
		(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval);

	// Sum(L,R,&A);
	friend void Sum
	(BzzMatrixLeft& lval, BzzMatrixRight& rval,
		BzzMatrix* result);

	// Sum(R,L,&A);
	friend void Sum
	(BzzMatrixRight& lval, BzzMatrixLeft& rval,
		BzzMatrix* result);

	// Sum(&L1,L2);
	friend void Sum
	(BzzMatrixLeft* lvalAndResult, const BzzMatrixLeft& rval);

	// L1 += L2;
	BzzMatrixLeft& operator +=
		(const BzzMatrixLeft& rval);

	// Sum(L1,&L2);
	friend void Sum
	(const BzzMatrixLeft& lval, BzzMatrixLeft* rvalAndResult);

	// Sum(&L);
	friend void Sum
	(BzzMatrixLeft* lvalRvalAndResult);

	// ============================================================================
	// *****************************< Difference >*********************************
	// ============================================================================
		// Difference(L1,L2,&L3);
	friend void Difference
	(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval,
		BzzMatrixLeft* result);

	// L3 = L1 - L2;
	friend BzzMatrixLeft operator -
		(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval);

	// Difference(&L1,L2);
	friend void Difference
	(BzzMatrixLeft* lvalAndResult, const BzzMatrixLeft& rval);

	// L1 -= L2;
	BzzMatrixLeft& operator -=
		(const BzzMatrixLeft& rval);

	// Difference(L1,&L2);
	friend void Difference
	(const BzzMatrixLeft& lval, BzzMatrixLeft* rvalAndResult);

	// ============================================================================
	// ********************************< Minus >***********************************
	// ============================================================================
		// Minus(L1,&L2);
	friend void Minus
	(const BzzMatrixLeft& rval, BzzMatrixLeft* result);

	// L2 = -L1;
	friend BzzMatrixLeft operator - // unary
		(const BzzMatrixLeft& rval);

	// Minus(&L1);
	friend void Minus
	(BzzMatrixLeft* rvalAndResult);

	// ============================================================================
	// ********************************< Product >*********************************
	// ============================================================================
		// L3 = L1*L2
	friend BzzMatrixLeft operator *
		(const BzzMatrixLeft& lval, const BzzMatrixLeft& rval);

	// Product(L,x,&y); y = L*x;
	friend void Product
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(L,&x); x = L*x;
	friend void Product(const BzzMatrixLeft& lval,
		BzzVector* rvalAndResult);

	// y = L*x;
	friend BzzVector operator *
		(const BzzMatrixLeft& lval, const BzzVector& rval);

	// L2 = 3.*L1;
	friend BzzMatrixLeft operator *
		(double lval, const BzzMatrixLeft& rval);

	// L *= 3.;
	BzzMatrixLeft& operator *=
		(double rval);

	// ============================================================================
	// ********************************< TProduct >********************************
	// ============================================================================
		// L.SelfTProduct(L,&S);
	void SelfTProduct(BzzMatrixSymmetric* S);

	// y = L%x;
	friend BzzVector operator %
		(BzzMatrixLeft& lval, BzzVector& rval);

	// TProduct(L,x,&y);
	friend void TProduct
	(BzzMatrixLeft& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(L,&x);
	friend void TProduct
	(BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ******************************< SelfProductT >******************************
	// ============================================================================
		// L.SelfProductT(&S);
	void SelfProductT(BzzMatrixSymmetric* S);

	// ============================================================================
	// ********************************< IProduct >********************************
	// ============================================================================

		// IProduct(L,x,&y);
	friend void IProduct
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(L,&x);
	friend void IProduct
	(const BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// *******************************< ITProduct >********************************
	// ============================================================================
		// ITProduct(L,x,&y);
	friend void ITProduct
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// ITProduct(L,&x);
	friend void ITProduct
	(const BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ********************************< Division >********************************
	// ============================================================================
		// A = B/3.;
	friend BzzMatrixLeft operator /
		(const BzzMatrixLeft& lval, double rval);

	// A /= 3.;
	BzzMatrixLeft& operator /=
		(double rval);

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// *********************************< Save >***********************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

// ******************************< Max and Min >*******************************
	double Max(int* imax = 0, int* jmax = 0);
	double MaxAbs(int* imax = 0, int* jmax = 0);

	double Min(int* imin = 0, int* jmin = 0);
	double MinAbs(int* imin = 0, int* jmin = 0);

	// ********************************< Norms >***********************************
	double NormT(void);
	double NormR(void);
	double NormC(void);
	double NormF(void);
	double NormI(void) { return NormR(); }
	double Norm1(void) { return NormC(); }
	double Norm2(void);

	// ============================================================================
	// ========================< Modifying Functions >=============================
	// ============================================================================
	friend void Delete(BzzMatrixLeft* L); // eliminates BzzMatrixLeft
	friend void ChangeDimensions(int rows, int columns,
		BzzMatrixLeft* result, char);
	// from formatted file
	friend void Load(BzzMatrixLeft* A, char* filematrix);
	// from binary file
	friend void Load(BzzMatrixLeft* A, char, char* filematrix);
	friend char Inverse(BzzMatrixLeft* L);
	friend void Swap(BzzMatrixLeft* lval, BzzMatrixLeft* rval);

	// ============================================================================
	// =======================< Functions for linear algebra >=====================
	// ============================================================================
	friend void Solve(const BzzMatrixLeft& L, BzzVector* bx);
	friend void Solve(const BzzMatrixLeft& L,
		const BzzVector& b, BzzVector* x);
	friend void Solve(const BzzMatrixLeft& L, BzzMatrix* BX);
	friend void Solve(const BzzMatrixLeft& L,
		const BzzMatrix& B, BzzMatrix* X);
	friend void Solve(BzzVectorInt& lowerBound, BzzMatrixLeft& L,
		BzzVector* bx);
	friend void Solve(BzzVectorInt& lowerBound, BzzMatrixLeft& L,
		BzzVector& b, BzzVector* x);
	friend void TransposeSolve(const BzzMatrixLeft& L,
		BzzVector* bx);
	friend void TransposeSolve(const BzzMatrixLeft& L,
		const BzzVector& b, BzzVector* x);
	friend void TransposeSolve(const BzzMatrixLeft& L,
		BzzMatrix* BX);
	friend void TransposeSolve(const BzzMatrixLeft& L,
		const BzzMatrix& B, BzzMatrix* X);

	double Determinant(void);
	double ConditionNumber();
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzMatrixLeft* result, char zero = 0);


#endif // BZZ_LEFT_DOUBLE_HPP