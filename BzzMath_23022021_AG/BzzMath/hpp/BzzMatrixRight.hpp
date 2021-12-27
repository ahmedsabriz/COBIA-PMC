// BZZMATH: Release 7.0

// ========================< BzzMatrixRight.hpp >========================
// * Class BzzMatrixRight for Right (Upper) square matrices					*
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
// *				MatrixRight\MatrixRight.cpp								*
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

// ============================================================================
// ****** BzzMatrixRight constructors:												*
// * BzzMatrixRight R; // default													*
// * BzzMatrixRight R = right; // copy-initializer								*
// * BzzMatrixRight L(3,3); // sizes and places at 0							*
// * BzzMatrixRight R(3,3,																*
// *		 1.,2.,3.,																				*
// *			 4.,5.,																				*
// *				 6.);// sizes and initialises												*
// * float x[]=																					*
// *		{																							*
// *		1.,2.,3.,																				*
// *			4.,5.,																				*
// *				6.																					*
// *		};																							*
// * BzzMatrixRight R(3,3,x); // from array										*
// * BzzMatrixRight R("MAT.DAT"); // Formatted File							*
// * BzzMatrixRight R('*',MAT.BIN"); // Binary File							*
// ****************************************************************************
// ***** Access functions :																	*
// *	i = R.Rows(); // numRows																*
// *	i = R.Columns(); // numColumns														*
// *	xf = R.GetValue(i);																		*
// *	xf = R(i,j);																				*
// *	R(i,j) = xf;																				*
// *	xf = R[i][j];																				*
// *	R[i][j] = xf;																				*
// *	R.SetValue(i,j,xf);																		*
// *	int who = R.WhoAmI();																	*
// *	int count = BzzMatrixRight::ObjectCount();								*
// *	int countInScope = BzzMatrixRight::ObjectCountInScope();				*
// ****************************************************************************
// ***** Assignments:																			*
// *	R = right; // right BzzMatrixRight											*
// ****************************************************************************
// *****	Operators for tests:										 							*
// * if(R1 == R2)																					*
// * if(R1 != R2)																					*
// ****************************************************************************
// ***** BzzMessage and BzzPrint																*
// * R.BzzPrint("Comment"); // always														*
// * R.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
// ****************************************************************************
// ***** Norms																						*
// * xf = R.NormT();																				*
// * xf = R.NormR();																				*
// * xf = R.NormC();																				*
// * xf = R.NormF();																				*
// * xf = R.NormI();																				*
// * xf = R.Norm1();																				*
// * xf = R.Norm2();																				*
// ****************************************************************************
// ***** Max and Min																				*
// * xf = R.Max(&imax,&jmax);																	*
// * xf = R.MaxAbs(&imax,&jmax);																*
// * xf = R.Min(&imin,&jmin);																	*
// * xf = R.MinAbs(&imin,&jmin);																*
// ****************************************************************************
// ***** Save and Load																			*
// * R.Save("RIGHT.DAT");	// formatted													*
// * R.Save('*',"RIGHT.DAT");	// unformatted												*
// * Load(&R,"RIGHT.DAT");																		*
// * Load(&R,'*',"RIGHT.DAT");																*
// ****************************************************************************
// ***** Delete, ChangeDimensions and Swap												*
// * Delete(&R);																					*
// * ChangeDimensions(newr,newc,&R);														*
// * Swap(&R1,&R2);																				*
// ****************************************************************************
// ***** Implemented operations :															*
//	* Sum(A,R,&B);				// B = A + R													*
//	* Sum(&A,R);				// A = A + R													*
//	* Sum(R,A,&B);				// B = R + A													*
//	* Sum(R,&A);				// A = R + A													*

//	* Difference(A,R,&B);	// B = A - R													*
//	* Difference(&A,R);		// A = A - R													*
//	* Difference(R,A,&B);	// B = R - A													*
//	* Difference(R,&A);		// A = R - A													*

//	* Product(A,R,&B);		// B = AR														*
//	* Product(&A,R);			// A = AR														*
//	* Product(R,A,&B);		// B = RA														*
//	* Product(R,&A);			// A = RA														*

//	* TProduct(A,R,&B);		// B = ATR														*
//	* TProduct(&A,R);			// A = ATR														*
//	* TProduct(R,A,&B);		// B = RTA														*
//	* TProduct(R,&A);			// A = RTA														*

//	* ProductT(A,R,&B);		// B = ART														*
//	* ProductT(&A,R);			// A = ART														*
//	* ProductT(R,A,&B);		// B = RAT														*
//	* ProductT(R,&A);			// A = RAT														*

//	* IProduct(A,R,&B);		// B = A-1R														*
//	* IProduct(&A,R);			// A = A-1R														*
//	* IProduct(R,A,&B);		// B = R-1A														*
//	* IProduct(R,&A);			// A = R-1A														*

//	* ProductI(A,R,&B);		// B = AR-1														*
//	* ProductI(&A,R);			// A = AR-1														*
//	* ProductI(R,A,&B);		// B = RA-1														*
//	* ProductI(R,&A);			// A = RA-1														*

//	* ITProduct(A,R,&B);		// B = A-TR														*
//	* ITProduct(&A,R);		// A = A-TR														*
//	* ITProduct(R,A,&B);		// B = R-TA														*
//	* ITProduct(R,&A);		// A = R-TA														*

//	* ProductIT(A,R,&B);		// B = AR-T														*
//	* ProductIT(&A,R);		// A = AR-T														*
//	* ProductIT(R,A,&B);		// B = RA-T														*
//	* ProductIT(R,&A);		// A = RA-T														*

//	* IProductT(A,R,&B);		// B = A-1RT													*
//	* IProductT(&A,R);		// A = A-1RT													*
//	* IProductT(R,A,&B);		// B = R-1AT													*
//	* IProductT(R,&A);		// A = R-1AT													*

//	* TProductI(A,R,&B);		// B = ATR-1													*
//	* TProductI(&A,R);		// A = ATR-1													*
//	* TProductI(R,A,&B);		// B = RTA-1													*
//	* TProductI(R,&A);		// A = RTA-1													*

// * Sum(R1,R2,&R3);			// R3 = R1 + R2;												*
// * R3 = R1 + R2;			 // R3 = R1 + R2;												*
// * R1 += R2;					// R1 = R1 + R2;												*
// * R1 = R1 + R2;			 // R1 = R1 + R2;												*
// * Sum(R1,R2,&R1);			// R1 = R1 + R2;												*
// * Sum(&R1,R2);				// R1 = R1 + R2;												*
// * Sum(R1,R2,&R2);			// R2 = R1 + R2;												*
// * Sum(R1,&R2);				// R2 = R1 + R2;												*
// * Sum(R1,R1,&R1);			// R1 = R1 + R1;												*
// * Sum(&R1);					// R1 = R1 + R1;												*
// * Difference(R1,R2,&R3); // R3 = R1 - R2;												*
// * R3 = R1 - R2;			 // R3 = R1 - R2;												*
// * Difference(R1,R2,&R1); // R1 = R1 - R2;												*
// * Difference(&R1,R2);	 // R1 = R1 - R2;												*
// * R1 -= R2;					// R1 = R1 - R2;												*
// * Difference(R1,R2,&R2); // R2 = R1 - R2;												*
// * Difference(R1,&R2);	 // R2 = R1 - R2;												*
// * Difference(R1,R1,&R1); // R1 = R1 - R1;												*
// * Difference(&R1);		 // R1 = R1 - R1;												*
// * Minus(R1,&R2);			// R2 = -R1;													*
// * R2 = -R1;					// R2 = -R1;													*
// * Minus(&R1);				// R1 = -R1;													*
// * R3 = R1*R2;																					*
// * y = R1*x;																						*
// * Product(R,x,&y);		// y = Rx;														*
// * Product(R,&x);			// x = Rx;														*
// * R2 = 3.*R1;																					*
// * R1 *= 3.;					// R1 = 3.*R1;													*
// * y = R%x;					// y = RTx;														*
// * TProduct(R,x,&y);		// y = RTx;														*
// * TProduct(R,&x);			// x = RTx;														*
// * R.SelfTProduct(&S);	// S = RTR														*
// * R.SelfProductT(&S);	// S = RRT														*
// * IProduct(R,x,&y);		// y = R-1x;													*
// * IProduct(R,&x);			// x = R-1x;													*
// * ITProduct(R,x,&y);		// y = R-Tx;													*
// * ITProduct(R,&x);			// x = R-Tx;												*
// * R2 = R1/3.;				// R2 = R1/3.													*
// * R1 /= 3.;					// R1 = R1/3.													*
// * Sum(L,R,&A);				// A = L + R;													*
// * Sum(R,L,&A);				// A = R + L;													*
// ****************************************************************************
// ***** Functions for linear algebra:														*
// * Solve(R,b,&x);																				*
// * Solve(R,&bx);	else Solve(R,b,&b);													*
// * Solve(R,B,&X);																				*
// * Solve(R,&BX);	else Solve(R,B,&B);													*
// * TranposeSolve(R,b,&x);																	*
// * TransposeSolve(R,&bx); else TransposeSolve(R,b,&b);								*
// * TranposeSolve(R,B,&X);																	*
// * TransposeSolve(R,&BX); else TransposeSolve(R,B,&B);								*
// * R.Determinant();																			*
// * R.ConditionNumber();																		*
// * Inverse(&R);																					*
// ****************************************************************************

#ifndef BZZ_RIGHT_DOUBLE_HPP
#define BZZ_RIGHT_DOUBLE_HPP

// preventive declarations
class BzzVector;
class BzzMatrix;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;

// ============================================================================
// ========================< class BzzMatrixRight >======================
// ============================================================================
class BzzMatrixRight : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrix;
	friend class BzzMatrixLeft;
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

	// BzzMatrixRight A('*',3,3);
	BzzMatrixRight(char, int rows, int columns);

public:
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixRight R
	BzzMatrixRight(void);

	// copy-initializer // BzzMatrixRight R = right;
	BzzMatrixRight(BzzMatrixRight& rval);

	// sizes and initialises to 0 // BzzMatrixRight R(3,3);
	BzzMatrixRight(int rows, int columns);

	// sizes and initialises
	// BzzMatrixRight R(2,2,1.,2.,3.);
	BzzMatrixRight(int rows, int columns, double a11, ...);

	// from array // BzzMatrixRight R(3,3,w)
	BzzMatrixRight(int rows, int columns, double* initvalues);

	// from formatted File // BzzMatrixRight R("RIGHT.DAT");
	BzzMatrixRight(char* filematrix);

	// from binary File // BzzMatrixRight R('*',"RIGHT.BIN");
	BzzMatrixRight(char, char* filematrix);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixRight(void);

	// ============================================================================
	// ********************< Non-modifying access functions >**********************
	// ============================================================================
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

	// receive the values with control
	double GetValue(int row, int col) const;

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// assigns values with control
	void SetValue(int row, int col, double val);

	// assigns and receives vector values with control
	double& operator () (int row, int col);

	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// ============================================================================
	// **************************< assignment operators >**************************
	// ============================================================================
	BzzMatrixRight& operator =
		(const BzzMatrixRight& rval);

	// ============================================================================
	// **************************< operators for tests >***************************
	// ============================================================================
	friend char operator ==
		(const BzzMatrixRight& lval, const BzzMatrixRight& rval);

	friend char operator !=
		(const BzzMatrixRight& lval, const BzzMatrixRight& rval);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	//	============================================================================
	//	********************************< Sum >*********************************
	//	============================================================================
		// Sum(A,D,&B);
	friend void Sum
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// Sum(&A,D);
	friend void Sum
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// Sum(D,A,&B);
	friend void Sum
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// Sum(D,&A);
	friend void Sum
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< Difference >*********************************
	//	============================================================================
		// Difference(A,D,&B);
	friend void Difference
	(BzzMatrix& lval, BzzMatrixRight& rval,
		BzzMatrix* result);

	// Difference(&A,D);
	friend void Difference
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// Difference(D,A,&B);
	friend void Difference
	(BzzMatrixRight& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Difference(D,&A);
	friend void Difference
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);
	// Difference(L,R,&A);

	friend void Difference
	(BzzMatrixLeft& lval, BzzMatrixRight& rval,
		BzzMatrix* result);

	// Difference(R,L,&A);
	friend void Difference
	(BzzMatrixRight& lval, BzzMatrixLeft& rval,
		BzzMatrix* result);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
		// Product(A,D,&B);
	friend void Product
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// Product(&A,D);
	friend void Product
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// Product(D,A,&B);
	friend void Product
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// Product(D,&A);
	friend void Product
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< TProduct >*********************************
	//	============================================================================
		// TProduct(A,D,&B);
	friend void TProduct
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// TProduct(&A,D);
	friend void TProduct
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// TProduct(D,A,&B);
	friend void TProduct
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// TProduct(D,&A);
	friend void TProduct
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,D,&B);
	friend void ProductT
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// ProductT(&A,D);
	friend void ProductT
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// ProductT(D,A,&B);
	friend void ProductT
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductT(D,&A);
	friend void ProductT
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< IProduct >*********************************
	//	============================================================================
		// IProduct(A,D,&B);
	friend void IProduct
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// IProduct(&A,D);
	friend void IProduct
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// IProduct(D,A,&B);
	friend void IProduct
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// IProduct(D,&A);
	friend void IProduct
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductI >*********************************
	//	============================================================================
		// ProductI(A,D,&B);
	friend void ProductI
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// ProductI(&A,D);
	friend void ProductI
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// ProductI(D,A,&B);
	friend void ProductI
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductI(D,&A);
	friend void ProductI
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< ITProduct >*******************************
	//	============================================================================
		// ITProduct(A,D,&B);
	friend void ITProduct
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// ITProduct(&A,D);
	friend void ITProduct
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// ITProduct(D,A,&B);
	friend void ITProduct
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// ITProduct(D,&A);
	friend void ITProduct
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< ProductIT >*******************************
	//	============================================================================
		// ProductIT(A,D,&B);
	friend void ProductIT
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// ProductIT(&A,D);
	friend void ProductIT
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// ProductIT(D,A,&B);
	friend void ProductIT
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// ProductIT(D,&A);
	friend void ProductIT
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< IProductT >*******************************
	//	============================================================================
		// IProductT(A,D,&B);
	friend void IProductT
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// IProductT(&A,D);
	friend void IProductT
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// IProductT(D,A,&B);
	friend void IProductT
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// IProductT(D,&A);
	friend void IProductT
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< TProductI >*********************************
	//	============================================================================
		// TProductI(A,D,&B);
	friend void TProductI
	(BzzMatrix& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// TProductI(&A,D);
	friend void TProductI
	(BzzMatrix* lvalAndResult, BzzMatrixRight& rval);

	// TProductI(D,A,&B);
	friend void TProductI
	(BzzMatrixRight& lval, BzzMatrix& rval, BzzMatrix* result);

	// TProductI(D,&A);
	friend void TProductI
	(BzzMatrixRight& rval, BzzMatrix* rvalAndResult);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================
		 // Sum(R1,R2,&R3);
	friend void Sum
	(const BzzMatrixRight& lval, const BzzMatrixRight& rval,
		BzzMatrixRight* result);

	// R3 = R1 + R2;
	friend BzzMatrixRight operator +
		(const BzzMatrixRight& lval, const BzzMatrixRight& rval);

	// Sum(L,R,&A);
	friend void Sum
	(BzzMatrixLeft& lval, BzzMatrixRight& rval, BzzMatrix* result);

	// Sum(R,L,&A);
	friend void Sum
	(BzzMatrixRight& lval, BzzMatrixLeft& rval, BzzMatrix* result);

	// Sum(&R1,R2);
	friend void Sum
	(BzzMatrixRight* lvalAndResult, const BzzMatrixRight& rval);

	// R1 += R2;
	BzzMatrixRight& operator +=
		(const BzzMatrixRight& rval);

	// Sum(R1,&R2);
	friend void Sum
	(const BzzMatrixRight& lval, BzzMatrixRight* rvalAndResult);

	// Sum(&R);
	friend void Sum
	(BzzMatrixRight* lvalRvalAndResult);

	// ============================================================================
	// *****************************< Difference >*********************************
	// ============================================================================
		// Difference(R1,R2,&R3);
	friend void Difference
	(const BzzMatrixRight& lval, const BzzMatrixRight& rval,
		BzzMatrixRight* result);

	// R3 = R1 - R2;
	friend BzzMatrixRight operator -
		(const BzzMatrixRight& lval, const BzzMatrixRight& rval);

	// Difference(&R1,R2);
	friend void Difference
	(BzzMatrixRight* lvalAndResult, const BzzMatrixRight& rval);

	// R1 -= R2;
	BzzMatrixRight& operator -=
		(const BzzMatrixRight& rval);

	// Difference(R1,&R2);
	friend void Difference
	(const BzzMatrixRight& lval, BzzMatrixRight* rvalAndResult);

	// ============================================================================
	// ********************************< Minus >***********************************
	// ============================================================================
		// Minus(R1,&R2);
	friend void Minus
	(const BzzMatrixRight& rval, BzzMatrixRight* result);

	// R2 = -R1;
	friend BzzMatrixRight operator -
		(const BzzMatrixRight& rval);

	// Minus(&R);;
	friend void Minus
	(BzzMatrixRight* rvalAndResult);

	// ============================================================================
	// ********************************< Product >*********************************
	// ============================================================================
		// R3 = R1*R2;
	friend BzzMatrixRight operator *
		(const BzzMatrixRight& lval, const BzzMatrixRight& rval);

	// y = R*x;
	friend BzzVector operator *
		(const BzzMatrixRight& lval, const BzzVector& rval);

	// Product(R,x,&y); y = R*x;
	friend void Product
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(R,&x); x = R*x;
	friend void Product
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// R2 = 3.*R1;
	friend BzzMatrixRight operator *
		(double lval, const BzzMatrixRight& rval);

	// R *= 3.;
	BzzMatrixRight& operator *=
		(double rval);

	// ============================================================================
	// ********************************< TProduct >********************************
	// ============================================================================
		// R.SelfTProduct(&S);
	void SelfTProduct(BzzMatrixSymmetric* S);

	// TProduct(R,x,&y);
	friend void TProduct
	(BzzMatrixRight& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(R,&x);
	friend void TProduct
	(BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// y = R%x;
	friend BzzVector operator %
		(BzzMatrixRight& lval, BzzVector& rval);

	// ============================================================================
	// *****************************< SelfProductT >*******************************
	// ============================================================================
		// R.SelfProductT(&S);
	void SelfProductT(BzzMatrixSymmetric* S);

	// ============================================================================
	// ******************************< IProduct >**********************************
	// ============================================================================
		// IProduct(R,x,&y);
	friend void IProduct
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(R,&x);
	friend void IProduct
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ******************************< ITProduct >*********************************
	// ============================================================================
		// ITProduct(R,x,&y);
	friend void ITProduct
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// ITProduct(R,&x);
	friend void ITProduct
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ********************************< Division >********************************
	// ============================================================================
		// A = B/3.;
	friend BzzMatrixRight operator /
		(const BzzMatrixRight& lval, double rval);

	// A /= 3.;
	BzzMatrixRight& operator /=
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
	friend void Delete(BzzMatrixRight* R);
	friend void ChangeDimensions(int rows, int columns,
		BzzMatrixRight* result, char);
	// from formatted file
	friend void Load(BzzMatrixRight* A, char* filematrix);
	// from binary file
	friend void Load(BzzMatrixRight* A, char, char* filematrix);
	friend char Inverse(BzzMatrixRight* R);
	friend void Swap(BzzMatrixRight* lval, BzzMatrixRight* rval);

	// ============================================================================
	// =======================< Functions for linear algebra >=====================
	// ============================================================================
	friend void Solve(const BzzMatrixRight& R, BzzVector* bx);
	friend void Solve(const BzzMatrixRight& R,
		const BzzVector& b, BzzVector* x);
	friend void Solve(const BzzMatrixRight& R, BzzMatrix* BX);
	friend void Solve(const BzzMatrixRight& R,
		const BzzMatrix& B, BzzMatrix* X);
	friend void TransposeSolve(const BzzMatrixRight& R,
		BzzVector* bx);
	friend void TransposeSolve(const BzzMatrixRight& R,
		const BzzVector& b, BzzVector* x);
	friend void TransposeSolve(const BzzMatrixRight& R,
		BzzMatrix* BX);
	friend void TransposeSolve(const BzzMatrixRight& R,
		const BzzMatrix& B, BzzMatrix* X);

	double Determinant(void);
	double ConditionNumber();
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzMatrixRight* result, char zero = 0);


#endif // BZZ_RIGHT_DOUBLE_HPP