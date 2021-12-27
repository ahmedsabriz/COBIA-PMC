// BZZMATH: Release 7.0

//	======================< BzzMatrixDiagonal.hpp >=======================
//	* Class BzzMatrixDiagonal for Diagonal square matrices					*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 3, 6)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixDiagonal\MatrixDiagonal.cpp						*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1995	Date Written

////////////////// Release 5.0

//	============================================================================
//	****** BzzMatrixDiagonal constructors:											*
//	* BzzMatrixDiagonal D; // default												*
//	* BzzMatrixDiagonal D2 = D1; // copy-initializer							*
//	* BzzMatrixDiagonal D(3); // sizes and places at 0							*
//	* BzzMatrixDiagonal D(3,1.,2.,3.);	// sizes and initialises			*
//	* BzzMatrixDiagonal D = v; // conversion from BzzVector			*
//	* double x[]={1.,2.,3.};																	*
//	* BzzMatrixDiagonal D(3,x); // from array										*
//	* BzzMatrixDiagonal D("DIAGONAL.DAT"); // Formatted File					*
//	* BzzMatrixDiagonal D('*',DIAGONAL.BIN"); // Binary File					*
//	****************************************************************************
//	***** Access functions :																	*
//	*	i = D.Rows(); // numRows																*
//	*	i = D.Columns(); // numColumns														*
//	*	xf = D.GetValue(i);																		*
//	*	xf = D(i);																					*
//	*	D(i) = xf;																					*
//	*	xf = D[i];																					*
//	*	D[i] = xf;																					*
//	*	D.SetValue(i,xf);																			*
//	*	int who = D.WhoAmI();																	*
//	*	int count = BzzMatrixDiagonal::ObjectCount();							*
//	*	int countInScope = BzzMatrixDiagonal::ObjectCountInScope();			*
//	****************************************************************************
//	***** Assignments:																			*
//	*	D = c; // c double																		*
//	*	D = d; // d BzzVector															*
//	*	D = Diag; // Diag BzzMatrixDiagonal											*
//	****************************************************************************
//	*****	Operators for tests:																	*
//	* if(D1 == D2)																					*
//	* if(D1 != D2)																					*
//	****************************************************************************
//	***** BzzMessage and BzzPrint																*
//	* D.BzzPrint("Comment"); // always														*
//	* D.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	****************************************************************************
//	***** Norms																						*
//	* xf = D.Norm2();																				*
//	* xf = D.NormT();																				*
//	* xf = D.NormR();																				*
//	* xf = D.NormC();																				*
//	* xf = D.NormF();																				*
//	* xf = D.NormI();																				*
//	* xf = D.Norm1();																				*
//	****************************************************************************
//	***** Max and Min																				*
//	* xf = D.Max(&imax);																			*
//	* xf = D.MaxAbs(&imax);																		*
//	* xf = D.Min(&imin);																			*
//	* xf = D.MinAbs(&imin);																		*
//	****************************************************************************
//	***** Save and Load																			*
//	* D.Save("DIAGONAL.DAT");	// formatted												*
//	* D.Save('*',"DIAGONAL.DAT");	// unformatted											*
//	* Load(&D,"DIAGONAL.DAT");																	*
//	* Load(&D,'*',"DIAGONAL.DAT");															*
//	****************************************************************************
//	***** Delete, ChangeDimensions and Swap												*
//	* Delete(&D);																					*
//	* ChangeDimensions(newrc,&D);																*
//	* Swap(&D1,&D2);																				*
//	****************************************************************************
//	***** Implemented operations :															*
//	* Sum(A,D,&B);				// B = A + D													*
//	* Sum(&A,D);				// A = A + D													*
//	* Sum(D,A,&B);				// B = D + A													*
//	* Sum(D,&A);				// A = D + A													*

//	* Difference(A,D,&B);	// B = A - D													*
//	* Difference(&A,D);		// A = A - D													*
//	* Difference(D,A,&B);	// B = D - A													*
//	* Difference(D,&A);		// A = D - A													*

//	* Product(A,D,&B);		// B = AD														*
//	* Product(&A,D);			// A = AD														*
//	* Product(D,A,&B);		// B = DA														*
//	* Product(D,&A);			// A = DA														*

//	* TProduct(A,D,&B);		// B = ATD														*
//	* TProduct(&A,D);			// A = ATD														*
//	* TProduct(D,A,&B);		// B = DTA														*
//	* TProduct(D,&A);			// A = DTA														*

//	* ProductT(A,D,&B);		// B = ADT														*
//	* ProductT(&A,D);			// A = ADT														*
//	* ProductT(D,A,&B);		// B = DAT														*
//	* ProductT(D,&A);			// A = DAT														*

//	* IProduct(A,D,&B);		// B = A-1D														*
//	* IProduct(&A,D);			// A = A-1D														*
//	* IProduct(D,A,&B);		// B = D-1A														*
//	* IProduct(D,&A);			// A = D-1A														*

//	* ProductI(A,D,&B);		// B = AD-1														*
//	* ProductI(&A,D);			// A = AD-1														*
//	* ProductI(D,A,&B);		// B = DA-1														*
//	* ProductI(D,&A);			// A = DA-1														*

//	* ITProduct(A,D,&B);		// B = A-TD														*
//	* ITProduct(&A,D);		// A = A-TD														*
//	* ITProduct(D,A,&B);		// B = D-TA														*
//	* ITProduct(D,&A);		// A = D-TA														*

//	* ProductIT(A,D,&B);		// B = AD-T														*
//	* ProductIT(&A,D);		// A = AD-T														*
//	* ProductIT(D,A,&B);		// B = DA-T														*
//	* ProductIT(D,&A);		// A = DA-T														*

//	* IProductT(A,D,&B);		// B = A-1DT													*
//	* IProductT(&A,D);		// A = A-1DT													*
//	* IProductT(D,A,&B);		// B = D-1AT													*
//	* IProductT(D,&A);		// A = D-1AT													*

//	* TProductI(A,D,&B);		// B = ATD-1													*
//	* TProductI(&A,D);		// A = ATD-1													*
//	* TProductI(D,A,&B);		// B = DTA-1													*
//	* TProductI(D,&A);		// A = DTA-1													*

//	* y = D*x;					// y = Dx														*
//	* Product(D,x,&y);		// y = Dx														*
//	* Product(D,&x);			// x = Dx														*

//	* Product(D1,D2,&D3);	// D3 = D1*D2													*
//	* Product(&D1,D2);		// D1 = D1*D2													*
//	* Product(D1,&D2);		// D2 = D1*D2													*
//	* Product(&D);				// D = D*D														*
//	* IProduct(D,x,&y);		// y = D-1x														*
//	* IProduct(D,&x);			// x = D-1x														*
//	****************************************************************************
//	***** Functions for linear systems:														*
//	* Solve(D,b,&x);																				*
//	* Solve(D,&bx);	else Solve(D,b,&b);													*
//	* Solve(D,B,&X);																				*
//	* Solve(D,&BX);	else Solve(D,B,&B);													*
//	* D.Determinant();																			*
//	* D.ConditionNumber();																		*
//	* Inverse(&D);																					*
//	****************************************************************************

#ifndef BZZ_DIAGONAL_DOUBLE_HPP
#define BZZ_DIAGONAL_DOUBLE_HPP

// preventive declarations
class BzzVector;
class BzzMatrix;

//	============================================================================
//	======================< class BzzMatrixDiagonal >=====================
//	============================================================================
class BzzMatrixDiagonal : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrix;
	friend class BzzSave;
	friend class BzzLoad;

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	double* matrix;
	int numRows;
	int size;
	int whoAmI;
	char shadow;

	// initialisation constructors
	void Initialize(int rows);

	// deinitialisation
	void Deinitialize(void);

	// BzzMatrixDiagonal A('*',3);
	BzzMatrixDiagonal(char, int rows);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default // BzzMatrixDiagonal D
	BzzMatrixDiagonal(void);

	// copy-initializer // BzzMatrixDiagonal D = Diagonal;
	BzzMatrixDiagonal(BzzMatrixDiagonal& rval);

	// sizes and initialises at 0
	// BzzMatrixDiagonal D(3);
	BzzMatrixDiagonal(int rows);

	// sizes and initialises
	// BzzMatrixDiagonal D(3,1.,2.,3.);
	BzzMatrixDiagonal(int rows, double a11, ...);

	// from BzzVector D = v;
	BzzMatrixDiagonal(BzzVector& v);

	// from array // BzzMatrixDiagonal D(3,w)
	BzzMatrixDiagonal(int rows, double* initvalues);

	// from formatted File // BzzMatrixDiagonal D("DIAGONAL.DAT");
	BzzMatrixDiagonal(char* filematrix);

	// from binary File // BzzMatrixDiagonal D('*',"DIAGONAL.BIN");
	BzzMatrixDiagonal(char, char* filematrix);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMatrixDiagonal(void);

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================
		// row number
	int Rows(void) const
	{
		return numRows;
	}

	// column number
	int Columns(void) const
	{
		return numRows;
	}

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// receives the values with control
	double GetValue(int row) const;

	//	============================================================================
	//	**********************< Modifying access functions >************************
	//	============================================================================
		// assigns values with control
	void SetValue(int row, double val);

	// assigns and receives vector values with control
	double& operator () (int row);

	// assigns and receives vector values without control
	double& operator [](int i)
	{
		return matrix[i];
	}

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	BzzMatrixDiagonal& operator =
		(double d);

	BzzMatrixDiagonal& operator =
		(BzzVector& d);

	BzzMatrixDiagonal& operator =
		(const BzzMatrixDiagonal& rval);

	//	============================================================================
	//	**************************< operators for tests >***************************
	//	============================================================================
	friend char operator ==
		(const BzzMatrixDiagonal& lval,
			const BzzMatrixDiagonal& rval);

	friend char operator !=
		(const BzzMatrixDiagonal& lval,
			const BzzMatrixDiagonal& rval);

	//	============================================================================
	//	=============================< OPERATIONS >=================================
	//	============================================================================
	//	============================================================================
	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================
		// Sum(A,D,&B);
	friend void Sum
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// Sum(&A,D);
	friend void Sum
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// Sum(D,A,&B);
	friend void Sum
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Sum(D,&A);
	friend void Sum
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*****************************< Difference >*********************************
	//	============================================================================
		// Difference(A,D,&B);
	friend void Difference
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// Difference(&A,D);
	friend void Difference
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// Difference(D,A,&B);
	friend void Difference
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Difference(D,&A);
	friend void Difference
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
		// Product(A,D,&B);
	friend void Product
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// Product(&A,D);
	friend void Product
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// Product(D,A,&B);
	friend void Product
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Product(D,&A);
	friend void Product
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	friend BzzVector operator *
		(BzzMatrixDiagonal& lval, BzzVector& rval);

	// Product(D,x,&y);
	friend void Product
	(BzzMatrixDiagonal& lval, BzzVector& rval,
		BzzVector* result);

	// Product(D,&x);
	friend void Product
	(BzzMatrixDiagonal& rval, BzzVector* rvalAndResult);

	// Product(D1,D2,&D3);
	friend void Product
	(BzzMatrixDiagonal& lval, BzzMatrixDiagonal& rval,
		BzzMatrixDiagonal* result);

	// Product(&D1,D2);
	friend void Product
	(BzzMatrixDiagonal* lvalAndResult, BzzMatrixDiagonal& rval);

	// Product(D1,&D2);
	friend void Product
	(BzzMatrixDiagonal& lval, BzzMatrixDiagonal* rvalAndResult);

	// Product(&D);
	friend void Product
	(BzzMatrixDiagonal* lvalAndRvalAndResult);

	//	============================================================================
	//	********************************< TProduct >********************************
	//	============================================================================
		// TProduct(A,D,&B);
	friend void TProduct
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// TProduct(&A,D);
	friend void TProduct
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// TProduct(D,A,&B);
	friend void TProduct
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// TProduct(D,&A);
	friend void TProduct
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,D,&B);
	friend void ProductT
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// ProductT(&A,D);
	friend void ProductT
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// ProductT(D,A,&B);
	friend void ProductT
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductT(D,&A);
	friend void ProductT
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< IProduct >*********************************
	//	============================================================================
		// IProduct(A,D,&B);
	friend void IProduct
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// IProduct(&A,D);
	friend void IProduct
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// IProduct(D,A,&B);
	friend void IProduct
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// IProduct(D,&A);
	friend void IProduct
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	// IProduct(D,x,&y);
	friend void IProduct
	(BzzMatrixDiagonal& lval, BzzVector& rval,
		BzzVector* result);

	// IProduct(D,&x);
	friend void IProduct
	(BzzMatrixDiagonal& rval, BzzVector* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductI >*********************************
	//	============================================================================
		// ProductI(A,D,&B);
	friend void ProductI
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// ProductI(&A,D);
	friend void ProductI
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// ProductI(D,A,&B);
	friend void ProductI
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductI(D,&A);
	friend void ProductI
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< ITProduct >*********************************
	//	============================================================================
		// ITProduct(A,D,&B);
	friend void ITProduct
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// ITProduct(&A,D);
	friend void ITProduct
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// ITProduct(D,A,&B);
	friend void ITProduct
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ITProduct(D,&A);
	friend void ITProduct
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< ProductIT >*********************************
	//	============================================================================
		// ProductIT(A,D,&B);
	friend void ProductIT
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// ProductIT(&A,D);
	friend void ProductIT
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// ProductIT(D,A,&B);
	friend void ProductIT
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductIT(D,&A);
	friend void ProductIT
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< IProductT >*********************************
	//	============================================================================
		// IProductT(A,D,&B);
	friend void IProductT
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// IProductT(&A,D);
	friend void IProductT
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// IProductT(D,A,&B);
	friend void IProductT
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// IProductT(D,&A);
	friend void IProductT
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< TProductI >*********************************
	//	============================================================================
		// TProductI(A,D,&B);
	friend void TProductI
	(BzzMatrix& lval, BzzMatrixDiagonal& rval,
		BzzMatrix* result);

	// TProductI(&A,D);
	friend void TProductI
	(BzzMatrix* lvalAndResult, BzzMatrixDiagonal& rval);

	// TProductI(D,A,&B);
	friend void TProductI
	(BzzMatrixDiagonal& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// TProductI(D,&A);
	friend void TProductI
	(BzzMatrixDiagonal& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	=====================< Non-modifying functions >============================
	//	============================================================================

	//	********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	//	*********************************< Save >***********************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

//	******************************< Max and Min >*******************************
	double Max(int* imax = 0);
	double MaxAbs(int* imax = 0);

	double Min(int* imin = 0);
	double MinAbs(int* imin = 0);

	//	********************************< Norms >***********************************
	double Norm2(void);
	double NormT(void);
	double NormR(void);
	double NormC(void);
	double NormF(void);
	double NormI(void) { return NormR(); }
	double Norm1(void) { return NormC(); }

	//	============================================================================
	//	========================< Modifying Functions >=============================
	//	============================================================================
	friend void Delete(BzzMatrixDiagonal* D);
	friend void ChangeDimensions(int rows,
		BzzMatrixDiagonal* result, char);
	// from formatted file
	friend void Load(BzzMatrixDiagonal* D, char* filematrix);
	// from binary file
	friend void Load(BzzMatrixDiagonal* D, char, char* filematrix);
	friend char Inverse(BzzMatrixDiagonal* D);
	friend void Swap(BzzMatrixDiagonal* lval, BzzMatrixDiagonal* rval);

	//	============================================================================
	//	=======================< Functions for linear algebra >=====================
	//	============================================================================
	friend void Solve(BzzMatrixDiagonal& D, BzzVector* bx);
	friend void Solve(BzzMatrixDiagonal& D,
		BzzVector& b, BzzVector* x);
	friend void Solve(BzzMatrixDiagonal& D, BzzMatrix* BX);
	friend void Solve(BzzMatrixDiagonal& D,
		BzzMatrix& B, BzzMatrix* X);

	double Determinant(void);
	double ConditionNumber();
};

// Friend functions with default arguments
void ChangeDimensions(int rows,
	BzzMatrixDiagonal* result, char zero = 0);


#endif // BZZ_DIAGONAL_DOUBLE_HPP