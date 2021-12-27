// BZZMATH: Release 7.0

// ===================< BzzMatrixSymmetric.hpp >=========================
// * Class BzzMatrixSymmetric for symmetric matrices							*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 12)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 12)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *																									*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixSymmetric\MatrixSymmetric.cpp						*
// *																									*
// * Eigenvalues problems:																		*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 12,13)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearAlgebra\Eigenvalues\							*
// *				MatrixSymmetric\MatrixSymmetric.cpp						*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
// 04-1991	Date Written
// 11-1992	English version
// 03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
// 04-1993	Added GetEigenvalues and GetEigenvaluesAndEigenvectors functions.
// 09-1994	Added shadow variable for returning object.
// 09-1994	Added functions ObjectCount, ObjectCountInScope.
// 10-1994	Minor changes throughout.
// 11-1994	Conversion to double precision done.
// 02-1995	Added the Bzz prefix to the names of the classes.
// 09-1995	Added TProduct, PtT,IPt,PtI, ITPt, PtIT, IPtT, TPtI
//	03-1999	Added GetRow, GetColumn, GetDiagonal.
// 03-1999	Added GetHandle function.
// 03-1999	Added GetPointer function.
// 05-2011	Added BzzMatrixSymmetricBand class.

// Release 7.0
// 04_2013	Added ConvertToPositiveDefinite function.

// ============================================================================
// ****** Constructions for BzzMatrixSymmetric:									*
// * BzzMatrixSymmetric S; // default												*
// * BzzMatrixSymmetric S = symm; // copy-initializer							*
// * BzzMatrixSymmetric S(3,3); // sized and placed at 0						*
// * BzzMatrixSymmetric S(3,3,														*
// *		1.,																						*
// *		2.,3.,																					*
// *		4.,5.,6.);// sizes and initialises												*
// * float x[]=																					*
// *		{																							*
// *		1.,																						*
// *		2.,3.,																					*
// *		4.,5.,6.																					*
// *		};																							*
// * BzzMatrixSymmetric S(3,3,x); // from array									*
// * BzzMatrixSymmetric S("SYMM.DAT"); // Formatted File						*
// * BzzMatrixSymmetric S('*',SYMM.BIN"); // Binary File						*
// ****************************************************************************
// ***** Access functions :																	*
// *	i = S.Rows(); // numRows																*
// *	i = S.Columns(); // numColumns														*
// *	xf = S.GetValue(i);																		*
// *	xf = S(i,j);																				*
// *	S(i,j) = xf;																				*
// *	xf = S[i][j]; // NOT implemented														*
// *	S[i][j] = xf; // NOT implemented														*
// *	S.SetValue(i,j,xf);																		*
//	*	S.GetRow(i,&v);																			*
//	*	S.GetColumn(j,&v);																		*
//	*	S.GetDiagonal(j,&v); // j = 0 principal positive at right					*
// *	int who = S.WhoAmI();																	*
// *	int count = BzzMatrixSymmetric::ObjectCount();							*
// *	int countInScope = BzzMatrixSymmetric::ObjectCountInScope();		*
// ****************************************************************************
// ***** Assignments:																			*
// *	S = symm; // symm BzzMatrixSymmetric										*
// ****************************************************************************
// *****	Operators for tests:																	*
// * if(S1 == S2)																					*
// * if(S1 != S2)																					*
// ****************************************************************************
// ***** BzzMessage and BzzPrint																*
// * S.BzzPrint("comments");																	*
// * S.BzzMessage("comments");																*
// ****************************************************************************
// ***** Save and Load																			*
// * S.Save("SYMM.DAT");																		*
// * S.Save('*',"SYMM.BIN");																	*
// * Load(&S,"SYMM.DAT");																		*
// * Load(&S,'*',"SYMM.BIB");																	*
// ****************************************************************************
// ***** Delete, ChangeDimensions and Swap												*
// * Delete(&S); // eliminates BzzMatrixSymmetric								*
// * ChangeDimensions(rows,columns,&S);													*
// * Swap(&S1,&S2);																				*
// ****************************************************************************
// ***** Implemented operations :															*
//	* Sum(A,S,&B);				// B = A + S													*
//	* Sum(&A,S);				// A = A + S													*
//	* Sum(S,A,&B);				// B = S + A													*
//	* Sum(S,&A);				// A = S + A													*

//	* Difference(A,S,&B);	// B = A - S													*
//	* Difference(&A,S);		// A = A - S													*
//	* Difference(S,A,&B);	// B = S - A													*
//	* Difference(S,&A);		// A = S - A													*

//	* Product(A,S,&B);		// B = AS														*
//	* Product(&A,S);			// A = AS														*
//	* Product(S,A,&B);		// B = SA														*
//	* Product(S,&A);			// A = SA														*

//	* TProduct(A,S,&B);		// B = ATS														*
//	* TProduct(&A,S);			// A = ATS														*
//	* TProduct(S,A,&B);		// B = STA														*
//	* TProduct(S,&A);			// A = STA														*

//	* ProductT(A,S,&B);		// B = AST														*
//	* ProductT(&A,S);			// A = AST														*
//	* ProductT(S,A,&B);		// B = SAT														*
//	* ProductT(S,&A);			// A = SAT														*

//	* IProduct(A,S,&B);		// B = A-1S														*
//	* IProduct(&A,S);			// A = A-1S														*
//	* IProduct(S,A,&B);		// B = S-1A														*
//	* IProduct(S,&A);			// A = S-1A														*

//	* ProductI(A,S,&B);		// B = AS-1														*
//	* ProductI(&A,S);			// A = AS-1														*
//	* ProductI(S,A,&B);		// B = SA-1														*
//	* ProductI(S,&A);			// A = SA-1														*

//	* ITProduct(A,S,&B);		// B = A-TS														*
//	* ITProduct(&A,S);		// A = A-TS														*
//	* ITProduct(S,A,&B);		// B = S-TA														*
//	* ITProduct(S,&A);		// A = S-TA														*

//	* ProductIT(A,S,&B);		// B = AS-T														*
//	* ProductIT(&A,S);		// A = AS-T														*
//	* ProductIT(S,A,&B);		// B = SA-T														*
//	* ProductIT(S,&A);		// A = SA-T														*

//	* IProductT(A,S,&B);		// B = A-1ST													*
//	* IProductT(&A,S);		// A = A-1ST													*
//	* IProductT(S,A,&B);		// B = S-1AT													*
//	* IProductT(S,&A);		// A = S-1AT													*

//	* TProductI(A,S,&B);		// B = ATS-1													*
//	* TProductI(&A,S);		// A = ATS-1													*
//	* TProductI(S,A,&B);		// B = STA-1													*
//	* TProductI(S,&A);		// A = STA-1													*

// * S3 = S1 + S2;			 // S3 = S1 + S2;												*
// * S3 = S1 - S2;			 // S3 = S1 - S2;												*
// * y = S1*x;																						*
// * Product(S,x,&y);																			*
// * Product(S,&x);																				*
// * S2 = 3.*S1;																					*
// * S1 *= 3.;					// S1 = 3.*S1;													*
// * L.SelfTProduct(&S);	// S = LTL														*
// * L.SelfProductT(&S);	// S = LLT														*
// * R.SelfTProduct(&S);	// S = RTR														*
// * R.SelfProductT(&S);	// S = RRT														*
// * IProduct(S,x,&y);																			*
// * IProduct(S,&x);																				*
// * S *= .01;																						*
// * double *matrix = A.GetHandle();														*
// * double **matrix = A.GetPointer();														*
// ****************************************************************************
// ***** Special functions for eigenvalues problems:									*
// ***** Eigenvalues and eigenvectors for the symmetric matrix S:					*
// * S.GetEigenvalues(&d); // Only Eigenvalues											*
// * S.GetEigenvaluesAndEigenvectors(&d,&U); // Eigenvalues and Eigenvectors	*
// * The columns of U contain the eigenvectors of matrix S.							*
// ***** Eigenvalues and eigenvectors for rank two matrices S = avvT + bwwT	*
// * Input: in BzzMatrix U(n,2) the vectors v and w, in BzzVector d(2)	*
// * the constants a and b.																	*
// * RankTwoEigenvaluesAndEigenvectors(&d,&U);											*
// * Output: the columns of U contain the eigenvectors of matrix S.				*
// ****************************************************************************

#ifndef BZZ_SYMM_DOUBLE_HPP
#define BZZ_SYMM_DOUBLE_HPP

// preventive declarations
class BzzVector;
class BzzMatrix;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;
class BzzFactorizedSymmetricBase;

// ============================================================================
// =========================< class BzzMatrixSymmetric >=================
// ============================================================================
class BzzMatrixSymmetric : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrix;
	friend class BzzMatrixRight;
	friend class BzzMatrixLeft;
	friend class BzzFactorizedSymmetricBase;
	friend void RankTwoEigenvaluesAndEigenvectors(BzzVector* d,
		BzzMatrix* U);
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

	void Initialize(int rows, int columns);
	// private constructor R('*',3,3)

	// deinitialisation
	void Deinitialize(void);

	BzzMatrixSymmetric(char, int rows, int columns); //OK

	// for eigenvalue problems
	void SymmTridiagonal(BzzMatrix& U, BzzVector& d,
		BzzVector& z, int iautovect);
	void SymmLQ(BzzMatrix& U, BzzVector& d, BzzVector& z,
		int iautovect);

public:
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft S
	BzzMatrixSymmetric(void);

	// copy-initializer // BzzMatrixSymmetric S = symm;
	BzzMatrixSymmetric(BzzMatrixSymmetric& rval);

	// sized and initialised to 0
	// BzzMatrixSymmetric S(3,3)
	BzzMatrixSymmetric(int rows, int columns);

	// sized and initialised
	// BzzMatrixSymmetric S(2,2,1.,2.,3.);
	BzzMatrixSymmetric(int rows, int columnsn, double a11, ...);

	// from array // BzzMatrixSymmetric S(3,3,w)
	BzzMatrixSymmetric(int rows, int columns, double* initvalues);

	// from formatted File // BzzMatrixSymmetric S("SYMM.DAT");
	BzzMatrixSymmetric(char* filematrix);

	// from binary File // BzzMatrixSymmetric L('*',"SYMM.BIN");
	BzzMatrixSymmetric(char, char* filematrix);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixSymmetric(void);

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

	// receives the values with control
	double GetValue(int row, int col) const;
	//takes the row i of a matrix
	void GetRow(int i, BzzVector* v);

	//takes the column j of a matrix
	void GetColumn(int j, BzzVector* v);

	// takes the diagonal i of a squared matrix
	// i = 0 principal
	// i > 0 right diagonal
	// i < 0 left diagonal
	void GetDiagonal(int i, BzzVector* v);

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// assigns values with control
	void SetValue(int row, int col, double val);

	// assigns and receives vector values with control
	double& operator () (int row, int col);

	// ============================================================================
	// **************************< assignment operators >**************************
	// ============================================================================
	BzzMatrixSymmetric& operator =
		(const BzzMatrixSymmetric& rval);

	// ============================================================================
	// **************************< operators for tests >***************************
	// ============================================================================
	friend char operator ==
		(const BzzMatrixSymmetric& lval,
			const BzzMatrixSymmetric& rval);

	friend char operator !=
		(const BzzMatrixSymmetric& lval,
			const BzzMatrixSymmetric& rval);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================
		// Sum(A,D,&B);
	friend void Sum
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// Sum(&A,D);
	friend void Sum
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// Sum(D,A,&B);
	friend void Sum
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Sum(D,&A);
	friend void Sum
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*****************************< Difference >*********************************
	//	============================================================================
		// Difference(A,D,&B);
	friend void Difference
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// Difference(&A,D);
	friend void Difference
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// Difference(D,A,&B);
	friend void Difference
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Difference(D,&A);
	friend void Difference
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
		// Product(A,D,&B);
	friend void Product
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// Product(&A,D);
	friend void Product
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// Product(D,A,&B);
	friend void Product
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// Product(D,&A);
	friend void Product
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< TProduct >*********************************
	//	============================================================================
		// TProduct(A,D,&B);
	friend void TProduct
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// TProduct(&A,D);
	friend void TProduct
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// TProduct(D,A,&B);
	friend void TProduct
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// TProduct(D,&A);
	friend void TProduct
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,D,&B);
	friend void ProductT
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// ProductT(&A,D);
	friend void ProductT
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// ProductT(D,A,&B);
	friend void ProductT
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductT(D,&A);
	friend void ProductT
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< IProduct >*********************************
	//	============================================================================
		// IProduct(A,D,&B);
	friend void IProduct
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// IProduct(&A,D);
	friend void IProduct
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// IProduct(D,A,&B);
	friend void IProduct
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// IProduct(D,&A);
	friend void IProduct
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	*******************************< ProductI >*********************************
	//	============================================================================
		// ProductI(A,D,&B);
	friend void ProductI
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// ProductI(&A,D);
	friend void ProductI
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// ProductI(D,A,&B);
	friend void ProductI
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductI(D,&A);
	friend void ProductI
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< ITProduct >*********************************
	//	============================================================================
		// ITProduct(A,D,&B);
	friend void ITProduct
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// ITProduct(&A,D);
	friend void ITProduct
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// ITProduct(D,A,&B);
	friend void ITProduct
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ITProduct(D,&A);
	friend void ITProduct
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< ProductIT >*********************************
	//	============================================================================
		// ProductIT(A,D,&B);
	friend void ProductIT
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// ProductIT(&A,D);
	friend void ProductIT
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// ProductIT(D,A,&B);
	friend void ProductIT
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// ProductIT(D,&A);
	friend void ProductIT
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< IProductT >*********************************
	//	============================================================================
		// IProductT(A,D,&B);
	friend void IProductT
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// IProductT(&A,D);
	friend void IProductT
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// IProductT(D,A,&B);
	friend void IProductT
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// IProductT(D,&A);
	friend void IProductT
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< TProductI >*********************************
	//	============================================================================
		// TProductI(A,D,&B);
	friend void TProductI
	(BzzMatrix& lval, BzzMatrixSymmetric& rval,
		BzzMatrix* result);

	// TProductI(&A,D);
	friend void TProductI
	(BzzMatrix* lvalAndResult, BzzMatrixSymmetric& rval);

	// TProductI(D,A,&B);
	friend void TProductI
	(BzzMatrixSymmetric& lval, BzzMatrix& rval,
		BzzMatrix* result);

	// TProductI(D,&A);
	friend void TProductI
	(BzzMatrixSymmetric& rval, BzzMatrix* rvalAndResult);

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Sum >*************************************
	// ============================================================================
	friend BzzMatrixSymmetric operator +
		(const BzzMatrixSymmetric& lval, const BzzMatrixSymmetric& rval);

	// ============================================================================
	// *****************************< Difference >*********************************
	// ============================================================================
	friend BzzMatrixSymmetric operator -
		(const BzzMatrixSymmetric& lval, const BzzMatrixSymmetric& rval);

	// ============================================================================
	// ********************************< Product >*********************************
	// ============================================================================
		// S *= 3.;
	BzzMatrixSymmetric& operator *= (double lval);

	// S2 = 3.*S1;
	friend BzzMatrixSymmetric operator *
		(double lval, const BzzMatrixSymmetric& rval);

	// y = S*x;
	friend BzzVector operator *
		(BzzMatrixSymmetric& lval, const BzzVector& rval);

	// Product(S,x,&y);
	friend void Product
	(BzzMatrixSymmetric& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(S,&x);
	friend void Product
	(BzzMatrixSymmetric& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ******************************< IProduct >**********************************
	// ============================================================================
		// IProduct(S,x,&y);
	friend void IProduct
	(const BzzMatrixSymmetric& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(S,&x);
	friend void IProduct
	(const BzzMatrixSymmetric& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// *********************************< Save >***********************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

// ============================================================================
// ========================< Modifying Functions >=============================
// ============================================================================
	void ConvertToPositiveDefinite(void);
	friend void Delete(BzzMatrixSymmetric* S);
	friend void ChangeDimensions(int rows, int columns,
		BzzMatrixSymmetric* result, char);
	// from formatted file
	friend void Load(BzzMatrixSymmetric* A, char* filematrix);
	// from binary file
	friend void Load(BzzMatrixSymmetric* A, char, char* filematrix);
	friend void Swap(BzzMatrixSymmetric* lval,
		BzzMatrixSymmetric* rval);

	// transforms a BzzMatrixSymmetric in BzzFactorizedSymmetricBase
	// BzzMatrixSymmetric gets destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSymmetric* lval, BzzFactorizedSymmetricBase* rval);
	double* GetHandle(void) { return matrix[0] + 1; }
	double** GetPointer(void) { return matrix; }

	// ============================================================================
	// ===============< Functions for eigenvalue problems >========================
	// ============================================================================
	void GetEigenvalues(BzzVector* d);
	void GetEigenvaluesAndEigenvectors(BzzVector* d, BzzMatrix* V);
};

// ===================< BzzMatrixSymmetricBand.hpp >===========================
// * Class BzzMatrixSymmetricBand for symmetric band matrices						*
// ============================================================================
class BzzFactorizedSymmetricBand;
class BzzMatrixSymmetricBand : public BzzBaseClass
{
	friend class BzzFactorizedSymmetricBand;
	friend void ReplaceBzzMatrixWithBzzFactorized(BzzMatrixSymmetricBand* M,
		BzzFactorizedSymmetricBand* F);
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	double** matrix;
	int band, numRows, numColumns;
	int size;
	int whoAmI;
	char shadow;

	void Initialize(int rows, int columns, int b);
	void Deinitialize(void);

	// private constructor R('*',3,3)
public:
	// ============================================================================
	// ****************************< constructors >********************************
	// ============================================================================
		// default // BzzMatrixLeft S
	BzzMatrixSymmetricBand(void);
	void operator()(int r, int c, int b);

	// copy-initializer // BzzMatrixSymmetricBand S = symm;
	BzzMatrixSymmetricBand(BzzMatrixSymmetricBand& rval);

	// sized and initialised to 0
	// BzzMatrixSymmetricBand S(3,3,1)
	BzzMatrixSymmetricBand(int rows, int columns, int b);

	// sized and initialised
	// BzzMatrixSymmetricBand S(2,2,1,1.,2.,3.);
	BzzMatrixSymmetricBand(int rows, int columns, int b, double a11, ...);

	// from formatted File // BzzMatrixSymmetricBand S("SYMM.DAT");
	BzzMatrixSymmetricBand(char* filematrix);

	// from binary File // BzzMatrixSymmetricBand L('*',"SYMM.BIN");
	BzzMatrixSymmetricBand(char, char* filematrix);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixSymmetricBand(void);

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

	// band
	int Band(void) const
	{
		return band;
	}

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// assigns and receives vector values with control
	double& operator () (int row, int col);
	// assigns and receives vector values without control
	double* operator [] (int r)
	{
		return matrix[r];
	}

	// ============================================================================
	// **************************< assignment operators >**************************
	// ============================================================================
	//	BzzMatrixSymmetricBand &operator =
	//			(const BzzMatrixSymmetricBand &rval);

	// ============================================================================
	// **************************< operators for tests >***************************
	// ============================================================================

	// ============================================================================
	// =============================< OPERATIONS >=================================
	// ============================================================================

	// ============================================================================
	// ********************************< Product >*********************************
	// ============================================================================
		// S *= 3.;
	BzzMatrixSymmetricBand& operator *= (double lval);

	// S2 = 3.*S1;
	friend BzzMatrixSymmetricBand operator *
		(double lval, const BzzMatrixSymmetricBand& rval);

	// y = S*x;
	friend BzzVector operator *
		(BzzMatrixSymmetricBand& lval, const BzzVector& rval);

	// Product(S,x,&y);
	friend void Product
	(BzzMatrixSymmetricBand& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(S,&x);
	friend void Product
	(BzzMatrixSymmetricBand& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// =====================< Non-modifying functions >============================
	// ============================================================================

	// ********************************< BzzPrint >********************************
	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// ========================< Modifying Functions >=============================
	// ============================================================================
	friend void Delete(BzzMatrixSymmetricBand* S);
	friend void ChangeDimensions(int b, int rows, int columns,
		BzzMatrixSymmetricBand* result, char);
	// from formatted file
//	friend void Load(BzzMatrixSymmetricBand *A,char *filematrix);
	// from binary file
//	friend void Load(BzzMatrixSymmetricBand *A,char,char *filematrix);
//	friend void Swap(BzzMatrixSymmetricBand *lval,
//		BzzMatrixSymmetricBand *rval);

	// transforms a BzzMatrixSymmetricBand in BzzFactorizedSymmetricBase
	// BzzMatrixSymmetricBand gets destroyed
//	friend void ReplaceBzzMatrixWithBzzFactorized
//	(BzzMatrixSymmetricBand *lval,BzzFactorizedSymmetricBase *rval);
	double* GetHandle(void) { return matrix[0] + 1; }
	double** GetPointer(void) { return matrix; }
	void SetDiagonal(int j, BzzVector& rval);
	void SetDiagonal(int j, double xf);

	// Product(S,x,&y);
	friend void Product
	(BzzMatrixSymmetricBand& S, BzzVector& x,
		BzzVector* y);

	// ============================================================================
	// ===============< Functions for eigenvalue problems >========================
	// ============================================================================
	//	void GetEigenvalues(BzzVector *d);
	//	void GetEigenvaluesAndEigenvectors(BzzVector *d,BzzMatrix *V);
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzMatrixSymmetric* result, char zero = 0);
void ChangeDimensions(int b, int rows, int columns,
	BzzMatrixSymmetricBand* result, char zero = 0);


#endif // BZZ_SYMM_DOUBLE_HPP