// BZZMATH: Release 7.0

//	=========================< BzzMatrixBand >============================
//	* BzzMatrixBand: Class for operations with banded matrices				*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exmaband.cpp							 				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-1995	Date Written.
//	11-1997	Added TProduct(B,x,&y);
//	09-1998	Added LowerBand and UpperBand functions.
//	11-1999	Added Save and Load functions.

////////////////// Release 5.0
//	11-2003	Added SetDimensions function.

////////////////// Release 7.0
//	12-2013	Added operator(row,columns,low,up) = SetDimensions function.

#ifndef BZZ_MATRIX_DOUBLE_BAND_HPP
#define BZZ_MATRIX_DOUBLE_BAND_HPP

class BzzMatrixSparse;
class BzzVectorInt;

//	============================================================================
//	======================< class BzzMatrixBand >=========================
//	============================================================================

class BzzMatrixBand : public BzzBaseClass
{
	friend class BzzFactorizedBandGauss;
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzDaeSparse;
	friend class BzzOdeSparseStiff;
	friend class BzzDaeSparseObject;
	friend class BzzOdeSparseStiffObject;

private:
	static const char* const BZZ_ERROR;
	static int count; // per whoAmI
	int numRows, numColumns;
	int size;
	int	lowerBand,
		upperBand,
		totalBand;

	double** matrix;

	int whoAmI;

	// initialising constructors
	void Initialize(int m, int n, int low, int up);

	// re-initialising
	void ReInitialize(int m, int n, int low, int up);
	friend void ChangeDimensions(int m, int n, int low, int up,
		BzzMatrixBand* B);

	// deinitialisation
	void Deinitialize(void);

	// preparing assignments
	void PrepCopy(int rRows, int rColumns);

	// calculating the norm
	void Norm(void);

	void Copy(BzzMatrixSparse& rval);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixBand(void);

	// copy constructor
	BzzMatrixBand(const BzzMatrixBand& rval);

	// constructor from BzzMatrixSparse
	BzzMatrixBand(BzzMatrixSparse& rval);

	// sizing constructor
	BzzMatrixBand(int rows, int columns, int low, int up);
	void SetDimensions(int rows, int columns, int low, int up);
	void operator()(int rows, int columns, int low, int up);

	// initialises from formatted File
// BzzMatrixBand A("MAT.DAT");
	BzzMatrixBand(char* filematrix);

	// initialises from binary File
	// BzzMatrixBand A('*',"MAT.BIN"); See Save
	BzzMatrixBand(char, char* filematrix);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzMatrixBand(void);

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

	// assigns and receives vector values with control
	double& operator ()
		(int row, int column);

	// assigns and receives vector values without control
	double* operator [] (int r)
	{
		return matrix[r];
	}

	int LowerBand(void) { return lowerBand; }
	int UpperBand(void) { return upperBand; }
	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	BzzMatrixBand& operator =
		(BzzMatrixBand& rval);

	BzzMatrixBand& operator =
		(BzzMatrixSparse& rval);

	// transforms a BzzMatrixBand in BzzFactorizedBandGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixBand* lval, BzzFactorizedBandGauss* rval);

	// transforms a BzzMatrixBand in BzzFactorizedBandGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzMatrix
	(BzzMatrixBand* lval, BzzFactorizedBandGauss* rval);

	/*
	//	============================================================================
	//	==============================< OPERATIONS >================================
	//	============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================

		friend BzzMatrixBand operator +
			 (const BzzMatrixBand &lval,const BzzMatrixBand &rval);

		BzzMatrixBand &operator +=
			 (const BzzMatrixBand &rval);

	//	============================================================================
	//	****************************< Difference >**********************************
	//	============================================================================

		friend BzzMatrixBand operator -
			 (const BzzMatrixBand &lval,const BzzMatrixBand &rval);

		BzzMatrixBand &operator -=
			 (const BzzMatrixBand &rval);

	*/
	//	============================================================================
	//	*******************************< Product >**********************************
	//	============================================================================
		// Product(c,B,&C);
	friend void Product
	(double lval, BzzMatrixBand& rval, BzzMatrixBand* result);

	// Product(c,&B);
	friend void Product
	(double lval, BzzMatrixBand* rvalAndResult);

	// Product(A,B,&C);
	friend void Product
	(BzzMatrix& lval, BzzMatrixBand& rval, BzzMatrix* result);

	// Product(&A,B);
//	friend void Product
//		(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

	// Product(B,A,&C);
	friend void Product
	(BzzMatrixBand& lval, BzzMatrix& rval, BzzMatrix* result);

	// Product(B,&A);
//	friend void Product
//		(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

//	friend BzzVector operator *
//	 (BzzMatrixBand &lval,BzzVector &rval);

	// Product(B,x,&y);
	friend void Product
	(BzzMatrixBand& lval, BzzVector& rval, BzzVector* result);

	// Product(B,&x);
	friend void Product
	(BzzMatrixBand& rval, BzzVector* rvalAndResult);

	// Product(B1,B2,&B3);
	friend void Product
	(BzzMatrixBand& lval, BzzMatrixBand& rval,
		BzzMatrixBand* result);

	// Product(&B1,B2);
	friend void Product
	(BzzMatrixBand* lvalAndResult, BzzMatrixBand& rval);

	// Product(B1,&B2);
	friend void Product
	(BzzMatrixBand& lval, BzzMatrixBand* rvalAndResult);

	// Product(&B);
	friend void Product
	(BzzMatrixBand* lvalAndRvalAndResult);

	//	============================================================================
	//	********************************< TProduct >*********************************
	//	============================================================================
		// TProduct(A,B,&C);
	//	friend void TProduct
	//		(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// TProduct(&A,B);
	//	friend void TProduct
	//		(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// TProduct(B,A,&C);
	//	friend void TProduct
	//		(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// TProduct(B,&A);
	//	friend void TProduct
	//		(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

		// TProduct(B,x,&y);
	friend void TProduct
	(BzzMatrixBand& B, BzzVector& x, BzzVector* y);

	// TProduct(B,&y);
	friend void TProduct
	(BzzMatrixBand& B, BzzVector* y);

	/*
	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,B,&C);
		friend void ProductT
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// ProductT(&A,B);
		friend void ProductT
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// ProductT(B,A,&C);
		friend void ProductT
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// ProductT(B,&A);
		friend void ProductT
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

	//	============================================================================
	//	*******************************< IProduct >*********************************
	//	============================================================================
		// IProduct(A,B,&C);
		friend void IProduct
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// IProduct(&A,B);
		friend void IProduct
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// IProduct(B,A,&C);
		friend void IProduct
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// IProduct(B,&A);
		friend void IProduct
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

		// IProduct(B,x,&y);
		friend void IProduct
			(BzzMatrixBand &lval,BzzVector &rval,BzzVector *result);

		// IProduct(B,&x);
		friend void IProduct
			(BzzMatrixBand &rval,BzzVector *rvalAndResult);

	//	============================================================================
	//	********************************< ProductI >********************************
	//	============================================================================
		// ProductI(A,B,&C);
		friend void ProductI
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// ProductI(&A,B);
		friend void ProductI
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// ProductI(B,A,&C);
		friend void ProductI
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// ProductI(B,&A);
		friend void ProductI
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

	//	============================================================================
	//	********************************< ITProduct >*******************************
	//	============================================================================
		// ITProduct(A,B,&C);
		friend void ITProduct
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// ITProduct(&A,B);
		friend void ITProduct
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// ITProduct(B,A,&C);
		friend void ITProduct
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// ITProduct(B,&A);
		friend void ITProduct
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

	//	============================================================================
	//	********************************< ProductIT >*******************************
	//	============================================================================
		// ProductIT(A,B,&C);
		friend void ProductIT
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// ProductIT(&A,B);
		friend void ProductIT
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// ProductIT(B,A,&C);
		friend void ProductIT
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// ProductIT(B,&A);
		friend void ProductIT
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

	//	============================================================================
	//	********************************< IProductT >*******************************
	//	============================================================================
		// IProductT(A,B,&C);
		friend void IProductT
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// IProductT(&A,B);
		friend void IProductT
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// IProductT(B,A,&C);
		friend void IProductT
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// IProductT(B,&A);
		friend void IProductT
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);

	//	============================================================================
	//	********************************< TProductI >*******************************
	//	============================================================================
		// TProductI(A,B,&C);
		friend void TProductI
			(BzzMatrix &lval,BzzMatrixBand &rval,BzzMatrix *result);

		// TProductI(&A,B);
		friend void TProductI
			(BzzMatrix *lvalAndResult,BzzMatrixBand &rval);

		// TProductI(B,A,&C);
		friend void TProductI
			(BzzMatrixBand &lval,BzzMatrix &rval,BzzMatrix *result);

		// TProductI(B,&A);
		friend void TProductI
			(BzzMatrixBand &rval,BzzMatrix *rvalAndResult);
	*/

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void);

	//	********************************< Save >************************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

//	*****************************< Max and Min >********************************
	//to have the position Max(im,jc)
	double Max(int* imax = 0, int* jmax = 0);
	double MaxAbs(int* imax = 0, int* jmax = 0);

	double Min(int* imin = 0, int* jmin = 0);
	double MinAbs(int* imin = 0, int* jmin = 0);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	friend void Delete(BzzMatrixBand* A);
	void SetToZeroAndReset(void);// eliminates all elements
	friend void Swap(BzzMatrixBand* A, BzzMatrixBand* B);
	void SetDiagonal(int j, BzzVector& rval);
	void SetDiagonal(int j, float xf);
	void SetDiagonal(int j, double xf);
	// recovery of formatted Save
	friend void Load(BzzMatrixBand* A, char* filemat);
	// recovery of binary Save
	friend void Load(BzzMatrixBand* A, char, char* filemat);
};

#endif // BZZ_MATRIX_DOUBLE_BAND_HPP