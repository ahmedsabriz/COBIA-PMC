// TODO: Condition Number
// TODO: SPARSE_LQ

// BZZMATH: Release 3.1

//	==============================< FASPAUND.HPP >==============================
//	* BzzFactorizedSparseUnspecified: Class for solution of linear sparse systems*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exfaspun.cpp							 				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	03-1997	Date Written.
//	11-1998	Bug fixed in ReplaceBzzMatrixWithBzzFactorized function.

//	============================================================================
//	******* Functions for solving linear sparse systems:								*
//	* Solve(&A,&bx);																				*
//	****************************************************************************

#ifndef BZZ_FACTORED_DOUBLE_SPARSE_UNSPECIFIED
#define BZZ_FACTORED_DOUBLE_SPARSE_UNSPECIFIED

class BzzMatrixSparse;
struct ElementBzzMatrixSparse;

//	============================================================================
//	==================< class BzzFactorizedSparseUnspecified >====================
//	============================================================================
class BzzFactorizedSparseUnspecified : public BzzBaseClass
{
	friend void Solve
	(BzzFactorizedSparseUnspecified* A, BzzVector* bx);

private:
	enum SparseFactorizationMethod
	{
		UNKNOWN,
		MATRIX_SPARSE,
		FULL_LEFT,
		MATRIX_SPARSE_LEFT,
		FULL_RIGHT,
		MATRIX_SPARSE_RIGHT,
		BAND,
		FULL_GAUSS,
		SPARSE_GAUSS,
		FULL_QR,
		SPARSE_QR,
		FULL_LQ,
		SPARSE_LQ
	}solutionMethod;
	static const char* const BZZ_ERROR;
	static int count; // per whoAmI

	int whoAmI;
	char analysis, ordered;
	int	numRows,
		numColumns,
		lowerBand,
		upperBand,
		band,
		countElements;

	BzzVectorInt	ordRows,
		ordColumns;

	void Initialize(void);
	friend void Delete(BzzFactorizedSparseUnspecified* A);
	void Solution(BzzVector* bx);

	BzzMatrixSparse matrixSparse;
	BzzMatrixLeft matrixLeft;
	BzzMatrixRight matrixRight;
	BzzFactorizedBandGauss factoredBand;
	BzzFactorizedGauss factoredGauss;
	BzzFactorizedSparseGauss factoredSparseGauss;
	BzzFactorizedQR factoredQR;
	//	BzzFactorizedSparseQR factoredSparseQR;
	BzzFactorizedLQ factoredLQ;
	//	BzzFactorizedSparseLQ factoredSparseLQ;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseUnspecified(void);

	// copy constructor
	BzzFactorizedSparseUnspecified(const BzzFactorizedSparseUnspecified& rval);

	// from BzzMatrixSparse
	BzzFactorizedSparseUnspecified(BzzMatrixSparse& A);
	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedSparseUnspecified(void);

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
	int WhoAmI(void) const { return whoAmI; }
	virtual void ObjectBzzPrint(void) {};
	void PrintSelectedMethod(void);
	void Analysis(void);
	double ConditionNumber(void);

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	// from a BzzMatrixSparse
	void operator = (BzzMatrixSparse& A);

	// transforms a BzzMatrixSparse in BzzFactorizedSparseUnspecified
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseUnspecified* rval,
		char);
};

// Friend functions with default arguments
void ReplaceBzzMatrixWithBzzFactorized
(BzzMatrixSparse* lval, BzzFactorizedSparseUnspecified* rval,
	char first = 0);


#endif // BZZ_FACTORED_DOUBLE_SPARSE_UNSPECIFIED