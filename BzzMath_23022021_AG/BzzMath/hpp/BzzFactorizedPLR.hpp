// BZZMATH: Release 7.0

//	========================< BzzFactorizedPLR.hpp >========================
//	* Class BzzFactorizedPLR derived from BzzFactorized						*
//	* Class BzzFactorizedGauss derived from BzzFactorizedPLR				*
//	* Class BzzFactorizedCrout derived from BzzFactorizedPLR				*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 15)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 15)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\LinearSystems\							*
// *				FactorizedGauss\FactorizedGauss.cpp							*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1991	Date Written
//	11-1992	English version
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.
//	06-1996	Improved GaussFactorization efficiency.
//	11-1999	Modified ConditionNumber.

////////////////// Release 4.0
//	07-2000	Added GetLinearCombinations.

////////////////// Release 5.0
//	08-2003	Added SolveWithOrdering.

////////////////// Release 6.0
//	07-2009	BzzFactorizedGauss parallelized wth OpenMP.

////////////////// Release 7.0
//	06-2013	Added GetPivot function in BzzFactorizedGauss class.

//	============================================================================
//	* These classes permit linear system solutions and solutions to related		*
// * problems only with square matrices													*
//	* These are the ones that require the least time for calculation.				*
//	****************************************************************************

#ifndef BZZ_FACTORIZEDPLR_DOUBLE_HPP
#define BZZ_FACTORIZEDPLR_DOUBLE_HPP

//	============================================================================
//	========================< prototype PLR functions >=========================
//	============================================================================
char GaussFactorization
(int n, double** a, int* indx, int* signd, BzzVectorInt& bzzPivot);
char CroutFactorization
(int n, double** a, int* indx, int* signd);
void PLRSolution
(int n, double** a, double* b, int* indx);
void PLRTransposeSolution
(int n, double** a, double* b, int* indx);
double PLRNormInvEst
(int n, double** a, int* indx);

//	============================================================================
//	=========================< class BzzFactorizedPLR >=====================
//	============================================================================
class BzzVector;
class BzzMatrix;
class BzzFactorizedPLR : public BzzFactorized
{
	friend class BzzSave;
	friend class BzzLoad;

protected:
	int* indx;
	int signd;

	// initialisation of constructors
	void FurtherInit(void);

	//	initialisation of special vectors (indx)
	virtual void SpecificInitialize(void);

	//	deinitialisation of special vectors (indx)
	virtual void SpecificDeinitialize(void);

	// constructor A('*',3,5);
	BzzFactorizedPLR(char ch, int rows, int columns)
		: BzzFactorized(ch, rows, columns) {
		FurtherInit();
	}

	//	============================================================================
	//	**************************< Functions for linear algebra >******************
	//	============================================================================
	virtual void	 Factorization(void) = 0;
	virtual void	 Solution(BzzVector* bx);
	virtual void	 Solution(const BzzVector& b, BzzVector* x);
	virtual void	 TransposeSolution(BzzVector* bx);
	virtual void	 TransposeSolution
	(const BzzVector& b, BzzVector* x);
	virtual void	 Solution(BzzMatrix* BX);
	virtual void	 Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void	 TransposeSolution(BzzMatrix* BX);
	virtual void	 TransposeSolution(const BzzMatrix& B, BzzMatrix* X);

	virtual double	Condition(void);
	virtual double	DeterminantEvaluation(void);

public:

	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorizedPLR(void)
		: BzzFactorized() {
		FurtherInit();
	}

	// copy constructor
	BzzFactorizedPLR(const BzzFactorized& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// constructor A(3,3);
	BzzFactorizedPLR(int rows, int columns)
		: BzzFactorized(rows, columns) {
		FurtherInit();
	}

	// constructor A(3,3,w);
	BzzFactorizedPLR(int rows, int columns, double* initvalues)
		: BzzFactorized(rows, columns, initvalues) {
		FurtherInit();
	}

	// constructor from BzzMatrix
	BzzFactorizedPLR(const BzzMatrix& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// makes a submatrix with rows,columns
	BzzFactorizedPLR(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, rval) {
		FurtherInit();
	}

	// as above, commencing from irow,jcol
	BzzFactorizedPLR(int rows, int columns, int irow,
		int jcol, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, irow, jcol, rval) {
		FurtherInit();
	}

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedPLR(void);
	virtual void GetLinearCombinations(BzzVectorInt* bx);
};

//	============================================================================
//	=======================< class BzzFactorizedGauss >=====================
//	============================================================================
class BzzFactorizedGauss : public BzzFactorizedPLR
{
	friend class BzzSave;
	friend class BzzLoad;
	BzzVectorInt ordRows, ordColumns;
	int orderingGauss;
	BzzVectorInt bzzPivot;
public:
public:
	void GetPivot(BzzVectorInt* pivot)
	{
		*pivot = bzzPivot;
	}

	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorizedGauss(void)
		: BzzFactorizedPLR() {}

	// copy constructor
	BzzFactorizedGauss(const BzzFactorized& rval)
		: BzzFactorizedPLR(rval) {}

	// constructor A(3,3);
	BzzFactorizedGauss(int rows, int columns)
		: BzzFactorizedPLR(rows, columns) {}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedGauss(int rows, int columns, double a11, ...);

	// constructor A(3,3,w);
	BzzFactorizedGauss(int rows, int columns, double* initvalues)
		: BzzFactorizedPLR(rows, columns, initvalues) {}

	// constructor from BzzMatrix
	BzzFactorizedGauss(const BzzMatrix& rval)
		: BzzFactorizedPLR(rval) {}

	// makes a submatrix with rows,columns
	BzzFactorizedGauss(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorizedPLR(rows, columns, rval) {}

	// as above, commencing from irow,jcol
	BzzFactorizedGauss(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval)
		: BzzFactorizedPLR(rows, columns, irow, jcol, rval) {}

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedGauss(void);

	//	============================================================================
	//	*****************************< assignment operators >***********************
	//	============================================================================
	void operator = (const BzzFactorizedGauss& rval);
	BzzFactorizedGauss& operator = (const BzzMatrix& rval);

	//	============================================================================
	//	*******************< Functions for linear algebra >*************************
	//	============================================================================
	friend void SolveWithOrdering
	(BzzFactorizedGauss* A, BzzVector* bx);
	virtual void Factorization(void);
	friend void Factorize(BzzFactorizedGauss* A)
	{
		A->PrepOnlySolve();
	}

	//	============================================================================
	//	*********************< Special function for ODE >***************************
	//	============================================================================
	friend void ProductForOde(double ch, BzzFactorizedGauss* A);

	//	============================================================================
	//	*********************< Special function for DAE >***************************
	//	============================================================================
	friend void ProductForOde(double ch, BzzFactorizedGauss* A);
	friend void ProductForDae(double ch, double c, BzzVectorInt& iDer,
		BzzFactorizedGauss* A);
};

//	============================================================================
//	=======================< class BzzFactorizedCrout >=====================
//	============================================================================

class BzzFactorizedCrout : public BzzFactorizedPLR
{
	friend class BzzSave;
	friend class BzzLoad;

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorizedCrout(void)
		: BzzFactorizedPLR() {}

	// copy constructor
	BzzFactorizedCrout(const BzzFactorized& rval)
		: BzzFactorizedPLR(rval) {}

	// constructor A(3,3);
	BzzFactorizedCrout(int rows, int columns)
		: BzzFactorizedPLR(rows, columns) {}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedCrout(int rows, int columns, double a11, ...);

	// constructor A(3,3,w);
	BzzFactorizedCrout(int rows, int columns, double* initvalues)
		: BzzFactorizedPLR(rows, columns, initvalues) {}

	// constructor from BzzMatrix
	BzzFactorizedCrout(const BzzMatrix& rval)
		: BzzFactorizedPLR(rval) {}

	// make a submatrix with rows,columns
	BzzFactorizedCrout(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorizedPLR(rows, columns, rval) {}

	//as above, commencing from irow,jcol
	BzzFactorizedCrout(int rows, int columns, int irow,
		int jcol, const BzzMatrix& rval)
		: BzzFactorizedPLR(rows, columns, irow, jcol, rval) {}

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedCrout(void);

	//	============================================================================
	//	*****************************< assignment operators >***********************
	//	============================================================================
	void operator = (const BzzFactorizedCrout& rval);
	BzzFactorizedCrout& operator = (const BzzMatrix& rval);

	//	============================================================================
	//	*******************< Functions for linear algebra >*************************
	//	============================================================================
	virtual void Factorization(void);
};

#endif // BZZ_FACTORIZEDPLR_DOUBLE_HPP