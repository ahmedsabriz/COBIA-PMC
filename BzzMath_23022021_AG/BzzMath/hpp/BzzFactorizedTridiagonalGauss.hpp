// BZZMATH: Release 7.0

//	================================< TRIDIAGD.HPP >============================
//	* class BzzFactorizedTridiagonalWithoutPivoting and BzzFactorizedTridiagonalGauss	*
//	* These classes are used for solving linear systems								*
//	* with tridiagonal structure																*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\exampled\dxtri.cpp												*
//	============================================================================

// Revision History (MM-YYYY)
//	05-1993	Date Written
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.

//	============================================================================
//	****** BzzFactorizedTridiagonalWithoutPivoting constructors:									*
//	* BzzFactorizedTridiagonalWithoutPivoting A; // default constructor							*
//	* BzzFactorizedTridiagonalWithoutPivoting A(&a,&b,&c); // BzzVector a,b,c					*
//	* BzzFactorizedTridiagonalWithoutPivoting A("TRID.DAT"); // From file ACSII				*
//	****************************************************************************
//	****** BzzFactorizedTridiagonalGauss constructors:							*
//	* BzzFactorizedTridiagonalGauss A; // default constructor					*
//	* BzzFactorizedTridiagonalGauss A(&a,&b,&c); // BzzVector a,b,c			*
//	* BzzFactorizedTridiagonalGauss A("TRID.DAT"); // From file ACSII		*
//	****************************************************************************
//	****** Functions for solving systems:													*
//	* Solve(A,b,&x);																				*
//	* Solve(A,&bx);	or	Solve(A,b,&b);														*
//	* Solve(A,B,&X);																				*
//	* Solve(A,&BX);	or	Solve(A,BX,&BX);													*
//	****************************************************************************
//	* Other functions:																			*
//	* double det = A.Determinant(void)														*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef TRIDIAGONAL_DOUBLE_HPP
#define TRIDIAGONAL_DOUBLE_HPP

//	============================================================================
//	==================< class BzzFactorizedTridiagonalBase >================
//	============================================================================

class BzzFactorizedTridiagonalBase : public BzzBaseClass
{
	//	============================================================================
	//	=====================< Functions for linear algebra >=======================
	//	============================================================================

	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector& b, BzzVector* x);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, const BzzMatrix& B,
		BzzMatrix* X);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzMatrix* BX);

protected:
	enum BzzFactorizedTridiagonalFactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int whoAmI;
	int numRows;
	BzzVector l, d, u;
	double norm;
	int sign;

	void Initialize(void);	// initialising constructors
	void InitializeFile(char* fileBzzFactorizedTridiagonal);
	void PrepSolve(void);

	virtual void Factorization(void) = 0;
	virtual void Solution(BzzVector& b, BzzVector* x) = 0;
	virtual void Solution(BzzVector* bx) = 0;
	virtual void Solution(const BzzMatrix& B,
		BzzMatrix* X) = 0;
	virtual void Solution(BzzMatrix* BX) = 0;
	virtual void operator ()
		(BzzVector* ll, BzzVector* dd,
			BzzVector* uu) = 0;

	//	virtual double ConditionNumber(void) = 0;

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedTridiagonalBase(void);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	virtual ~BzzFactorizedTridiagonalBase(void) {};

	//	============================================================================
	//	===============================< Functions >================================
	//	============================================================================
	virtual void ObjectBzzPrint(void);
	double Determinant(void);
};

//	============================================================================
//	=====================< class BzzFactorizedTridiagonalWithoutPivoting >=================
//	============================================================================

class BzzFactorizedTridiagonalWithoutPivoting : public BzzFactorizedTridiagonalBase
{
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector& b,
		BzzVector* x);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, const BzzMatrix& B,
		BzzMatrix* X);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzMatrix* BX);

public:
	// Default
	BzzFactorizedTridiagonalWithoutPivoting(void) {};
	// Copy-Initializer
	BzzFactorizedTridiagonalWithoutPivoting(const BzzFactorizedTridiagonalWithoutPivoting& R)
	{
		BzzError("%s %d", BZZ_ERR_IMPLEMENTATION, R.whoAmI);
	}
	// Generic
	BzzFactorizedTridiagonalWithoutPivoting(BzzVector* ll, BzzVector* dd,
		BzzVector* uu);
	BzzFactorizedTridiagonalWithoutPivoting(char* fileBzzFactorizedTridiagonal);

	virtual void operator ()
		(BzzVector* ll, BzzVector* dd,
			BzzVector* uu);
	virtual void Factorization(void);
	virtual void Solution(BzzVector& b, BzzVector* x);
	virtual void Solution(BzzVector* bx);
	virtual void Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void Solution(BzzMatrix* BX);
};

//	============================================================================
//	==================< class BzzFactorizedTridiagonalGauss >===============
//	============================================================================

class BzzFactorizedTridiagonalGauss : public BzzFactorizedTridiagonalBase
{
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector& b,
		BzzVector* x);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, const BzzMatrix& B,
		BzzMatrix* X);
	friend void Solve
	(BzzFactorizedTridiagonalBase& A, BzzMatrix* BX);

private:
	BzzVector r;
	int* indx;

public:
	// Default
	BzzFactorizedTridiagonalGauss(void) {};
	// Copy-Initializer
	BzzFactorizedTridiagonalGauss(const BzzFactorizedTridiagonalGauss& R)
	{
		BzzError("%s %d", BZZ_ERR_IMPLEMENTATION, R.whoAmI);
	}
	// Generic
	BzzFactorizedTridiagonalGauss(BzzVector* ll,
		BzzVector* dd, BzzVector* uu);
	BzzFactorizedTridiagonalGauss(char* fileBzzFactorizedTridiagonal);
	virtual ~BzzFactorizedTridiagonalGauss(void)
	{
		delete indx;
	}

	virtual void operator ()
		(BzzVector* ll, BzzVector* dd, BzzVector* uu);
	virtual void Factorization(void);
	virtual void Solution(BzzVector& b, BzzVector* x);
	virtual void Solution(BzzVector* bx);
	virtual void Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void Solution(BzzMatrix* BX);
};

#endif // TRIDIAGONAL_DOUBLE_HPP