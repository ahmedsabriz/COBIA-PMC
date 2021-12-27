// BZZMATH: Release 3.1

//	=============================< FACTBLOD.HPP >===============================
//	* class BzzFactorizedLQMatrixBlocks													*
//	* This class is used for solving linear systems with block structure			*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6,7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exfatblo.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	03-1995	Date Written
//	04-1997	Added Reinitialize.
//	05-1997	Added ConditionNumber.

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedLQMatrixBlocks A(numVariables,numDiagonalMatrices);			*
//	****************************************************************************
//	****** Insert SubMatrices																	*
//	* A.Insert(i,j,&aij);																		*
//	****************************************************************************
//	****** Functions for solving systems:													*
//	* Solve(&A,&bx);																				*
//	* A.SolveAndDEstroy(&A,&bx);																*
//	****************************************************************************
//	* Other functions:																			*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef FACTBLOCLQD_HPP
#define FACTBLOCLQD_HPP

class ElementBzzFactorizedLQMatrixBlocks
{
	friend class BzzFactorizedLQMatrixBlocks;
private:
	int column;
	BzzMatrix alk;
	ElementBzzFactorizedLQMatrixBlocks* next;
};

//	============================================================================
//	=================< class BzzFactorizedLQMatrixBlocks >=====================
//	============================================================================

class BzzMatrix;
class BzzVector;
class BzzVectorInt;
class BzzFactorized;
class BzzFactorizedPLR;
class BzzFactorizedLQ;

class BzzFactorizedLQMatrixBlocks : public BzzBaseClass
{
	friend class ElementBzzFactorizedLQMatrixBlocks;
private:
	enum BzzFactorizedLQMatrixBlocksStatus
	{
		AVAILABLE,
		DESTROYED,
		FACTORIZED
	}matrixBlockStatus;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int whoAmI;
	int	numDiagonalMatrices,
		numVariables;

	BzzVector b;
	BzzVectorInt nb;
	BzzFactorizedLQ* All;
	ElementBzzFactorizedLQMatrixBlocks** elementRow;

	void Initialize(int nvar, int numMat);	// initialising constructors
	BzzMatrix* GetSubBzzMatrix(int l, int k);
	void DeleteSubBzzMatrix(int l, int k);
	void DeleteBzzMatrix(void);

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedLQMatrixBlocks(void);
	BzzFactorizedLQMatrixBlocks(int numVar, int numMat);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedLQMatrixBlocks(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	friend void SolveAndDestroy(BzzFactorizedLQMatrixBlocks* ab,
		BzzVector* bx);
	friend void Solve(BzzFactorizedLQMatrixBlocks* ab, BzzVector* bx);
	double ConditionNumber(void);
	void GetBzzVectorD(BzzVector* dd);

	virtual void ObjectBzzPrint(void);
	int CountSubMatrices(void);
	void Insert(int l, int k, BzzMatrix* alk);
	int GetSubMatrixExistence(int row, int column);
	void Reinitialize(int nvar, int numMat);
	BzzFactorizedLQMatrixBlocks& operator = (BzzMatrixBlocks& m);
	void ReplaceForOde(double c, BzzMatrixBlocks& rval);
};

#endif // FACTBLOCLQD_HPP