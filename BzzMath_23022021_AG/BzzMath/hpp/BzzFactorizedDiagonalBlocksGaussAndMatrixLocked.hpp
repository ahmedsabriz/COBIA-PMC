// BZZMATH: Release 7.0

/*

*/

//	================< BzzFactorizedDiagonalBlocksGaussAndMatrixLocked.hpp >==============
//	* class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked: 									*
//	* This class is used for solving linear systems with block diagonal 			*
//	* structure and other sparse element													*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	11-2004	Date Written

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedDiagonalBlocksGaussAndMatrixLocked A("FileName",B);						*
//	****************************************************************************
//	****** Function for solving systems:													*
//	* Solve(&A,&X);																				*

#ifndef FACTORIZED_DOUBLE_DIAGONAL_GAUSS_AND_MATRIX_LOCKED_HPP
#define FACTORIZED_DOUBLE_DIAGONAL_GAUSS_AND_MATRIX_LOCKED_HPP

//	============================================================================
//	===============< class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked >=============
//	============================================================================
class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int	whoAmI,
		numDiagonalMatrices,
		numVariables,
		blocksDimensions,
		numIterations,
		maxIterations;

	BzzVectorInt diag;

	double	tollAbs,
		tollRel,
		maxDx,
		maxX;

	enum BzzFactorizedDiagonalBlocksGaussAndMatrixLockedType
	{
		UNKNOWN,
		MATRIX_FILE,
		MATRIX_MEMO,
		VECTOR_FILE,
		VECTOR_MEMO
	}bzzFactorizedDiagonalBlocksGaussAndMatrixLockedType;

	char fileF[20];

	//	BzzMatrix X0;
	BzzVector x0;
	//	BzzMatrix A;
	BzzFactorizedGauss G;
	BzzMatrixSparseLockedByRows Q;
	BzzFactorizedDiagonalBlocksGauss DA;

	//	void Initialize(int numvar,int nummat);	// initialising constructors

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked(void);
	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked(const char* fileM,
		BzzMatrixSparseLockedByRows* QQ);
	void operator()(const char* fileM, BzzMatrixSparseLockedByRows* QQ);
	void operator()(const char* fileM);

	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked(BzzMatrixDiagonalBlocks* AA,
		BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzMatrixDiagonalBlocks* AA,
		BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzMatrixDiagonalBlocks* AA);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedDiagonalBlocksGaussAndMatrixLocked(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	//	int GaussIterativeSolve(BzzMatrix *b);
	int GaussIterativeSolve(BzzVector* b);

	virtual void ObjectBzzPrint(void);

	void SetToleranceAbs(double tolA) { tollAbs = tolA; }
	void SetToleranceRel(double tolR) { tollRel = tolR; }
	void SetMaxIterations(int max) { maxIterations = max; }
	int GetIterations(void) { return numIterations; }
	//	void SetStartX(BzzMatrix *XX);
	void SetStartX(BzzVector* xx);
	void SetStartX(double v);
	double GetMaxX(void) { return maxX / double(numVariables); }
	double GetMaxDx(void) { return maxDx / double(numVariables); }

	//	friend void Delete(BzzFactorizedDiagonalBlocksGaussAndMatrixLocked *A);
	//	friend void ChangeDimensions(int numv,int mat,
	//			 BzzFactorizedDiagonalBlocksGaussAndMatrixLocked *result);
};

#endif // FACTORIZED_DOUBLE_DIAGONAL_GAUSS_AND_MATRIX_LOCKED_HPP