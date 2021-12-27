// BZZMATH: Release 7.0

/*
Linear systems with band structure
and sparse matrix
Example:
*xx......+......
x*xx............
xx*xx......+....
xxx*xx...+.....
.xxx*xx........
..xxx*xx....+..
+..xxx*xx......
....xxx*xx....+
..+..xxx*xx....
......xxx*xx...
.+.....xxx*xx..
........xxx*xx.
.....+...xxx*xx
..........xxx*x
...........xxx*
*/

//	==============< BzzFactorizedBandGaussAndMatrixLocked.hpp >============
//	* class BzzFactorizedBandGaussAndMatrixLocked: 								*
//	* This class is used for solving linear systems band structure					*
//	* and other sparse element																	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	11-2004	Date Written

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedBandGaussAndMatrixLocked A("FileName",B);						*
//	****************************************************************************
//	****** Function for solving systems:													*
//	* Solve(&A,&X);																				*

#ifndef FACTORIZED_DOUBLE_BAND_GAUSS_AND_MATRIX_LOCKED_HPP
#define FACTORIZED_DOUBLE_BAND_GAUSS_AND_MATRIX_LOCKED_HPP

//	============================================================================
//	==============< class BzzFactorizedBandGaussAndMatrixLocked >===========
//	============================================================================
class BzzFactorizedBandGaussAndMatrixLocked : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int	whoAmI,
		lowerBand,
		upperBand,
		numVariables,
		numIterations,
		maxIterations;

	double	tollAbs,
		tollRel,
		maxDx,
		maxX;

	BzzVector x0;
	BzzMatrixSparseLockedByRows Q;
	BzzFactorizedBandGauss DA;

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedBandGaussAndMatrixLocked(void);
	BzzFactorizedBandGaussAndMatrixLocked(BzzMatrixBand* AA,
		BzzMatrixSparseLockedByRows* QQ);
	void operator()(BzzMatrixBand* AA,
		BzzMatrixSparseLockedByRows* QQ);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedBandGaussAndMatrixLocked(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	void GaussIterativeSolve(BzzVector* b);

	virtual void ObjectBzzPrint(void);

	void SetToleranceAbs(double tolA) { tollAbs = tolA; }
	void SetToleranceRel(double tolR) { tollRel = tolR; }
	void SetMaxIterations(int max) { maxIterations = max; }
	int GetLowerBand(void) { return lowerBand; }
	int GetUpperBand(void) { return upperBand; }
	int GetIterations(void) { return numIterations; }
	double GetMaxX(void) { return maxX / double(numVariables); }
	double GetMaxDx(void) { return maxDx / double(numVariables); }
	void SetStartX(BzzVector* xx);
};

#endif // FACTORIZED_BAND_DIAGONAL_GAUSS_AND_MATRIX_LOCKED_HPP