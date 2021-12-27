// BZZMATH: Release 7.0

//	===================< BzzNonLinearRegression.hpp >=====================
//	* BzzNonLinearRegression class for linear regression						*
//	* Description: 																				*
//	*					by G. Buzzi-Ferraris														*
// *																									*
// * Examples: c:\bzzmath\examples\BzzMathAdvanced\Regressions\					*
// *           NonLinearRegression\NonLinearRegression.cpp			*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	02-2000	Date Written.
//	10-2005	Non deterministic Independent variables.
//	02-2010	New numTij.
// 08-2011	Added new class BzzFindTrends

//	============================================================================
// S(b) = S {[w*(ys - yt)]**2}																*
//	******	BzzNonLinearRegression constructors:										*
//	* BzzNonLinearRegression r(numModels,X,Y,W,ModelEx1);								*
//	* BzzNonLinearRegression r(numModels,X,Y,W,ModelEx1,ConstraintsEx1);			*
//	* BzzVector b(3,11.,3000.,-40.);															*
//	* r.InitializeModel(1,b);																	*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = r.NumExperiments(); 																	*
//	* j = r.NumParameters(); 																	*
//	* who = r.WhoAmI();																			*
//	* int count = BzzNonLinearRegression::ObjectCount();								*
//	* int countInScope = BzzNonLinearRegression::ObjectCountInScope();			*
//	****************************************************************************
//	* Functions:																					*
//	* BzzVector b,bl,bu,x0,y0,w0,v;															*
//	* BzzVectorInt iv;																			*
//	* int i,m;																						*
//	* r.BzzPrint("Comment"); // always														*
//	* r.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	* r.InsertExperiment(x0,y0,w0);															*
//	* r.DeleteExperiment(i);																	*
//	* r.DeleteExperiments(iv);																	*
//	* b = r.GetParameters(m);																	*
//	* bl = r.Get95LowerConfidenceBounds(m);												*
//	* bu = r.Get95UpperConfidenceBounds(m);												*
//	* bl = r.Get99LowerConfidenceBounds(m);												*
//	* bu = r.Get99UpperConfidenceBounds(m);												*
//	* v = r.GetResiduals(m);																	*

//	* iv = r.GetOutliers(m);																	*
//	* b = r.GetRobustParameters(m);															*
//	* v = r.GetRobustResiduals(m);															*
//	* Delete(&r);																					*
//	****************************************************************************

#ifndef NON_LIN_REG_DOUBLE_HPP
#define NON_LIN_REG_DOUBLE_HPP

//	============================================================================
//	=====================< class BzzNonLinearRegression >=======================
//	============================================================================
class BzzMatrixSparseLockedByRows;
class BzzNonLinearRegression : public BzzBaseClass
{
	friend double BzzMinimumNonLinearRegression(BzzVector& b);
	friend double BzzMinimumNonLinearRegressionRobust(BzzVector& b);
	friend double BzzMinimumXNonLinearRegression(BzzVector& xx);

	enum BzzNonLinearRegressionIdependentVariablesType
	{
		TYPE_DETERMINISTIC = 1,
		TYPE_NON_DETERMINISTIC = 2
	}nonLinearRegressionIdependentVariablesType;

	enum BzzNonLinearRegressionVarianceType
	{
		TYPE_WEIGHT = 1,
		TYPE_VARIANCE = 2,
		TYPE_W_CONST = 3
	}nonLinearRegressionVarianceType;

	enum BzzNonLinearRegressionFindExperimentsType
	{
		TYPE_XMIN_XMAX = 1,
		TYPE_XX = 2,
		TYPE_XX_WXX = 3
	}nonLinearRegressionFindExperimentsType;

private:
	static const char* const BZZ_ERROR;
	static const char* const WARNING;

	static int count; // for whoAmI
	static int countInScope;

	static double	t95_1_20[],
		t95i,
		t99_1_20[],
		t99i;

	char* bzzFileSave;
	char saveResults;

	int	numModels,
		numModelsPresent,
		numExperiments,
		numX,
		numY,
		numEY,
		whoAmI,
		numGrid,
		numExperimentsXX,
		numRowsFx,
		maxIterations,
		printTasks,
		printSubTasks;
	int missingExperiments;

	char war;

	double	thresholdForMaxMinX,
		tijTolerance;
	BzzVectorInt maxNumTij;

	int selectedCriterion, usedCriterion;
	//	nonLinReg.SelectCriterion(0); // default
	//	nonLinReg.SelectCriterion(1); // determinant
	//	nonLinReg.SelectCriterion(2); // dMIn
	//	nonLinReg.SelectCriterion(3); // fTFTF-1f
	//	nonLinReg.SelectCriterion(4); // fTFTF-1f
	//	nonLinReg.SelectCriterion(5); // Sum(di)
	//	nonLinReg.SelectCriterion(6); // Inverso condizionamaneto
	//	nonLinReg.SelectCriterion(7); // sphericity
	//	nonLinReg.SelectCriterion(8); // maxXDistance
	//	nonLinReg.SelectCriterion(9); // maxFDistance
	//	nonLinReg.SelectCriterion(10); // maxXDistance in principal axes
	//	nonLinReg.SelectCriterion(11); // maxFDistance in principal axes

	BzzMatrix	X, Xx,
		Z,
		W, Wx,
		Y,
		J,
		VbInv,
		XX,
		WXX,
		UWXX2,
		Fs, // J non pesato
		Fx, // F in corrispondenza di X new (per i vari modelli)
		Fsx, // J + numY esperimento i-esimo preso da Fx
		Yt, // TODO
		S2; // TODO!!!!!!!!!!!!

	BzzVector	x,
		yy,
		w,
		uw2, // s2 = 1./(w*w)
		z,
		zw,
		* y,
		* yw,
		* r,
		* rw,
		* g, // for constraints
		s2,
		b,
		bs2,
		db,
		* mb0,
		* mb,
		bRobust,
		* mbRobust,
		* mbMin,
		* mbMax,
		mFInitial,
		mFFinal,
		mcondition,
		mbNorm2,
		bStandardDeviation,
		* mbStandardDeviation,
		* myStandardDeviation,
		mrNorm2,
		msse,
		mmse,
		mr2,
		mstandardDeviation,
		mt95,
		mt99,
		* yXX;

	BzzVector sigma2, sigma2x;
	int varianceDegreeOfFreedom;
	int numNonDetX;

	double	sst,
		sigma;
	//	BzzMatrixSymmetric *mVb,Vb;
	BzzMatrixSymmetric* mVb;

	BzzVectorInt	numParameters,
		dof,
		leastParametersEvaluated,
		robustParametersEvaluated,
		highlyProbableOutliers,
		possibleOutliers,
		highlyProbableOutliersY,
		possibleOutliersY,
		presenceModel;
	//						initializeModel,
	//						includeModel;

	BzzVectorInt	iDistanceFxfromFs,
		iDistanceXfromXs,
		iDistanceFxfromFsPrincipal,
		iDistanceXfromXsPrincipal,
		ivolumeFunction;

	BzzVector	volumeFunction,
		volumeFunctionIm,
		mseRatio,
		d,
		normColumn,
		Tijm,
		distanceFxfromFs,
		distanceXfromXs,
		distanceFxfromFsPrincipal,
		distanceXfromXsPrincipal,
		bestNewExperiment,
		bestNewExperiment8,
		bestNewExperiment9,
		bestNewExperiment10,
		bestNewExperiment11;

	BzzMatrix globalBestNewExperiment;

	BzzFactorizedSVD svdF, svdX, svdFB;
	BzzFactorizedQR qr;

	BzzMatrix* mV, // SVD
		* mU,	// SVD
		* mP; // Principal components

	BzzVector* md,	  // singular values
		mdMax,
		xMin,
		xMax;

	//	BzzVector	dMin,
	//					condition,
	//					minMaxS,
	//					determinant,
	//					dSum,
	//					dMax;

	//	BzzVectorInt	idMin,
	//						icondition,
	//						iminMaxS,
	//						ideterminant,
	//						idSum,
	//						idMax;

	void (*ptrModel)(int m, int e, BzzVector& bb, BzzVector& xx, BzzVector& yy);
	void (*ptrConstraints)(int m, BzzVector& bb, BzzVector& xx, BzzVector& gg);

	//	BzzFactorizedQR G;
	//	BzzFactorizedSVD S;

	//	BzzVector	y,	  // observed dependent variable
	//					yex,	// estimate dependent variable
	//					b,	  // parametres
	//					r,	  // residuals
	//					d,	  // singular values
	//					fMed,
	//					yexRobust,
	//					bRobust,
	//					rRobust,
	//					rRobustStandardized,
	//					bStandardDeviation,
	//					fMedian,
	//					fMad,
	//					yStandardized,
	//					mseBest;

	//	double	sigma2,
	//			mse,
	//			sse,
	//			sst,
	//			ssm,
	//			r2,
	//			t95,
	//			t99,
	//			F95,
	//			F99,
	//			dMax,
	//			condition,
	//			bNorm2,
	//			rNorm2,
	//			standardDeviation,
	//			s2,
	//			yMed,
	//			yMedian,
	//			yMad,
	//			rRobustMedian,
	//			rRobustMad,
	//			robustScale,
	//			fMinRobust;

	//	char	normalState,
	//			robustState,stepWiseState;

	//	void LeastSquaresRegression(void);
	void Initialize(void);
	//	void Reinitialize(void);
	void LeastSquaresRegression(int m);
	void RobustRegression(int m);
	void Jacobian(int m);
	void MatrixFxFxs(int m);
	void ParametersEvaluation(int m);
	void NonDeterministicXEvaluation(int m);
	void RobustParametersEvaluation(int m);
	void SortXX(void);
	void PrintExperiments(int iM);
	void PrintExperiments(void);
	void GetBTolerance(int m);
	BzzVector bTolerance;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzNonLinearRegression(void);
	BzzNonLinearRegression(const BzzNonLinearRegression& rval);
	// constructor
	BzzNonLinearRegression(int nm, BzzMatrix& XX, BzzMatrix& YY, BzzMatrix& WW,
		void (*ptrM)(int m, int e, BzzVector& bb, BzzVector& XX, BzzVector& YY));
	BzzNonLinearRegression(int nm, BzzMatrix& XX, BzzMatrix& YY,
		void (*ptrM)(int m, int e, BzzVector& bb, BzzVector& XX, BzzVector& YY));

	BzzNonLinearRegression(int m, BzzMatrix& XX,
		BzzMatrixSparseLockedByRows& YY, int dof, BzzVector& s2,
		void (*ptrM)(int m, int e, BzzVector& bb, BzzVector& xx, BzzVector& yy));
	void operator()(int m, BzzMatrix& XX,
		BzzMatrixSparseLockedByRows& YY, int dof, BzzVector& s2,
		void (*ptrM)(int m, int e, BzzVector& bb, BzzVector& xx, BzzVector& yy));

	void operator()(int nm, BzzMatrix& XX, BzzMatrix& YY,
		void (*ptrM)(int m, int e, BzzVector& bb, BzzVector& XX, BzzVector& YY));

	void operator()(BzzMatrix& XX);

	// TODO ?????
//	void operator()(int nm,BzzMatrix &XX,BzzMatrix &YY,
//		void (*ptrM)(int m,int e,BzzVector &bb,BzzVector &xx,BzzVector &yy),
//		void (*ptrC)(int m,BzzVector &bb,BzzVector &xx,BzzVector &gg));

	void InitializeModel(int m, BzzVector& bb);
	void InitializeModel(int m, BzzVector& bb, BzzVector& bbMin, BzzVector& bbMax);
	//	void InitializeModel(int m,BzzVector &bb,BzzVector &bbMin,BzzVector &bbMax,
	//		BzzVector &gg);

	void SelectCriterion(int c)
	{
		if (c < 0 || c > 8)
			selectedCriterion = 0;
		else
			selectedCriterion = c;
	}
	//	void IncludeModel(int m){includeModel[m] = 1;}
	//	void ExcludeModel(int m){includeModel[m] = 0;}

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzNonLinearRegression(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	void Save(char* sv)
	{
		bzzFileSave = sv;
		saveResults = 1;
	}
	// number of experiments
	int NumExperiments(void) const { return numExperiments; }

	// number of parameters
	int NumParameters(int m) { return numParameters(m); }

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	void SetVariance(int df, BzzVector& s2);
	void SetVariance(BzzVectorInt df, BzzVector& s2);
	void SetTijTolerance(double tijTol)
	{
		if (tijTol > 1.2)tijTolerance = tijTol;
		else tijTolerance = 1.2;
	}

	// ============================================================================
	// ******************************< Setting functions >*************************
	// ============================================================================
	void StopTaskPrint(void) { printTasks = 0; }
	void StopSubTasksPrint(void) { printSubTasks = 0; }
	void SetTasksPrint(void) { printTasks = 1; }
	void SetSubTasksPrint(int psb = 1)
	{
		if (psb > 0)
			printSubTasks = psb;
		else
			printSubTasks = 1;
	}

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================
	BzzVector GetParameters(int m);
	BzzVector Get95LowerConfidenceBounds(int m);
	BzzVector Get95UpperConfidenceBounds(int m);
	BzzVector Get99LowerConfidenceBounds(int m);
	BzzVector Get99UpperConfidenceBounds(int m);
	BzzVector GetResiduals(int m);
	BzzVector GetResidualsStandardized(int m);
	double GetMeanSquareError(int m);

	BzzVectorInt GetCandidateOutliers(int m);
	BzzVectorInt GetHighlyProbableOutliers(int m);
	BzzVectorInt GetPossibleOutliers(int m);

	//	BzzVector GetRobustParameters(void);
	BzzVector GetRobustParameters(int m);
	BzzVector GetRobustResiduals(int m);
	BzzVector GetRobustResidualsStandardized(int m);

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================
	virtual void ObjectBzzPrint(void) {};////////???????????????????!!!!!!!!!
	void DataPrint(void);
	void LeastSquaresAnalysis(void);
	void LeastSquaresAnalysisAndExperimentsSearch(BzzVector& xMi, BzzVector& xMa);
	void LeastSquaresAnalysisAndExperimentsSearch(BzzMatrix& XX);
	void LeastSquaresAnalysisAndExperimentsSearch(BzzMatrix& XX, BzzMatrix& WXX);
	void RobustAnalysis(void);
	void ExtensiveAnalysis(void);
	void LeastSquaresAnalysis(int m);
	void LeastSquaresAnalysisAndExperimentsSearch(int m);
	void RobustAnalysis(int m);
	void ExtensiveAnalysis(int m);
	void CleverLeastSquaresAnalysis(int m);
	void RemoveModel(int iM)
	{
		presenceModel[iM] = 0;
	}
	void SetBalance1011AndDiscrimination(double threshold)
	{
		if (threshold < 0. || threshold > 1.)
			thresholdForMaxMinX = .5;
		else
			thresholdForMaxMinX = threshold;
	}
	void SetMaxIterations(int mx = 0)
	{
		maxIterations = mx;
	}

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void ModifyExperiment(int i, BzzVector& x0, BzzVector& y0, BzzVector& w0);
	void AddExperiment(BzzVector& x0, BzzVector& y0, BzzVector& w0);
	void InsertExperiment(int i, BzzVector& x0, BzzVector& y0, BzzVector& w0);
	void DeleteExperiment(int i);
	void DeleteExperiments(BzzVectorInt& iv);
};

//	============================================================================
//	========================< class BzzFindTrends >=============================
//	============================================================================

class BzzFindTrends : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int	noisePresent, // default 1 (noise present)
		noOutliersPresent; // default 1 (NO outliers)

	BzzVector	pmc,
		abscissas;
	BzzVectorInt	localTrend;
	BzzVector	lowerBoundStrongLocalAndGlobalPositiveTrend,
		lowerBoundStrongLocalAndGlobalNegativeTrend,
		lowerBoundDoubiousLocalAndGlobalTrend,
		lowerBoundDoubiousLocalAndMaximumGlobal,
		lowerBoundDoubiousLocalAndMinimumGlobal,
		lowerBoundDoubiousLocalAndPositiveGlobalTrend,
		lowerBoundDoubiousLocalAndNegativeGlobalTrend,

		upperBoundStrongLocalAndGlobalPositiveTrend,
		upperBoundStrongLocalAndGlobalNegativeTrend,
		upperBoundDoubiousLocalAndGlobalTrend,
		upperBoundDoubiousLocalAndMaximumGlobal,
		upperBoundDoubiousLocalAndMinimumGlobal,
		upperBoundDoubiousLocalAndPositiveGlobalTrend,
		upperBoundDoubiousLocalAndNegativeGlobalTrend;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzFindTrends(void);
	void operator ()(char* fileBzz,
		BzzVector& x, BzzVector& y,
		BzzVector* l1,
		BzzVector* l2,
		BzzVector* l3,
		BzzVector* l4,
		BzzVector* l5,
		BzzVector* l6,
		BzzVector* l7,
		BzzVector* u1,
		BzzVector* u2,
		BzzVector* u3,
		BzzVector* u4,
		BzzVector* u5,
		BzzVector* u6,
		BzzVector* u7,
		BzzVectorInt* io);
	void GetParameters(BzzVector* ppmc)
	{
		Swap(ppmc, &pmc);
	}
	void GetAbscissas(BzzVector* absc);
	double CompareTrends(char* fileBzz, BzzVector& absc, BzzVector& pmc1, BzzVector& pmc2, int perc = 0);
	void CompareTrends(BzzVector& absc, BzzVector& trend, BzzMatrix& Trends,
		BzzVector* e, int perc = 0);
	virtual void ObjectBzzPrint(void);
	void NoisePresent(void)
	{
		noisePresent = 1;
	} // default 1 (noise present)
	void NoNoisePresent(void)
	{
		noisePresent = 0;
	} // default 1 (noise present)
	void OutliersPresent(void)
	{
		noOutliersPresent = 0;
	} // default 1 (NO outliers)
	void NoOutliersPresent(void)
	{
		noOutliersPresent = 0;
	} // default 1 (NO outliers)
};

#endif // NON_LIN_REG_DOUBLE_HPP