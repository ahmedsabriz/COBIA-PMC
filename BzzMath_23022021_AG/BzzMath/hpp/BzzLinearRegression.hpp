// BZZMATH: Release 7.0

//	==================< BzzLinearRegression.hpp >=========================
//	* BzzLinearRegression class for linear regression							*
// * Examples: c:\bzzmath\examples\BzzMathAdvabced\Regressions\					*
// *				LinearRegression\LinearRegression.cpp					*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	02-1995	Date Written.

////////////////// Release 4.0
//	03-2001	Added default constructor.
//	03-2001	Added operator(F,y).
//	04-2001	Added P-Value in statistical analysis.

////////////////// Release 5.0
//	10-2003	Added GetR2.
//	10-2003	Added GetAdjustedR2.
//	10-2003	Added GetBTolerance
//	10-2003	Added GetBzzColumnsCorrelation.

////////////////// Release 6.0
//	10-2003	Added selection costructor.
//	03-2009	Added GetStudentizedResiduals dunction.
//	03-2009	Added CleverLeastSquaresAnalysis function.
//	04-2009	Added HomoscedasticityAnalysis function.

//	============================================================================
//	******	BzzLinearRegression constructors:									*
//	* BzzLinearRegression r(F,y); 													*
// * BzzLinearRegression r("REGR.DAT");											*
//	* BzzLinearRegression r;r(F,y);													*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = r.NumExperiments(); 																	*
//	* j = r.NumParameters(); 																	*
//	* who = r.WhoAmI();																			*
//	* int count = BzzLinearRegression::ObjectCount();									*
//	* int countInScope = BzzLinearRegression::ObjectCountInScope();				*
//	****************************************************************************
//	* Functions:																					*
//	* BzzVector b,bl,bu,f0,v;																	*
//	* BzzVectorInt iv;																			*
//	* double y0;																						*
//	* int i;																							*
//	* r.BzzPrint("Comment"); // always														*
//	* r.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	* r.InsertExperiment(f0,y0);																*
//	* r.DeleteExperiment(i);																	*
//	* r.DeleteExperiments(iv);																	*
//	* b = r.GetParameters();																	*
//	* bl = r.Get95LowerConfidenceBounds();													*
//	* bu = r.Get95UpperConfidenceBounds();													*
//	* bl = r.Get99LowerConfidenceBounds();													*
//	* bu = r.Get99UpperConfidenceBounds();													*
//	* v = r.GetResiduals();																		*

//	* iv = r.GetOutliers();																		*
//	* b = r.GetRobustParameters();															*
//	* v = r.GetRobustResiduals();																*
//	* Delete(&r);																					*
//	****************************************************************************

#ifndef LINREG_DOUBLE_HPP
#define LINREG_DOUBLE_HPP

//	============================================================================
//	======================< class BzzLinearRegression >=========================
//	============================================================================

class BzzLinearRegression : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const char* const WARNING;

	static int count; // for whoAmI
	static int countInScope;

	static double	t95_1_20[],
		t95i,
		t99_1_20[],
		t99i;
public:
	BzzFactorizedSVD svdF;
	BzzFactorizedSVD svdS;

	BzzMatrix	F, // Original data
		V, // SVD
		U,	// SVD
		FStandardized,
		P,	// Principal Component
		bBest, // StepWise Analysis
		tBest, // StepWise Analysis
		VbInv; // variance inflation factor

	BzzFactorizedGauss G;

	BzzMatrixSymmetric Vb;

	BzzVector	y,	  // observed dependent variable
		yex,	// estimate dependent variable
		b,	  // parametres
		r,	  // residuals
		d,	  // singular values
		fMed,
		bTolerance,
		bzzColumnsCorrelation,
		yexRobust,
		bRobust,
		rRobust,
		rRobustStandardized,
		bStandardDeviation,
		yStandardDeviation,
		fMedian,
		fMad,
		yStandardized,
		mseBest,
		bNormSVD,
		rNormSVD;

	BzzVectorInt	outliersCandidate,
		highlyProbableOutliers,
		possibleOutliers;

	int	numExperiments,	  // data number
		numParameters,  // parameters number
		dof,
		intercept, //  = 1 model with intercept
		whoAmI,
		iQRHigh, iQRProb,
		varianceDegreeOfFreedom,
		numWarningD;

	double	sigma2,
		mse,
		sse,
		sst,
		ssm,
		r2,
		adjustedR2,
		t95,
		t99,
		F95,
		F99,
		dMax,
		condition,
		systemConditioning,
		bNorm2,
		rNorm2,
		standardDeviation,
		s2,
		yMed,
		yMedian,
		yMad,
		rRobustMedian,
		rRobustMad,
		robustScale,
		fMinRobust,
		dMaxRob,
		dMinRob;

	char	normalState,
		warCorr,
		conditionWarning,
		principalWarning,
		parametersWarning,
		robustState, stepWiseState;

	void LeastSquaresRegression(void);
	void RobustRegression(void);
	void StepWiseRegression(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzLinearRegression(void);
	BzzLinearRegression(const BzzLinearRegression& rval);
	// constructor linReg(F,yy);
	BzzLinearRegression(const BzzMatrix& FF,
		const BzzVector& yy);
	//da FILE ASCII
	BzzLinearRegression(char* filereg);
	BzzLinearRegression(BzzMatrix& FF, BzzVector& yy,
		BzzVectorInt& jF, BzzVectorInt& iE);
	BzzLinearRegression(BzzMatrix& FF, int jy,
		BzzVectorInt& jF, BzzVectorInt& iE);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzLinearRegression(void) { countInScope--; }

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
		// number of experiments
	int NumExperiments(void) const { return numExperiments; }

	// number of parameters
	int NumParameters(void) const { return numParameters; }

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	void SetVariance(double s2, int degfr)
	{
		if (s2 > 0. && degfr > 0)
		{
			sigma2 = s2;
			varianceDegreeOfFreedom = degfr;
		}
	}

	void SetVariance(int degfr, double s2)
	{
		if (s2 > 0. && degfr > 0)
		{
			sigma2 = s2;
			varianceDegreeOfFreedom = degfr;
		}
	}

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================
	BzzVector GetParameters(void);
	BzzVector Get95LowerConfidenceBounds(void);
	BzzVector Get95UpperConfidenceBounds(void);
	BzzVector Get99LowerConfidenceBounds(void);
	BzzVector Get99UpperConfidenceBounds(void);
	BzzVector GetResiduals(void);
	BzzVector GetResidualsStandardized(void);
	double GetR2(void);
	double GetAdjustedR2(void);
	double GetMeanSquareError(void);
	BzzVector GetBTolerance(void);
	BzzVector GetBzzColumnsCorrelation(void);

	BzzVectorInt GetCandidateOutliers(void);
	BzzVectorInt GetHighlyProbableOutliers(void);
	BzzVectorInt GetPossibleOutliers(void);

	BzzVector GetRobustParameters(void);
	BzzVector GetRobustResiduals(void);
	BzzVector GetRobustResidualsStandardized(void);
	void GetStudentizedResiduals(BzzVector* str);
	void StudentizedResiduals(void);

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================
	virtual void ObjectBzzPrint(void);
	void DataPrint(void);
	void LeastSquaresAnalysis(void);
	void RobustAnalysis(void);
	void SVDAnalysis(void);
	void ExtensiveAnalysis(void);
	void StepWiseAnalysis(void);
	void StepwiseAnalysis(void) { StepWiseAnalysis(); }
	void CleverLeastSquaresAnalysis(int nMse = 0);
	void HomoscedasticityAnalysis(void);
	void InfluentialObservations(BzzVectorInt* ifo = 0,
		BzzVector* infO = 0, BzzVector* inf = 0);
	void GetSortedBzzVectorDAndBzzMatricesVU(BzzVector* d, BzzMatrix* V,
		BzzMatrix* U);
	void GetPrincipalComponents(BzzMatrix* P);
	void GetBestStepWiseModels(BzzVector* mseBest, BzzMatrix* bBest);
	void GetBestStepwiseModels(BzzVector* mseBest, BzzMatrix* bBest)
	{
		GetBestStepWiseModels(mseBest, bBest);
	}

	double FunzMinimum(void);
	double MinimumSimplex(BzzVector& x);

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void operator()(BzzMatrix& FF, BzzVector& yy);
	void Initialize(BzzMatrix& FF, BzzVector& yy);
	void ModifyExperiment(int i, BzzVector& f0, double y0);
	void InsertExperiment(int i, BzzVector& f0, double y0);
	void DeleteExperiment(int i);
	void DeleteExperiments(BzzVectorInt& iv);

	void ModifyFunction(int j, BzzVector& fj);
	void InsertFunction(int j, BzzVector& fj);
	void DeleteFunction(int j);
};

#endif // LINREG_DOUBLE_HPP