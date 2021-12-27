// BZZMATH: Release 7.0

// =============================< BZZREACTION.HPP >=================================
// * Class BzzReaction for *
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
/////////// Release 4.1
//	11-2002	Date Written.

// ============================================================================
// ****** Constructors for BzzReaction:														*

// ****************************************************************************
// ***** Access functions :																	*
// **
// ****************************************************************************
// ***** Assignment:																				*
// **
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// * v.BzzPrint();																				*
// * v.BzzMessage();																				*
// **
// ****************************************************************************
// ***** Other functions:																		*
// **
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_BZZREACTION_HPP
#define BZZ_BZZREACTION_HPP

// ============================================================================
// ============================< class BzzReaction >=============================
// ============================================================================

class BzzReaction : public BzzBaseClass
{
private:
	static const double R_KIN;
	static const double R_CTOT;
	static const double UR_CTOT;
	static const double D1;
	static const double D2;
	static const double D3;
	static const double DSQ;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	// initialise constructors
//	void Initialize(int nc);

	int	whoAmI,
		printTasks,
		printSubTasks;
	int	sameTemperature;
	int	numComponents, // num Rows nuij (200)
		numReactions,	// num columns nuij (3000)
		numEquil; // num reaction with equilibrium

	char** strReaction;
	char** strComponent;

	BzzVectorInt	numDir1,
		numDir2,
		numDir3,
		numDir4,
		numDir5,
		numInvTot1,
		numInvTot2,
		numInvTot3,
		numInvTot4,
		numInvTot5,
		numInvEq1,
		numInvEq2,
		numInvEq3,
		numInvEq4,
		numInvEq5,
		jDir1,
		jDir2,
		jDir3,
		jDir4,
		jDir5,
		jInvTot1,
		jInvTot2,
		jInvTot3,
		jInvTot4,
		jInvTot5,
		jInvEq1,
		jInvEq2,
		jInvEq3,
		jInvEq4,
		jInvEq5,
		jEquil,	// 0 se non c'è 1 se c'è equilibrio
		reactionWithEquil, // indice j della reazione // dimensionato per jEquil=1
		jCorpoMorto, // 0,1,2,3..
		jNumEfficiency,
		* iEfficiency;
	BzzVector valDir5, valInvTot5, valInvEq5;

	BzzVector* efficiency;
	BzzVector	beta1, beta2,
		E1, E2,
		k01, k02,
		aFall, bFall, cFall, dFall, eFall,
		coeffM;

	BzzVector* aDH; // alta T
	BzzVector* bDH; // bassa T
	BzzVector* aDS; // alta T
	BzzVector* bDS; // bassa T
	// TODO a e b per derivate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	BzzVector T1, T2, T3;

	BzzVector	componentDH,
		componentDS,
		reactionDH,
		reactionDS,
		sumNuij;
	//							molecularWeight;
		// r[j] = rDirT[j] *(rDirC[j] + rInvC[j] / kEq[j]);
	BzzVector	r,
		k1, k2,
		rDirT,
		rDirC,
		//							kEq,
		uKeq,
		rInvC;
	int numPoints;
	// Memo for Derivatives
		// One point
	double	wT, mT,
		wP, mP;
	BzzVector	wc, mc,
		wr, mr,
		wrDirT, mrDirT,
		wk1, mk1,
		wk2, mk2,
		wrDirC, mrDirC,
		//wkEq,mkEq,
		wuKeq, muKeq,
		wrInvC, mrInvC,
		wR, mR,
		wcoeffM, mcoeffM,
		//							wF;
		wF, mF,
		logFcent, mlogFcent;

	// Several points
	BzzVector ac;
	BzzVector aR;
	BzzVector	wpT, mpT,
		wpP, mpP;
	BzzMatrix	wpc, mpc,
		wpr, mpr,
		wprDirT, mprDirT,
		wpk1, mpk1,
		wpk2, mpk2,
		wprDirC, mprDirC,
		//wpkEq,mpkEq,
		wpuKeq, mpuKeq,
		wprInvC, mprInvC,
		wpR, mpR,
		wpcoeffM, mpcoeffM,
		wpF, mpF,
		plogFcent, mplogFcent;
	// Derivatives
	//BzzVector	drT;
	//void GetKeqForDerivatives(double T);
	void GetMolecularWeight(int i);
public:
	BzzVector molecularWeight;

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzReaction v;
	BzzReaction(void);

	// copy-initializer
//	BzzReaction(BzzReaction &rval){BzzError("TODO)";}

	// from FILE
	BzzReaction(char* fileSt, char* fileKin, char* fileEq);
	void operator()(char* fileSt, char* fileKin, char* fileEq);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzReaction(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

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

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzReaction& operator =
		(const BzzReaction& rval);

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);

	//	********************************< Save >************************************
	//	void Save(char *file***); // formatted
	//	void Save(char,char *file***);// binary

	// ============================================================================
	// ======================< Modifying Functions >===============================
	// ============================================================================
	friend void Delete(BzzReaction* result); // eliminates BzzReaction

	// recovery from Save
//	friend void Load
//			(BzzReaction *result,char *file***); // formatted
//	friend void Load
//			(BzzReaction *result,char,char *file***);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
	void operator()(double T, double P, BzzVector& c, BzzVector* R);
	void operator()(BzzVector pT, BzzVector pP,
		BzzMatrix& pc, BzzMatrix* pR);
	void ResetNumPoints(int nP);
	void SetBaseForDerivatives(void);

	void DerivativesT(BzzVector* dRT);
	void DerivativesT(int iP, BzzVector* dRT);

	void DerivativesC(BzzMatrix* dRC);
	void DerivativesC(int iP, BzzMatrix* dRC);

	void GetMolarFractionsFromMassFractions(BzzVector& molarFractions,
		BzzVector& massFractions, double wM);
	void GetMassFractionsFromMolarFractions(BzzVector& massFractions,
		BzzVector& molarFractions, double wM);
	void GetConcentrationsFromMolarFractions(BzzVector& concentrations,
		BzzVector& molarFractions, double cTot);
	double GetMeanWeightFromMolarFractions(BzzVector& molarFractions);
	double GetMeanWeightFromMassFractions(BzzVector& massFractions);
};

#endif // BZZ_BZZREACTION_HPP