
// BZZMATH: Release 7.0

	extern BzzMatrix O37_8,P37_8;
	extern BzzVector TStar,deltaStar;

// ==========================< BzzGasPropoerties.HPP >=========================
// * Class BzzGasProperties for *
// * Examples: c:\bzzmath\examples\.cpp*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05 - 2003	Date Written.

// ============================================================================
// ****** Constructors for BzzGasProperties:														*
// * BzzGasProperties v; // default																	*
// * BzzGasProperties v = x; // copy-initializer													*
// * BzzGasProperties v(n); // 			*
// ****************************************************************************
// ***** Access functions :																	*
// **
// ****************************************************************************
// ***** Assignment:																				*
// **
// ****************************************************************************
// ****************************************************************************
// ***** Implemented operations :															*
// **
// ****************************************************************************
// ****************************************************************************
// ***** Other functions:																		*
// **
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_GAS_PROPERTIES_HPP
#define BZZ_GAS_PROPERTIES_HPP

// ============================================================================
// ============================< class BzzGasProperties >=============================
// ============================================================================

class BzzGasProperties : public BzzBaseClass
	{
private:
	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double	AVOGADRO,
								BOLTZMANN,
								BOLTZMANN3,
								GAS_PERFETTI,
								P_ATMOSFERICA,
								X_MIN,
								MU_MIN;
//	static BzzMatrix O37_8,P37_8;

	int numComponents;

	BzzVector	epsylonSuKb,
							sigma,
							mu,
							cP,
							cV,
							cVtrans,
							cVrot,
							cVvib,
							alfa,
							zRot298,
							M, // pesi molecolari
							uM, // 1. / pesi molecolari
							m, // masse molecolari
							epsylon,
							alfaStar,
							muStar,
							omega11k,
							omega22k,
							TkStar,
							deltakStar,
							Dkm,
							f298,
							fT,
							zRot,
							rho;

	BzzMatrix	mjk,
							epsylonjkSuKb,
							epsylonjk,
							sigmajk,
							mujk,
							mujkStar,
							deltajkStar,
							TjkStar,
							omega11jk,
							Djk,
							D12,
							Dkk,
							eta;

	BzzVectorInt	forma;

	double	T,
				P;
	BzzVector	Y, // massive
							X; // molare

	double massaMolareMedia;


	// initialise constructors
	void Initialize(void);

	int	whoAmI,
			printTasks,
			printSubTasks;
	void MassaMolareMedia(void);
	void FrazioniMolariXDateLeFrazioniMassiveYNotaLaMassaMolareMedia(void);
	void FrazioniMolariXDateLeFrazioniMassiveY(void);
	void TemperaturaRidotta(void);
	void TemperaturaRidottaComponentePuro(void);
//	void TemperaturaRidotta(double TT);
//	void TemperaturaRidottaComponentePuro(double TT);
	void Omega11jk(void);
	void Omega11k(void);
	void Omega22k(void);
	void CoefficienteDiDiffusioneBinario(void);
	void FattoriCorrettiviD12(void);
	void CoefficienteDiDiffusioneMedio(void);
	void CoefficienteDiAutoDiffusione(void);
	void Viscosit‡DelComponente(void);
	void F298(void);
	void FT(void);
	void Zrot(void);
	void Cv(void);
	void CvVib(void);
	void Densit‡Componente(void);

public:

// ============================================================================
// ***************************< constructors >*********************************
// ============================================================================
	// default constructor BzzGasProperties v;
//	BzzGasProperties(void);
	BzzGasProperties(void)
		{
		BzzError();
		}

	// copy-initializer
	BzzGasProperties(BzzGasProperties &rval);

	// other constructor
	BzzGasProperties(int nc);

// ============================================================================
// *****************************< destructor >*********************************
// ============================================================================
	~BzzGasProperties(void){}

// ============================================================================
// *******************< Non-modifying access functions >***********************
// ============================================================================
	int WhoAmI(void) const {return whoAmI;}
	static int ObjectCount(void){return count;}
	static int ObjectCountInScope(void){return countInScope;}

// ============================================================================
// ******************************< Setting functions >*************************
// ============================================================================
	void StopTaskPrint(void){printTasks = 0;}
	void StopSubTasksPrint(void){printSubTasks = 0;}
	void SetTasksPrint(void){printTasks = 1;}
	void SetSubTasksPrint(int psb = 1)
		{
		if(psb > 0)
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
	BzzGasProperties &operator =
		 (const BzzGasProperties &rval);

// ============================================================================
// *************************< operators for tests >****************************
// ============================================================================
	friend char operator ==
		 (const BzzGasProperties &lval,const BzzGasProperties &rval);

	friend char operator !=
		 (const BzzGasProperties &lval,const BzzGasProperties &rval);

// ============================================================================
// =============================< OPERATIONS >=================================
// ============================================================================


// ============================================================================
// ====================< Non-modifying functions >=============================
// ============================================================================

// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void){}

//	********************************< Save >************************************

// ============================================================================
// ======================< Modifying Functions >===============================
// ============================================================================
	friend void Delete(BzzGasProperties *result); // eliminates BzzGasProperties


// ============================================================================
// =========================< Other functions >================================
// ============================================================================

};

#endif // BZZ_GAS_PROPERTIES_HPP


