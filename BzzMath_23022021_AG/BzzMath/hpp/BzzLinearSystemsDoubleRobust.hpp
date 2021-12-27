// BZZMATH: Release 7.0

// ===============< BzzLinearSystemsRobust.HPP >=========================
// * Class BzzLinearSystemsRobust for sparse no structured linear			*
// * systems solution Sx = b																	*
// * Examples: c:\bzzmath\examples\.cpp													*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2006	Date Written.

// ============================================================================
// ****** Constructors for BzzLinearSystemsRobust:								*
// * BzzLinearSystemsRobust ls; // default										*
// ****************************************************************************
// ****************************************************************************
// ***** Assignment:																				*
// ** ls(&S,&b);																					*
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// * ls.BzzPrint();																				*
// ****************************************************************************
// ***** System solution:																		*
// ** ls(x0);																						*
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_LINEAR_SYSTEMS_DOUBLE_ROBUST
#define BZZ_LINEAR_SYSTEMS_DOUBLE_ROBUST

// ============================================================================
// =================< class BzzLinearSystemsRobust >=====================
// ============================================================================

class BzzLinearSystemsRobust : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int	whoAmI,
		numVariables,
		printTasks,
		printSubTasks,
		maxIterations,
		iter;
	char aitken;
	BzzMatrixSparseLockedByColumns A;
	BzzVector	b,
		norm2C,
		r0, s, w, r1,
		x, r, x0, x1, x2, x3, xI, xD, xB, xF, xMin;
	double f, f0, fI, f3, fD, fB, fF, fMin, dx;
	void FindNewxIxD(void);
	void FindNewxF(void);
	void Initialize(void);
public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzLinearSystemsRobust v;
	BzzLinearSystemsRobust(void);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzLinearSystemsRobust(void) { countInScope--; };

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// ============================================================================
	// ******************************< Setting functions >*************************
	// ============================================================================
	void operator()(BzzMatrixSparseLockedByColumns* AA, BzzVector* bb);
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
	void SetMaxIterations(int maxIt) { maxIterations = maxIt; }
	// ============================================================================
	// *********************************< System Solution >************************
	// ============================================================================
	void operator()(BzzVector* x00);

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void) {};
};

#endif // BZZ_LINEAR_SYSTEMS_DOUBLE_ROBUST