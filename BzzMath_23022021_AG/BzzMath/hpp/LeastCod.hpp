// BZZMATH: Release 3.1

//	=============================< LEASTCOD.HPP >===============================
//	* Class BzzLestCorrection															*
// * Examples: c:\bzzmath\exampled\dxLeastc.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-1996	Date Written

//	============================================================================
//	* Use BzzLeastCorrection class to find the minimum norm2 correction	*
// * dx subject to:																				*
//	* x = x0 + dx																					*
//	* xMin <= x <= xMax																			*
//	* Ax = a																							*
//	* Bx >= b																						*
//	****************************************************************************
//	* Constructors:																				*
//	* BzzLeastCorrection c; // default												*
//	* BzzLeastCorrection c(&xMin,&xMax,&A,&a,&B,&b);							*
//	****************************************************************************
//	* Utilization:																					*
//	* BzzLeastCorrection c(&xMin,&xMax,&A,&a,&B,&b);							*
//	* dx = c(x0); x = x0 + dx;																	*
//	* // or																							*
//	* c(&x);	// Input: x; Output: x + dx													*
//	* ************																					*
//	* BzzLeastCorrection c;																*
//	* c(&xMin,&xMax,&A,&a,&B,&b);																*
//	* dx = c(x0); x = x0 + dx;																	*
//	* // or																							*
//	* c(&x);																							*
//	* ************																					*
//	* BzzLeastCorrection c(&xMin1,&xMax1,&A1,&a1,&B1,&b1);					*
//	* c(&x);																							*
//	* Delete(&c);																					*
//	* c(&xMin2,&xMax2,&A2,&a2,&B2,&b2);														*
//	* c(&x);																							*
//	****************************************************************************
//	* Other function:																				*
//	* void GetBzzMatrixN(BzzMatrix *N);												*
//	*	 calculates the matrix N null space of A											*
//	****************************************************************************

#ifndef BZZ_LEAST_CORRECTION_DOUBLE_HPP
#define BZZ_LEAST_CORRECTION_DOUBLE_HPP

//	============================================================================
//	=====================< class BzzLeastCorrection >========================
//	============================================================================
class BzzLeastCorrection : public BzzBaseClass
{
	friend void Delete(BzzLeastCorrection* c);
	friend class BzzConstrainedMinimum;
	friend double BzzPenalityAuxiliryFun(BzzVector& y);
	friend double BzzPenalityAuxiliryFunMono(double y);

private:
	char firstStep;
	int	numVariables,
		numXMin,
		numXMax,
		numLinearEquality,
		numLinearInequality;

	BzzVector	xMin,
		xMax,
		a,
		b,
		y,
		t,
		dx;

	BzzMatrix	A,
		B,
		T,
		N;

	BzzFactorizedLQ LQ;

	//	============================================================================
	//	*********************< Functions for linear algebra >***********************
	//	============================================================================
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
public:
	// default constructor
	BzzLeastCorrection(void);

	// copy constructor
	BzzLeastCorrection(const BzzFactorized& rval)
	{
		BzzError("No default constructor");
	}

	// constructor c(&xMin,&xMax,&A,&a,&B,&b);
	BzzLeastCorrection(BzzVector* xmi, BzzVector* xma,
		BzzMatrix* AA, BzzVector* aa,
		BzzMatrix* BB, BzzVector* bb);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzLeastCorrection(void) {};

	//	============================================================================
	//	******************************< Functions >*********************************
	//	============================================================================
	virtual void ObjectBzzPrint(void) {};
	BzzVector operator () (BzzVector& x);
	void operator () (BzzVector* x);
	void operator () (BzzVector* xmi, BzzVector* xma,
		BzzMatrix* AA, BzzVector* aa,
		BzzMatrix* BB, BzzVector* bb);

	void GetBzzMatrixN(BzzMatrix* N); // Null Space of Matrix A
};

#endif // BZZ_LEAST_CORRECTION_DOUBLE_HPP