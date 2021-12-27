// BZZMATH: Release 7.0

//	==========================< BzzBVP.hpp >==========================
//	* BzzBVP class: class for BVP				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2014	Date Written

#ifndef BZZ_BVP_HPP
#define BZZ_BVP_HPP

//	===============================================================
//	======================< class BzzBVP >=========================
//	===============================================================
class BzzBVP : public BzzBaseClass
{
private:
	int	numVariables,
		numEquations;

public:
	BzzBVP(void);
};

#endif // BZZ_BVP_HPP