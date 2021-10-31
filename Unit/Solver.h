#pragma once
#include <COBIA.h>
#include <boost/numeric/odeint.hpp>
#include "MaterialPort.h"
#include "Reaction.h"

using namespace boost::numeric::odeint;
using namespace COBIA;

class Solver :
	public CapeOpenObject<Solver>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Solver> {

	// Members
	CAPEOPEN_1_2::CapeThermoMaterial feed, product;
	std::vector<ReactionPtr> reactions;
	CapeArrayRealImpl molarFlow, totalMolarFlow, T, P, X, Cp;

	// The type of container used to hold the state vector
	typedef /*std::vector<double>*/ CapeArrayRealImpl state_type;

	// TODO: REMOVE HARDCODING
	const CapeReal L = 6;			// [m]
	const CapeReal A = 0.77 / L;	// [m2]

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Solver");
	}

	Solver(CAPEOPEN_1_2::CapeThermoMaterial _feed, CAPEOPEN_1_2::CapeThermoMaterial _product,
		std::vector<ReactionPtr> _reactions) :
		feed(_feed), product(_product), reactions(_reactions) {
	}

	~Solver() {
	}

	void getInitialConditions() {

		// Get inlet molarFlow
		feed.GetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		feed.GetOverallProp(ConstCapeString(COBIATEXT("totalflow")),
			ConstCapeString(COBIATEXT("mole")), totalMolarFlow);

		// Get inlet T[K], P[Pa], X[mol/mol]
		T.resize(1);
		P.resize(1);
		feed.GetOverallTPFraction(T[0], P[0], X);

		// Get heat capacity at constant pressure [J/mol/K]
		CapeArrayStringImpl phaseLabels;
		CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;
		feed.GetPresentPhases(phaseLabels, phaseStatus);

		Cp.resize(phaseLabels.size());
		for (ConstCapeString phaseLabel : phaseLabels) {
			CapeArrayRealImpl phaseCp;
			feed.GetSinglePhaseProp(ConstCapeString(COBIATEXT("heatCapacityCp")),
				phaseLabel, ConstCapeString(COBIATEXT("mole")), phaseCp);
			Cp.emplace_back(phaseCp[0]);
		}
	}

	void pfr(state_type& y, state_type& dydz, CapeReal z) {
		/* Assumuptions:
		*				1D
		*				Homogenous
		*				Empty bed
		*				Adiabatic
		*				Constant heat capacity
		*				No axial dispersion
		*				No radial gradient
		*				No pressure drop
		*/
		// Focus on single component (Ethyl benzene at index 0] for Proof Of Concept
		size_t base = 0;

		// The following variables are overwritten in every step and elaborated only for readability
		// Change to state types to use them in profiling
		CapeReal k, r, R = 0, Qr;
		for (ReactionPtr reaction : reactions) {
			k = reaction->forwarRateConstant * std::exp(-reaction->forward_ArrheniusEnergy / 8.314 / y[1]);	// [mol/m3/Pa/s]
			r = k * y[0] / totalMolarFlow[0] * P[0];														// [mol/m3/s]
			R += reaction->stoichiometry[base] * r;															// [mol/m3/s]
			Qr = reaction->heatOfReaction * r;																// [J/m3/s]
		}
		dydz[0] = A * R;																					// [mol/s/m]
		dydz[1] = A * Qr / (Cp[0] * totalMolarFlow[0]);														// [K/m]
	}

	void odeSolver() {

		// Set initial conditions
		state_type y(2);
		y[0] = molarFlow[0];
		y[1] = T[0];

		// Define step size
		const CapeReal dz = 0.01;

		// Integrate with constant step
		// Perfromance optimisation ignored for Proof of Concept
		// integrate_const(runge_kutta4<state_type>(), &Solver::pfr, y, 0.0, this->L, dz);

		// Post processing
		printf("%f\n%f", y[0], y[1]);
		molarFlow[1] += molarFlow[0] - y[0];
		molarFlow[2] += molarFlow[0] - y[0];
		molarFlow[0] = y[0];
		T[0] = y[1];
	}

	void setProducts() {

		// Set product(s) overall props.
		product.SetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		product.SetOverallProp(ConstCapeString(COBIATEXT("fraction")),
			ConstCapeString(COBIATEXT("mole")), X);
		product.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		product.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), P);
	}
	
	void flashProducts() {

	}


	//CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = L"Solver";
	}
	
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = L"ODE Solver for mass and energy balance";
	}
	
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}

};

using SolverPtr = CapeOpenObjectSmartPointer<Solver>;
