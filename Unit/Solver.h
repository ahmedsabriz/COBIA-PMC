#pragma once
#include <COBIA.h>
#include <boost/numeric/odeint.hpp>
#include "ReactionPackage.h"
#include "Reaction.h"

using namespace boost::numeric::odeint;
using namespace COBIA;

class Solver :
	public CapeOpenObject<Solver>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Solver> {

	// Members
	CAPEOPEN_1_2::CapeThermoMaterial product;
	std::vector<ReactionPtr> reactions;

	const CapeInteger numComponents, numReactions;

	CapeArrayRealImpl molarFlow, totalMolarFlow, T, P, X, Cp;

	// The type of container used to hold the state vector
	typedef std::vector<CapeReal> state_type;// CapeArrayRealImpl state_type;

	// TODO: REMOVE HARDCODING
	const CapeReal L = 6;			// [m]
	const CapeReal A = 0.77 / L;	// [m2]

	// Focus on single component (Ethyl benzene at index 0] for Proof Of Concept
	size_t base = 0;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Solver");
	}

	Solver(CAPEOPEN_1_2::CapeThermoMaterial _product, ReactionPackagePtr _reactionPackage) :
		product(_product), reactions(_reactionPackage->getReactions()),
		numComponents(_reactionPackage->numComponents), numReactions(_reactionPackage->numReactions) {
	}

	~Solver() {
	}

	void getInitialConditions() {

		// Get inlet molarFlow
		product.GetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		product.GetOverallProp(ConstCapeString(COBIATEXT("totalflow")),
			ConstCapeString(COBIATEXT("mole")), totalMolarFlow);

		// Get inlet T[K], P[Pa], X[mol/mol]
		T.resize(1);
		P.resize(1);
		product.GetOverallTPFraction(T[0], P[0], X);

		// Get present phaseLabels to calculate their properties
		CapeArrayStringImpl phaseLabels;
		CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;
		product.GetPresentPhases(phaseLabels, phaseStatus);

		// Before getting a property from a material object, property must be calculated using CalcSinglePhaseProp.
		// Exceptions are for properties that are available at phase equilibrium:
		// flow rate, composition, pressure, temperature and phase fraction
		CapeArrayStringImpl prop(1);

		// Get heat capacity at constant pressure [J/mol/K]
		prop[0] = L"heatCapacityCp";
		Cp.resize(phaseLabels.size());

		// Iterate over present phases and calculate selected phase propertie
		for (ConstCapeString phaseLabel : phaseLabels) {
			CAPEOPEN_1_2::CapeThermoPropertyRoutine routine(product);
			routine.CalcSinglePhaseProp(prop, phaseLabel);

			// Get and store the properties in an array.
			CapeArrayRealImpl phaseCp;
			product.GetSinglePhaseProp(prop[0],	phaseLabel, ConstCapeString(COBIATEXT("mole")), phaseCp);
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
		
		// The following variables are overwritten in every step
		CapeReal n_tot, p, Treact, k, r, Qr;
		CapeArrayRealImpl R(numComponents);
		n_tot = 0;
		for (size_t i = 0; i < numComponents; i++) {
			n_tot += y[i];																						// [mol/s]
			R[i] = 0;
		}
		p = y[base] / n_tot * P[0];																				// [Pa]
		Treact = y[numComponents];																				// [K]

		for (ReactionPtr reaction : reactions) {
			// TODO: Implement Universal Constants
			k = reaction->forwarRateConstant * std::exp(-reaction->forward_ArrheniusEnergy / 8.314 / Treact);	// [mol/m3/Pa/s]
			r = k * p;																							// [mol/m3/s]
			for (size_t i = 0; i < numComponents; i++) {
				R[i] += reaction->stoichiometry[i] * r;															// [mol/m3/s]
				dydz[i] = A * R[i];
			}
			Qr = reaction->heatOfReaction * R[base];															// [J/m3/s]
		}
		dydz[numComponents] = A * Qr / (Cp[1] * totalMolarFlow[0]);												// [K/m]
	}

	void odeSolver() {
		
		// Set initial conditions
		state_type y(numComponents + 1);
		for (size_t i = 0; i < numComponents; i++) {
			y[i] = molarFlow[i];																				// [mol/s]
		}
		y[numComponents] = T[0];

		// Define step size
		const CapeReal dz = 0.01;

		// Integrate with constant step
		// Perfromance optimisation ignored for Proof of Concept
		using namespace std::placeholders;
		integrate_const(runge_kutta4<state_type>(),
			std::bind(&Solver::pfr, this, _1, _2, _3),
			y, 0.0, this->L, dz);

		// Post processing
		for (size_t i = 0; i < numComponents; i++) {
			molarFlow[i] = y[i];																				// [mol/s]
		}
		T[0] = y[numComponents];
	}

	void setProducts() {

		// Set product(s) overall props.
		product.SetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		product.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		product.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), P);
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
