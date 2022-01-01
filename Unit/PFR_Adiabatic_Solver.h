#pragma once
#include <COBIA.h>
#include <boost/numeric/odeint.hpp>
#include "ReactionPackage.h"
#include "Reaction.h"


using namespace COBIA;


class Solver :

	/*
		This solver uses boost's numeric library to solve the defined pfr
		Class is intiated with a single stream (product) and a reaction package
		Product stream must be identical to the PFR's single feed stream
	*/

	public CapeOpenObject<Solver>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Solver> {

	// Members
	MaterialPortPtr productStream;
	CAPEOPEN_1_2::CapeThermoMaterial product;

	std::vector<ReactionPtr> reactions;
	const CapeInteger numComponents, numReactions;

	ParameterCollectionPtr paramCollection;
	CAPEOPEN_1_2::CapeRealParameter reactorLength, reactorVolume;

	CapeArrayRealImpl molarFlow, totalMolarFlow, T, P, X;
	CapeReal Cp, heatFlow, inletBaseMolarFlow;

	// The type of container used to hold the state vector to solve odeint's runge_kutta
	typedef std::vector<CapeReal> state_type;

public:

	std::vector<state_type> profile;
	state_type grid;

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Solver");
	}

	Solver(MaterialPortPtr _productStream, ReactionPackagePtr _reactionPackage,
		ParameterCollectionPtr _paramCollection) :
		productStream(_productStream), product(productStream->getMaterial()),
		reactions(_reactionPackage->getReactions()),
		numComponents(_reactionPackage->numComponents),
		numReactions(_reactionPackage->numReactions),
		paramCollection(_paramCollection),
		reactorLength(paramCollection->RealItem(0)), 
		reactorVolume(paramCollection->RealItem(1)) {
		getInitialConditions();
		odeSolver();
	}

	~Solver() {
	}

	void pfr(state_type& y, state_type& dydz, CapeReal z) {
		/* Assumptions:
		*	1D
		*	Homogenous
		*	Empty bed
		*	Adiabatic
		*	Constant heat capacity
		*	No axial dispersion
		*	No radial gradient
		*	No pressure drop
		*/

		// The following variables are overwritten in every step
		CapeReal n_tot, p, Treact, k, r, Qr;
		CapeArrayRealImpl R(numComponents);
		n_tot = 0;
		for (size_t i = 0; i < numComponents; i++) {
			n_tot += y[i];													// [mol/s]
			R[i] = 0;
		}
		p = y[0] / n_tot * P[0];											// [Pa]
		Treact = y[numComponents];											// [K]

		for (ReactionPtr reaction : reactions) {
			k = reaction->forwarRateConstant * std::exp(-reaction->forward_ArrheniusEnergy
				/ 8.31446261815324 / Treact);								// [mol/m3/Pa/s]
			r = k * p;														// [mol/m3/s]
			for (size_t i = 0; i < numComponents; i++) {
				R[i] += reaction->stoichiometry[i] * r;						// [mol/m3/s]
				dydz[i] = reactorVolume.getValue() / reactorLength.getValue() * R[i];
			}
			// To Do: Calculate Qr as function of T
			Qr = reaction->heatOfReaction * R[0];							// [J/m3/s]
		}
		dydz[numComponents] = reactorVolume.getValue() / reactorLength.getValue() * Qr
			/ (Cp * totalMolarFlow[0]);										// [K/m]

		// TODO: Set product after each iteration to obtain accurate property
	}

	// In case of variation in HOR
	CapeReal heatOfReaction(CapeReal T) {

		CapeInterface connectedObject = productStream->getConnectedObject();
		CAPEOPEN_1_2::CapeThermoCompounds compounds(connectedObject);

		CapeArrayStringImpl compIDs, formulae, names, casNumbers;
		CapeArrayRealImpl boilTemps, molecularWeights;
		compounds.GetCompoundList(compIDs, formulae, names,
			boilTemps, molecularWeights, casNumbers);

		CapeArrayRealImpl enthalpies;
		CapeArrayStringImpl props(1);
		props[0] = COBIATEXT("idealGasEnthalpy");
		CapeBoolean containsMissingValues = true;

		compounds.GetTDependentProperty(props, T, compIDs, containsMissingValues, enthalpies);
		return enthalpies[1] + enthalpies[2] - enthalpies[0];
	}

	void getInitialConditions() {

		// Get inlet molarFlow
		product.GetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		inletBaseMolarFlow = molarFlow[0];

		product.GetOverallProp(ConstCapeString(COBIATEXT("totalflow")),
			ConstCapeString(COBIATEXT("mole")), totalMolarFlow);

		// Get inlet T[K], P[Pa], X[mol/mol]
		T.resize(1);
		P.resize(1);
		product.GetOverallTPFraction(T[0], P[0], X);

		Cp = calcOverallFromPhaseProps(product, COBIATEXT("heatCapacityCp"));	// [J/mol/K]
		heatFlow = calcOverallFromPhaseProps(product, COBIATEXT("enthalpyF")) *
			totalMolarFlow[0];													// [J/mol]
	}

	state_type odeSolver() {
		
		// Set initial conditions
		state_type y(numComponents + 1);
		for (size_t i = 0; i < numComponents; i++) {
			y[i] = molarFlow[i];												// [mol/s]
		}
		y[numComponents] = T[0];

		// Define step size
		CapeReal dz = 0.01;

		// Observer definition
		struct observer {
			std::vector<state_type>& states;
			std::vector<double>& times;

			observer(std::vector<state_type>& _states, state_type& _times)
				: states(_states), times(_times) { }

			void operator()(const state_type& y, CapeReal t)
			{
				states.push_back(y);
				times.push_back(t);
			}
		};

		// Integrate with constant step
		// Perfromance optimisation ignored for Proof of Concept
		boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta4<state_type>(),
			std::bind(&Solver::pfr, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3),
			y, 0.0, this->reactorLength.getValue(), dz, observer(profile, grid));
		
		setProduct(y, P);
		return y;
	}

	void setProduct(/*in*/ state_type y, /*in*/ CapeArrayRealImpl _P) {

		// Post processing
		totalMolarFlow[0] = 0;
		for (size_t i = 0; i < numComponents; i++) {
			molarFlow[i] = y[i];												// [mol/s]
			totalMolarFlow[0] += molarFlow[i];									// [mol/s]
		}
		T[0] = y[numComponents];

		// Set product(s) overall props.
		product.SetOverallProp(ConstCapeString(COBIATEXT("flow")),
			ConstCapeString(COBIATEXT("mole")), molarFlow);
		product.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		product.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), _P);
	}

	void flashProduct(/*out*/ CapeArrayStringImpl _phaseIDs,
		/*out*/ CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> _phaseStatus) {

		// Allow all phases to take part in product flash
		product.SetPresentPhases(_phaseIDs, _phaseStatus);

		// Flash product(s) at specified T & P
		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine equilibriumRoutine(product);

		// Prepare T & P flash specifications for product flash
		// specification format: CapeArrayRealImpl = { propertyIdentifier, basis, phaseLabel [, compoundIdentifier] }
		// basis is undefined when it is not a dependency of the property (e.g. T, P)
		CapeArrayStringImpl flashCondT(3), flashCondP(3);
		flashCondT[0] = COBIATEXT("temperature");
		flashCondT[2] = COBIATEXT("overall");
		flashCondP[0] = COBIATEXT("pressure");
		flashCondP[2] = COBIATEXT("overall");

		equilibriumRoutine.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());
	}

	CapeReal dH() {
		return heatFlow - calcOverallFromPhaseProps(product, COBIATEXT("enthalpyF")) *
			totalMolarFlow[0];												// [J/mol]
	}

	CapeReal calConversion() {
		return (inletBaseMolarFlow - molarFlow[0]) / inletBaseMolarFlow;	// [J/mol]
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
