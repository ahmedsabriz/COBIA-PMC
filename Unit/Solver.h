#pragma once
#include <COBIA.h>
#include <boost/numeric/odeint.hpp>
#include "ParameterCollection.h"
#include "MaterialPort.h"


using namespace boost::numeric::odeint;
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
	CAPEOPEN_1_2::CapeCollection<MaterialPort> portCollection;
	ParameterCollectionPtr paramCollection;

	/*CapeArrayRealImpl molarFlow, totalMolarFlow, T, P, X;
	CapeReal Cp, heatFlow, inletBaseMolarFlow;*/

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Solver");
	}

	Solver(CAPEOPEN_1_2::CapeCollection<MaterialPort> _portCollection,
		ParameterCollectionPtr _paramCollection) :
		portCollection(_portCollection),
		paramCollection(_paramCollection) {

		// CAPE-OPEN unit operations may not have side effects on material objects connected to feeds.
		// Product material object id copied from feed material object allowing to perfrom calculations
		// before setting its properties to override feed peoperties' values
		
		for(MaterialPort p : portCollection) {
			if (p.isConnected() && p.getDirection() == CAPEOPEN_1_2::CAPE_OUTLET) {
				CapeString portName;
				p.getComponentName(portName);	
				p.getMaterial().CopyFromMaterial(portCollection[COBIATEXT("Input " + portName.c_str()[portName.length() - 1])].getMaterial());
			}
		}
	}

	~Solver() {
	}
	
	void getInitialConditions() {

		// Get inlet molarFlow
		//product.GetOverallProp(ConstCapeString(COBIATEXT("flow")),
		//	ConstCapeString(COBIATEXT("mole")), molarFlow);
		//inletBaseMolarFlow = molarFlow[base];

		//product.GetOverallProp(ConstCapeString(COBIATEXT("totalflow")),
		//	ConstCapeString(COBIATEXT("mole")), totalMolarFlow);

		//// Get inlet T[K], P[Pa], X[mol/mol]
		//T.resize(1);
		//P.resize(1);
		//product.GetOverallTPFraction(T[0], P[0], X);

		//Cp = calcOverallFromPhaseProps(product, COBIATEXT("heatCapacityCp"));	// [J/mol/K]
		//heatFlow = calcOverallFromPhaseProps(product, COBIATEXT("enthalpyF")) *
		//	totalMolarFlow[0];													// [J/mol]
	}

	void setProduct() {
		// Set product(s) overall props.
	}

	void flashProduct(/*in*/ std::vector<CapeArrayStringImpl> phaseIDs,
		/*in*/ std::vector<CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus>> phaseStatus,
		/*in*/ CapeArrayStringImpl flashCond1, /*in*/ CapeArrayStringImpl flashCond2) {

		size_t i = 0;
		for (MaterialPort p : portCollection) {
			if (p.isConnected() && p.getDirection() == CAPEOPEN_1_2::CAPE_OUTLET) {

				// Allow all phases to take part in product flash
				CAPEOPEN_1_2::CapeThermoMaterial material = p.getMaterial();
				material.SetPresentPhases(phaseIDs[i], phaseStatus[i]);

				// Flash product(s) at specified Conditions
				CAPEOPEN_1_2::CapeThermoEquilibriumRoutine equilibriumRoutine(material);
				equilibriumRoutine.CalcEquilibrium(flashCond1, flashCond2, ConstCapeEmptyString());

				i++;
			}
		}
	}

	//CAPEOPEN_1_2::ICapeIdentification
	void getComponentName(/*out*/ CapeString name) {
		name = COBIATEXT("Solver");
	}
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("ODE Solver for mass and energy balance");
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}

};

using SolverPtr = CapeOpenObjectSmartPointer<Solver>;
