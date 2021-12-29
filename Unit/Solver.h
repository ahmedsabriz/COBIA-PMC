#pragma once
#include <COBIA.h>
#include "PortCollection.h"
#include "MaterialPort.h"

using namespace COBIA;

class Solver :

	/*
		TODO
	*/

	public CapeOpenObject<Solver>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Solver> {

	// Members
	PortCollectionPtr& portCollection;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Solver");
	}

	Solver(PortCollectionPtr& _portCollection) : portCollection(_portCollection) {

		// CAPE-OPEN unit operations may not have side effects on material objects connected to feeds.
		// Product material object id copied from feed material object allowing to perfrom calculations
		// before setting its properties to override feed peoperties' values
		MaterialPortPtr portPtr, inletPtr;
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			if (portPtr->isConnected() && portPtr->getDirection() == CAPEOPEN_1_2::CAPE_OUTLET) {
				inletPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i-1));
				portPtr->getMaterial().CopyFromMaterial(inletPtr->getMaterial());
			}
		}
	}

	~Solver() {
	}
	
	void getInitialConditions() {
	}

	void setProduct() {
		// Set product(s) overall props.
	}

	void flashProduct(/*in*/ std::vector<CapeArrayStringImpl> phaseIDs,
		/*in*/ std::vector<CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus>> phaseStatus,
		/*in*/ CapeArrayStringImpl flashCond1, /*in*/ CapeArrayStringImpl flashCond2) {

		size_t j = 0;
		MaterialPortPtr portPtr;
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			if (portPtr->isConnected() &&
				portPtr->getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL &&
				portPtr->getDirection() == CAPEOPEN_1_2::CAPE_OUTLET) {

				// Allow all phases to take part in product flash
				CAPEOPEN_1_2::CapeThermoMaterial material = portPtr->getMaterial();
				material.SetPresentPhases(phaseIDs[j], phaseStatus[j]);

				// Flash product(s) at specified Conditions
				CAPEOPEN_1_2::CapeThermoEquilibriumRoutine equilibriumRoutine(material);
				equilibriumRoutine.CalcEquilibrium(flashCond1, flashCond2, ConstCapeEmptyString());

				j++;
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