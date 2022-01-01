#pragma once
#include <COBIA.h>
#include "PortCollection.h"
#include "MaterialPort.h"

using namespace COBIA;


class Validator :

	/*
	* Validator Class
	* TODO: Description
	* Validate ICapeParameterSpecification & ICapeRealParameterSpecification
	* if (param1->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
	*	val = param1->Validate(message) && param1->Validate(param1->getValue(), message);
	* }
	*/

	public CapeOpenObject<Validator>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Validator> {
	
	// Members
	PortCollectionPtr& portCollection;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Validator");
	}

	Validator(PortCollectionPtr& _portCollection) : portCollection(_portCollection) {
	}

	~Validator() {
	}

	CapeBoolean validateMaterialPorts(/*out*/ CapeString message) {

		// Check whether all ports are connected, and connected to materials with equal compound lists
		CapeInterface connectedObject;
		CapeString portName;
		CapeArrayStringImpl refCompIDs, compIDs;
		CapeArrayStringImpl formulae, names, casNumbers;
		CapeArrayRealImpl boilTemps, molecularWeights;
		CapeBoolean sameCompList;

		MaterialPortPtr portPtr;
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			connectedObject = portPtr->getConnectedObject();
			// Check whether port is connected
			if (!connectedObject) {
				if (portPtr->isPrimary()) {
					portPtr->getComponentName(portName);
					CapeStringImpl _portName(portName);
					message = COBIATEXT("Port ") + _portName + COBIATEXT(" is not connected");
					return false;
				}
			}
			else {
				if (portPtr->getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL) {
					CAPEOPEN_1_2::CapeThermoCompounds compounds(connectedObject);
					// Get compound list for from first stream as reference
					if (refCompIDs.empty()) {
						compounds.GetCompoundList(refCompIDs, formulae, names,
							boilTemps, molecularWeights, casNumbers);
					}
					// Compare all other streams with it
					compounds.GetCompoundList(compIDs, formulae, names, boilTemps, molecularWeights, casNumbers);
					sameCompList = (refCompIDs.size() == compIDs.size());
					for (size_t i = 0; (i < compIDs.size()) && (sameCompList); i++) {
						sameCompList = (compIDs[i].compare(refCompIDs[i]) == 0);
					}
					if (!sameCompList) {
						message = COBIATEXT("Connected material streams expose inconsistent compound lists");
						return false;
					}
					else if (compIDs.empty()) {
						message = COBIATEXT("Connected material streams expose zero compound");
						return false;
					}
				}
			}
		}
		return true;
	}

	CapeBoolean validateHeatExchangerInletOutlet(/*out*/ CapeString message) {

		// 	Validate that each optional inlet has an outlet
		MaterialPortPtr portPtr;
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			if (
				!portPtr->isPrimary() &&
				portPtr->getConnectedObject() &&
				portPtr->getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL &&
				portPtr->getDirection() == CAPEOPEN_1_2::CAPE_INLET
				) {
				if (!portCollection->Item(i + 1).getConnectedObject()) {
					message = COBIATEXT("An inlet material stream is missing an outlet connection");
					return false;
				}
			}
		}

		// 	Validate that each optional outlet has an inlet
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			if (
				!portPtr->isPrimary() &&
				portPtr->getConnectedObject() &&
				portPtr->getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL &&
				portPtr->getDirection() == CAPEOPEN_1_2::CAPE_OUTLET
				) {
				if (!portCollection->Item(i - 1).getConnectedObject()) {
					message = COBIATEXT("An outlet material stream is missing an inlet connection");
					return false;
				}
			}
		}
		
		// TODO: Validate a minimum of 1 hot and 1 cold streams

		return true;
	}

	void preparePhaseIDs (/*out*/ std::vector<CapeArrayStringImpl> &productsPhaseIDs,
		/*out*/ std::vector<CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus>> &productsPhaseStatus) {

		// Clear vectors before requesting phaseIDs of products
		productsPhaseIDs.clear();
		productsPhaseStatus.clear();
		// Prepare lists of supported phase label and phase status for product flash
		// This remains constant between validations
		CAPEOPEN_1_2::CapeThermoMaterial material;
		CapeArrayStringImpl phaseIDs, stateOfAggregation, keyCompounds;
		CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;

		MaterialPortPtr portPtr;
		for (size_t i = 0, count = portCollection->getCount(); i < count; i++) {
			portPtr = static_cast<MaterialPort*>((CAPEOPEN_1_2::ICapeUnitPort*)portCollection->Item(i));
			if (portPtr->getConnectedObject() && portPtr->getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL &&
				portPtr->getDirection() == CAPEOPEN_1_2::CAPE_OUTLET) {
				material = portPtr->getMaterial();
				CAPEOPEN_1_2::CapeThermoPhases phases(material);
				phases.GetPhaseList(phaseIDs, stateOfAggregation, keyCompounds);
				phaseStatus.resize(phaseIDs.size());
				std::fill(phaseStatus.begin(), phaseStatus.end(), CAPEOPEN_1_2::CAPE_UNKNOWNPHASESTATUS);
				productsPhaseIDs.emplace_back(phaseIDs);
				productsPhaseStatus.emplace_back(phaseStatus);
			}
		}
	}

	//CAPEOPEN_1_2::ICapeIdentification
	void getComponentName(/*out*/ CapeString name) {
		name = COBIATEXT("Validator");
	}
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Validator Class");
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
};

using ValidatorPtr = CapeOpenObjectSmartPointer<Validator>;