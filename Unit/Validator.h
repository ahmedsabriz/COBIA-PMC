#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "PortCollection.h"

using namespace COBIA;


class Validator :

	/*
	* Validator Class
	* TODO: Description
	*/

	public CapeOpenObject<Validator>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Validator> {
	
	// Members
	CAPEOPEN_1_2::CapeCollection<MaterialPort> portCollection;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Validator");
	}

	Validator(CAPEOPEN_1_2::CapeCollection<MaterialPort> _portCollection) : portCollection(_portCollection) {
	}

	~Validator() {
	}

	// Validate ICapeParameterSpecification & ICapeRealParameterSpecification
	// TODO: Iterate over parameter collection instead
	/*if (param1->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
		val = param1->Validate(message) && param1->Validate(param1->getValue(), message);*/


	CapeBoolean validateMaterialPorts(/*out*/ CapeString message) {

		// Check whether all ports are connected, and connected to materials with equal compound lists
		CapeInterface connectedObject;
		CapeString portName;
		CapeArrayStringImpl refCompIDs, compIDs;
		CapeArrayStringImpl formulae, names, casNumbers;
		CapeArrayRealImpl boilTemps, molecularWeights;
		CapeBoolean sameCompList;

		for (MaterialPort p : portCollection) {
			connectedObject = p.getConnectedObject();

			// Check whether port is connected
			if (!connectedObject) {
				if (p.isPrimary()) {
					p.getComponentName(portName);
					CapeStringImpl _portName(portName);
					message = COBIATEXT("Port ") + _portName + COBIATEXT(" is not connected");
					return false;
				}
				else {
					p.ignoreOptionalStream();
				}
			}

			if (p.isConnected() && p.getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL) {
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
				else if (refCompIDs.empty()) {
					message = COBIATEXT("Connected material streams expose zero compound");
					return false;
				}
			}
		}
		return true;
	}
	// TODO (Unit Specific)
	// Validate that each input has output
	// Validate a minimum of 1 hot and 1 cold streams

	void preparePhaseIDs (/*out*/ std::vector<CapeArrayStringImpl> productsPhaseIDs,
		/*out*/ std::vector<CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus>> productsPhaseStatus) {

		// Prepare lists of supported phase label and phase status for product flash
		// This remains constant between validations
		CAPEOPEN_1_2::CapeThermoMaterial material;
		CapeArrayStringImpl phaseIDs, stateOfAggregation, keyCompounds;
		CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;

		for (MaterialPort p : portCollection) {
			if (p.isConnected() && p.getPortType() == CAPEOPEN_1_2::CAPE_OUTLET) {
				material = p.getMaterial();
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
