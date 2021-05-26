#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "PortCollection.h"

using namespace COBIA;

class Unit :
	public CapeOpenObject<Unit>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUnitAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUtilitiesAdapter<Unit>,
	public CAPEOPEN_1_2::CapePersistAdapter<Unit> {

	// Members
	CapeStringImpl name, description;
	CapeBoolean dirty;
	MaterialPortPtr feed, product;
	PortCollectionPtr portCollection;
	CAPEOPEN_1_2::CapeValidationStatus validationStatus;
	CapeArrayStringImpl flashCondT, flashCondP;
	CapeArrayStringImpl flashPhaseIDs;
	CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> flashPhaseStatus;
public:

	const CapeStringImpl getDescriptionForErrorSource() {
		// Returns a description of the current object for error handling
		// May be changed to any type (by value of reference) that implements const CapeCharacter * c_str()
		return COBIATEXT("Unit") + name;
	}

	//Registration info
	static const CapeUUID getObjectUUID() {
		//Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704A8
		return CapeUUID{{0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xa8}};
	}

	static void Register(CapePMCRegistrar registrar) {
		registrar.putName(COBIATEXT("Unit Name"));
		registrar.putDescription(COBIATEXT("Unit Description"));
		registrar.putCapeVersion(COBIATEXT("1.1"));
		registrar.putComponentVersion(COBIATEXT("0.1.0.0"));
		registrar.putAbout(COBIATEXT("About Unit"));
		registrar.putVendorURL(COBIATEXT("www.unitwebsite.com"));
		registrar.putProgId(COBIATEXT("Polimi.Unit"));
		registrar.addCatID(CAPEOPEN::categoryId_UnitOperation);
		registrar.addCatID(CAPEOPEN_1_2::categoryId_Component_1_2);
		//registrar.putCreationFlags(CapePMCRegistationFlag_None);
	}

	Unit()  :
		name(COBIATEXT("Unit")),
		feed(new MaterialPort(name, validationStatus, CAPEOPEN_1_2::CAPE_INLET, COBIATEXT("Feed"))),
		product(new MaterialPort(name, validationStatus, CAPEOPEN_1_2::CAPE_OUTLET, COBIATEXT("Product"))),
		portCollection(new PortCollection(name)),
		validationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED),
		dirty(false) {
		portCollection->addPort(feed);
		portCollection->addPort(product);
		
		// Prepare T & P flash conditions for product flash
		flashCondT.resize(3);
		flashCondT[0] = COBIATEXT("temperature");
		flashCondT[2] = COBIATEXT("overall");
		flashCondP.resize(3);
		flashCondP[0] = COBIATEXT("pressure");
		flashCondP[2] = COBIATEXT("overall");
	}

	~Unit() {
	}

	// CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->name;
	}
	void putComponentName(/*in*/ CapeString name) {
		this->name = name;
		dirty = true;
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = description;
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		description = desc;
		dirty = true;
	}

	// CAPEOPEN_1_2::ICapeUnit
	
	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> ports() {
		return portCollection;
	}
	CAPEOPEN_1_2::CapeValidationStatus getValStatus() {
		return validationStatus;
	}
	void Calculate() {

		// Check validation status
		if (validationStatus != CAPEOPEN_1_2::CAPE_VALID) {
			//the PME must validate the unit prior to calculation
			throw cape_open_error(COBIATEXT("Unit is not in a valid state"));
		}

		// Import material streams
		CAPEOPEN_1_2::CapeThermoMaterial feedMaterial = feed->getMaterial();
		CAPEOPEN_1_2::CapeThermoMaterial productMaterial = product->getMaterial();

		// This array of property values could be a class member; that would avoid most memory allocations
		// the array is re-used below for various properties.
		CapeArrayRealImpl propValues;

		// Copy total flow
		feedMaterial.GetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), propValues);
		productMaterial.SetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), propValues);

		// Copy P, X and T
		CapeReal T, P;
		feedMaterial.GetOverallTPFraction(T, P, propValues); // X = Propvalues
		productMaterial.SetOverallProp(ConstCapeString(COBIATEXT("fraction")),
			ConstCapeString(COBIATEXT("mole")), propValues);
		propValues.resize(1);
		propValues[0] = P;
		productMaterial.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), propValues);
		propValues[0] = T;
		productMaterial.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), propValues);

		// Needs more clarification
		// Allow all phases to take part in product flash
		productMaterial.SetPresentPhases(flashPhaseIDs, flashPhaseStatus);

		// Flash product at specified T & P
		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine productEquilibrium(productMaterial);
		productEquilibrium.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean res = true;

		// Check whether all ports are connected, and connected to materials with equal compound lists
		CapeArrayStringImpl refCompIDs, compIDs;
		CapeArrayStringImpl formulae, names, casNumbers;
		CapeArrayRealImpl boilTemps, molecularWeights;
		CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> collection(portCollection);
		for (CAPEOPEN_1_2::CapeUnitPort p : collection) {
			CapeInterface connectedObject = p.getConnectedObject();
			
			// Check whether port is connected
			if (!connectedObject) {
				CAPEOPEN_1_2::CapeIdentification portIdentification(p);
				CapeStringImpl portName;
				portIdentification.getComponentName(portName);
				message = COBIATEXT("Port ") + portName + COBIATEXT(" is not connected");
				res = false;
				break;
			}

			// Needs more clarification
			CAPEOPEN_1_2::CapeThermoCompounds materialCompounds(connectedObject);
			if (refCompIDs.empty()) {
				materialCompounds.GetCompoundList(refCompIDs, formulae, names,
					boilTemps, molecularWeights, casNumbers);
			}
			else {
				materialCompounds.GetCompoundList(compIDs, formulae, names,
					boilTemps, molecularWeights, casNumbers);
			}
		}
		if (res) {

			// Additional checks
			bool sameCompList = (refCompIDs.size() == compIDs.size());
			for (size_t i = 0; (i < compIDs.size()) && (sameCompList); i++) {
				sameCompList = (compIDs[i].compare(refCompIDs[i]) == 0);
			}
			if (!sameCompList) {
				message = COBIATEXT("Connected material streams expose inconsistent compound lists");
				res = false;
			}
			else if (refCompIDs.empty()) {
				message = COBIATEXT("Connected material streams expose zero compound");
				res = false;
			}
		}
		if (res) {
			// Prepare the flash phase list for the product flash
			// This remains constant between validations
			CAPEOPEN_1_2::CapeThermoMaterial material = product->getMaterial();
			CAPEOPEN_1_2::CapeThermoPhases materialPhases(material);
			CapeArrayStringImpl stateOfAggregation, keyCompounds;
			materialPhases.GetPhaseList(flashPhaseIDs, stateOfAggregation, keyCompounds);
			flashPhaseStatus.resize(flashPhaseIDs.size());
			std::fill(flashPhaseStatus.begin(), flashPhaseStatus.end(), CAPEOPEN_1_2::CAPE_UNKNOWNPHASESTATUS);
		}

		// Update validation status
		validationStatus = (res) ? CAPEOPEN_1_2::CAPE_VALID : CAPEOPEN_1_2::CAPE_INVALID;
		return res;
	}

	// CAPEOPEN_1_2::ICapeUtilities

	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeParameter> getParameters() {
		throw cape_open_error(COBIAERR_NotImplemented); // No parameters
	}
	void putSimulationContext(/*in*/ CAPEOPEN_1_2::CapeSimulationContext context) {
		// Do nothing
	}
	void Initialize() {
		// Do nothing
	}
	void Terminate() {
		// Disconnect ports
		// In case a reference to the simulation context is stored, it too must be released at Terminate
		CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> collection(portCollection);
		for (CAPEOPEN_1_2::CapeUnitPort p : collection) {
			p.Disconnect();
		}
	}
	CAPEOPEN_1_2::CapeEditResult Edit(CapeWindowId parent) {
		throw cape_open_error(COBIAERR_NotImplemented); //no GUI
	}

	//CAPEOPEN_1_2::ICapePersist

	void Save(/*in*/ CAPEOPEN_1_2::CapePersistWriter writer,/*in*/ CapeBoolean clearDirty) {
		writer.Add(ConstCapeString(COBIATEXT("name")), name);
		writer.Add(ConstCapeString(COBIATEXT("description")), description);
		if (clearDirty) {
			dirty = false;
		}
	}

	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};