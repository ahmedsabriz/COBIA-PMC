#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "PortCollection.h"
#include "RealParameter.h"
#include "ParameterCollection.h"

using namespace COBIA;

class Unit :
	public CapeOpenObject<Unit>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUnitAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUtilitiesAdapter<Unit>,
	public CAPEOPEN_1_2::CapePersistAdapter<Unit> {

	// Members
	CapeStringImpl name, description;

	// Ports and Parameters
	MaterialPortPtr feed, product1, product2;
	RealParameterPtr splitRatio;

	// Collections
	PortCollectionPtr portCollection;
	ParameterCollectionPtr paramCollection;

	// Material Properties
	CapeArrayRealImpl feedFlow, p1Flow, p2Flow;
	CapeArrayRealImpl T, P, X;

	// Thermo Properties
	CapeArrayStringImpl flashCondT, flashCondP;
	CapeArrayStringImpl flashPhaseIDs;
	CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> flashPhaseStatus;

	// Persistence and Validation
	CapeBoolean dirty;
	CAPEOPEN_1_2::CapeValidationStatus validationStatus;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		// Returns a description of the current object for error handling
		return COBIATEXT("Unit: ") + *&name;
	}

	// Registration info
	static const CapeUUID getObjectUUID() {
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704A8
		return CapeUUID{{0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xa8}};
	}

	static void Register(CapePMCRegistrar registrar) {
		registrar.putName(COBIATEXT("CO Splitter"));
		registrar.putDescription(COBIATEXT("Description"));
		registrar.putCapeVersion(COBIATEXT("1.1"));
		registrar.putComponentVersion(COBIATEXT("0.2.0.1"));
		registrar.putAbout(COBIATEXT("About Unit"));
		registrar.putVendorURL(COBIATEXT("www.polimi.it"));
		registrar.putProgId(COBIATEXT("Polimi.Unit"));
		registrar.addCatID(CAPEOPEN::categoryId_UnitOperation);
		registrar.addCatID(CAPEOPEN_1_2::categoryId_Component_1_2);
		//registrar.putCreationFlags(CapePMCRegistationFlag_None);
	}

	Unit()  :
		name(COBIATEXT("CO Splitter")), description(COBIATEXT("Description")),
		feed(new MaterialPort(COBIATEXT("Feed"), CAPEOPEN_1_2::CAPE_INLET, name, validationStatus)),
		product1(new MaterialPort(COBIATEXT("Product 1"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		product2(new MaterialPort(COBIATEXT("Product 2"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		portCollection(new PortCollection(name)), 
		splitRatio(new RealParameter(name, COBIATEXT("Split Ratio"), 0.5, 0.0, 1.0, validationStatus, dirty)),
		paramCollection(new ParameterCollection(name)), 
		validationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED),
		dirty(false) {

		// Set parameter dimensionality
		splitRatio->putDimensionality(8, 1);

		// Add ports and parameters to collections
		portCollection->addPort(feed);
		portCollection->addPort(product1);
		portCollection->addPort(product2);
		paramCollection->addParameter(splitRatio);
		
		// Prepare T & P flash conditions for product flash
		flashCondT.resize(3);
		flashCondT[0] = COBIATEXT("temperature");
		flashCondT[2] = COBIATEXT("overall");
		flashCondP.resize(3);
		flashCondP[0] = COBIATEXT("pressure");
		flashCondP[2] = COBIATEXT("overall");
	}

	~Unit() {
		delete(paramCollection);
		delete(splitRatio);
		delete(portCollection);
		delete(product1);
		delete(product2);
		delete(feed);
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
		
		// Check validation status before calculation
		if (validationStatus != CAPEOPEN_1_2::CAPE_VALID) {
			throw cape_open_error(COBIATEXT("Unit is not in a valid state"));
		}

		// Implement ICapeThermoMaterial Interface
		CAPEOPEN_1_2::CapeThermoMaterial feedMaterial = feed->getMaterial();
		CAPEOPEN_1_2::CapeThermoMaterial p1Material = product1->getMaterial();
		CAPEOPEN_1_2::CapeThermoMaterial p2Material = product2->getMaterial();
		
		// Set product(s) flowrate
		feedMaterial.GetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), feedFlow);

		p1Flow.resize(1);
		p2Flow.resize(1);

		p1Flow[0] = feedFlow[0] * splitRatio->getValue();
		p1Material.SetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), p1Flow);

		p2Flow[0] = feedFlow[0] - p1Flow[0];
		p1Material.SetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), p1Flow);

		// Set product(s) T, P, X
		T.resize(1);
		P.resize(1);
		feedMaterial.GetOverallTPFraction(T[0], P[0], X);

		p1Material.SetOverallProp(ConstCapeString(COBIATEXT("fraction")),
			ConstCapeString(COBIATEXT("mole")), X);
		p1Material.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		p1Material.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), P);

		p2Material.SetOverallProp(ConstCapeString(COBIATEXT("fraction")),
			ConstCapeString(COBIATEXT("mole")), X);
		p2Material.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		p2Material.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), P);


		// Allow all phases to take part in product flash
		p1Material.SetPresentPhases(flashPhaseIDs, flashPhaseStatus);
		p2Material.SetPresentPhases(flashPhaseIDs, flashPhaseStatus);

		// Flash product at specified T & P
		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine p1Equilibrium(p1Material);
		p1Equilibrium.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());

		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine p2Equilibrium(p2Material);
		p2Equilibrium.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean val = true;

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
				val = false;
				break;
			}

			// Get compound list for feed and product
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
		if (val) {

			// Additional checks
			bool sameCompList = (refCompIDs.size() == compIDs.size());
			for (size_t i = 0; (i < compIDs.size()) && (sameCompList); i++) {
				sameCompList = (compIDs[i].compare(refCompIDs[i]) == 0);
			}
			if (!sameCompList) {
				message = COBIATEXT("Connected material streams expose inconsistent compound lists");
				val = false;
			}
			else if (refCompIDs.empty()) {
				message = COBIATEXT("Connected material streams expose zero compound");
				val = false;
			}
		}
		if (val) {
			// Implement ICapeThermoMaterial Interface
			CAPEOPEN_1_2::CapeThermoMaterial p1Material = product1->getMaterial();

			// Phase arrays are identical for product 1 and 2. No need obtain a duplicate
			// CAPEOPEN_1_2::CapeThermoMaterial p2Material = product2->getMaterial();
			
			// Prepare the flash phase list for the product flash
			// This remains constant between validations
			CapeArrayStringImpl stateOfAggregation, keyCompounds;

			CAPEOPEN_1_2::CapeThermoPhases materialPhases(p1Material);
			materialPhases.GetPhaseList(flashPhaseIDs, stateOfAggregation, keyCompounds);
			flashPhaseStatus.resize(flashPhaseIDs.size());
			std::fill(flashPhaseStatus.begin(), flashPhaseStatus.end(), CAPEOPEN_1_2::CAPE_UNKNOWNPHASESTATUS);
		}

		// Update validation status
		validationStatus = (val) ? CAPEOPEN_1_2::CAPE_VALID : CAPEOPEN_1_2::CAPE_INVALID;
		return val;
	}

	// CAPEOPEN_1_2::ICapeUtilities

	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeParameter> getParameters() {
		if (paramCollection != NULL)
		{
			return paramCollection;
		}
		throw cape_open_error(COBIAERR_NotImplemented);
	}
	void putSimulationContext(/*in*/ CAPEOPEN_1_2::CapeSimulationContext context) {
		throw cape_open_error(COBIAERR_NotImplemented);
	}
	void Initialize() {
		// * The PME will order the PMC to get initialized through this method.
		// ** Any initialisation that could fail must be placed here.
		// *** Initialize is guaranteed to be the first method called by the client
		// (except low level methods such as class constructors or initialization persistence methods).
		// **** Initialize has to be called once when the PMC is instantiated in a particular flowsheet.
		// ***** When the initialization fails, before signalling an error,
		// the PMC must free all the resources that were allocated before the failure occurred.
		// When the PME receives this error, it may not use the PMC anymore.
		// The method terminate of the current interface must not either be called.
		// Hence, the PME may only release the PMC through the middleware native mechanisms.
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
		// Specification refrences ECapeUnknown if no GIU is implementd
		throw cape_open_error(COBIAERR_NotImplemented);
	}

	//CAPEOPEN_1_2::ICapePersist

	void Save(/*in*/ CAPEOPEN_1_2::CapePersistWriter writer,/*in*/ CapeBoolean clearDirty) {
		writer.Add(ConstCapeString(COBIATEXT("name")), name);
		writer.Add(ConstCapeString(COBIATEXT("description")), description);
		writer.Add(ConstCapeString(COBIATEXT("Split Ratio")), splitRatio->getValue());
		if (clearDirty) {
			dirty = false;
		}
	}

	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
		splitRatio->putValue(reader.GetReal(ConstCapeString(COBIATEXT("Split Ratio"))));
	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};