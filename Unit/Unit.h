#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "EnergyPort.h"
#include "PortCollection.h"
#include "RealParameter.h"
#include "ParameterCollection.h"

#define projectVersion COBIATEXT("0.3.0")
#define unitDescription COBIATEXT("Simple heater")

#ifdef _DEBUG
#define unitName COBIATEXT("CO PMC Debug")
// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704A8
#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xa8
#else
#define unitName COBIATEXT("CO PMC")
// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704FF
#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xff
#endif

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
	MaterialPortPtr feed1, product1;
	EnergyPortPtr energy;
	RealParameterPtr realParameter1;

	// Collections
	PortCollectionPtr portCollection;
	ParameterCollectionPtr paramCollection;

	// ThermoMaterial Properties
	CAPEOPEN_1_2::CapeThermoMaterial feed1Material, product1Material;
	CapeArrayRealImpl feed1Flow, product1Flow, feed1Enthalpy, product1Enthalpy;
	CapeArrayRealImpl feed1phaseFraction, product1phaseFraction, T, P, X;

	// Flash Conditions
	CapeArrayStringImpl flashCondT, flashCondP, product1PhaseIDs;
	CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> product1FlashPhaseStatus;

	// Persistence and Validation
	CapeBoolean dirty;
	CAPEOPEN_1_2::CapeValidationStatus validationStatus;

public:

	// Returns a description of the current object for error handling
	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Unit: ") + name;
	}

	// Registration info
	static const CapeUUID getObjectUUID() {
		return CapeUUID{unitUUID};
	}

	static void Register(CapePMCRegistrar registrar) {
		registrar.putName(unitName);
		registrar.putDescription(unitDescription);
		registrar.putCapeVersion(COBIATEXT("1.1"));
		registrar.putComponentVersion(projectVersion);
		registrar.putAbout(COBIATEXT("Sample Unit Operation using COBIA."));
		registrar.putVendorURL(COBIATEXT("www.polimi.it"));
		// registrar.putProgId(COBIATEXT("Polimi.Unit"));
		registrar.addCatID(CAPEOPEN::categoryId_UnitOperation);
		registrar.addCatID(CAPEOPEN_1_2::categoryId_Component_1_2);
		//registrar.putCreationFlags(CapePMCRegistationFlag_None);
	}

	Unit()  :
		name(unitName),	description(unitDescription), validationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED),
		feed1(new MaterialPort(COBIATEXT("Feed 1"), CAPEOPEN_1_2::CAPE_INLET, name, validationStatus)),
		product1(new MaterialPort(COBIATEXT("Product 1"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		energy(new EnergyPort(COBIATEXT("Energy"), CAPEOPEN_1_2::CAPE_INLET, name, validationStatus)),
		realParameter1(new RealParameter(name, COBIATEXT("Outlet Temperature"), 373.15, 273.15, 1273.15,
			CAPEOPEN_1_2::CAPE_INPUT, validationStatus, dirty)), dirty(false),
		portCollection(new PortCollection(name)), paramCollection(new ParameterCollection(name)) {

		// Set parameter dimensionality
		realParameter1->putDimensionality(4, 1); // K
		
		// Add ports and parameters to collections
		portCollection->addPort(feed1);
		portCollection->addPort(product1);
		portCollection->addPort(energy);
		paramCollection->addParameter(realParameter1);
		
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
		
		// Check validation status before calculation
		if (validationStatus != CAPEOPEN_1_2::CAPE_VALID) {
			throw cape_open_error(COBIATEXT("Unit is not in a valid state"));
		}

		// Implement ICapeThermoMaterial Interface for feed(s)
		// Implementation for product(s) is done previously in validation method
		feed1Material = feed1->getMaterial();
		
		// Get feed(s) overall props.
		feed1Material.GetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), feed1Flow);
		T.resize(1);
		P.resize(1);
		feed1Material.GetOverallTPFraction(T[0], P[0], X);

		// Material Balance
		product1Flow.resize(feed1Flow.size());
		product1Flow[0] = feed1Flow[0];

		// Energy Balance // TODO: Pressure drop
		T[0] = realParameter1->getValue();

		// Set product(s) overall props.
		product1Material.SetOverallProp(ConstCapeString(COBIATEXT("totalFlow")),
			ConstCapeString(COBIATEXT("mole")), product1Flow);
		product1Material.SetOverallProp(ConstCapeString(COBIATEXT("fraction")),
			ConstCapeString(COBIATEXT("mole")), X);
		product1Material.SetOverallProp(ConstCapeString(COBIATEXT("temperature")),
			ConstCapeEmptyString(), T);
		product1Material.SetOverallProp(ConstCapeString(COBIATEXT("pressure")),
			ConstCapeEmptyString(), P);

		// Allow all phases to take part in product flash
		product1Material.SetPresentPhases(product1PhaseIDs, product1FlashPhaseStatus);

		// Flash product(s) at specified T & P
		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine product1Equilibrium(product1Material);
		product1Equilibrium.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());


		// Enthalpy calculations
		CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;
		CapeArrayStringImpl phaseLabels, props;
		props.resize(1);
		props[0] = L"enthalpy";
		CapeArrayRealImpl enthalpy, phaseFraction, phaseEnthalpy;

		// Feed enthalpy
		CAPEOPEN_1_2::CapeThermoPropertyRoutine feed1Routine(feed1Material);
		feed1Material.GetPresentPhases(phaseLabels, phaseStatus);
		feed1Enthalpy.resize(phaseLabels.size() + 1);
		feed1phaseFraction.resize(phaseLabels.size());
		phaseEnthalpy.resize(phaseLabels.size());
		for (size_t i = 0; i < phaseLabels.size(); i++) {
			feed1Routine.CalcSinglePhaseProp(props, phaseLabels[i]);
			feed1Material.GetSinglePhaseProp(props[0], phaseLabels[i], 
				ConstCapeString(COBIATEXT("mole")), enthalpy);
			feed1Enthalpy[i] = enthalpy[0];
			feed1Material.GetSinglePhaseProp(ConstCapeString(COBIATEXT("phaseFraction")), phaseLabels[i],
				ConstCapeString(COBIATEXT("mole")), phaseFraction);
			feed1phaseFraction[i] = phaseFraction[0];
			phaseEnthalpy[i] = feed1Enthalpy[i] * feed1phaseFraction[i];
		}
		// Calculate feed overall enthalpy
		feed1Enthalpy[phaseLabels.size()] = 0;
		for (size_t i = 0; i < phaseLabels.size(); i++) {
			feed1Enthalpy[phaseLabels.size()] += phaseEnthalpy[i];
		}
		phaseLabels.clear();
		phaseStatus.clear();
		phaseEnthalpy.clear();

		// Product enthalpy
		CAPEOPEN_1_2::CapeThermoPropertyRoutine product1Routine(product1Material);
		product1Material.GetPresentPhases(phaseLabels, phaseStatus);
		product1Enthalpy.resize(phaseLabels.size() + 1);
		product1phaseFraction.resize(phaseLabels.size());
		phaseEnthalpy.resize(phaseLabels.size());
		for (size_t i = 0; i < phaseLabels.size(); i++) {
			product1Routine.CalcSinglePhaseProp(props, phaseLabels[i]);
			product1Material.GetSinglePhaseProp(props[0], phaseLabels[i],
				ConstCapeString(COBIATEXT("mole")), enthalpy);
			product1Enthalpy[i] = enthalpy[0];
			product1Material.GetSinglePhaseProp(ConstCapeString(COBIATEXT("phaseFraction")), phaseLabels[i],
				ConstCapeString(COBIATEXT("mole")), phaseFraction);
			product1phaseFraction[i] = phaseFraction[0];
			phaseEnthalpy[i] = product1Enthalpy[i] * product1phaseFraction[i];
		}
		// Calculate product overall enthalpy
		product1Enthalpy[phaseLabels.size()] = 0;
		for (size_t i = 0; i < phaseLabels.size(); i++) {
			product1Enthalpy[phaseLabels.size()] += phaseEnthalpy[i];
		}

		CapeReal deltaH = product1Enthalpy[product1Enthalpy.size() - 1]
			- feed1Enthalpy[feed1Enthalpy.size() - 1];
		CapeReal work = deltaH * feed1Flow[0];
		
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean val = true;

		// Validate ICapeParameterSpecification & ICapeRealParameterSpecification
		if (realParameter1->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = realParameter1->Validate(message) && realParameter1->Validate(realParameter1->getValue(), message);
		}

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

			if (p.getPortType() == CAPEOPEN_1_2::CAPE_MATERIAL) {
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
			
			// Prepare lists of supported phase label and phase status for product flash
			// This remains constant between validations

			product1Material = product1->getMaterial();
			CAPEOPEN_1_2::CapeThermoPhases product1Phases(product1Material);

			CapeArrayStringImpl stateOfAggregation, keyCompounds;		
			product1Phases.GetPhaseList(product1PhaseIDs, stateOfAggregation, keyCompounds);

			product1FlashPhaseStatus.resize(product1PhaseIDs.size());
			std::fill(product1FlashPhaseStatus.begin(), product1FlashPhaseStatus.end(), CAPEOPEN_1_2::CAPE_UNKNOWNPHASESTATUS);
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
		writer.Add(ConstCapeString(COBIATEXT("Split Ratio")), realParameter1->getValue());
		if (clearDirty) {
			dirty = false;
		}
	}

	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
		realParameter1->putValue(reader.GetReal(ConstCapeString(COBIATEXT("Split Ratio"))));
	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};