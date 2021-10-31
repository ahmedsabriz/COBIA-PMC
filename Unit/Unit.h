#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
// #include "EnergyPort.h"
#include "PortCollection.h"
#include "RealParameter.h"
#include "ParameterCollection.h"
#include "ReactionPackage.h"
#include "Helpers.h"
#include "Solver.h"

#define projectVersion COBIATEXT("0.4.0")
#define unitDescription COBIATEXT("PFR with integrated reaction package")

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

	// Members: Identification
	CapeStringImpl name, description;

	// Members: Persistence and Validation
	CapeBoolean dirty;
	CAPEOPEN_1_2::CapeValidationStatus validationStatus;

	// Members: Ports and Parameters
	MaterialPortPtr feed1, product1;
	// EnergyPortPtr energy;
	// RealParameterPtr work, temperatureLow, temperatureHigh;

	// Members: Collections
	PortCollectionPtr portCollection;
	ParameterCollectionPtr paramCollection;

	// Members: Reaction Packages
	ReactionPackagePtr reactionPackage;

	// Members: Simulation Context
	CAPEOPEN_1_2::CapeSimulationContext simulationContext; // Not Implemented

	// Members: ThermoMaterial Properties
	CAPEOPEN_1_2::CapeThermoMaterial feed1Material; // Migrate
	CAPEOPEN_1_2::CapeThermoMaterial product1Material;
	CapeArrayRealImpl feed1Flow;	// Migrate
	CapeArrayRealImpl product1Flow, feed1Enthalpy, product1Enthalpy;
	CapeArrayRealImpl feed1phaseFraction, product1phaseFraction, T, P, X; //Migrate

	// Members: Flash Conditions
	CapeArrayStringImpl product1PhaseIDs;
	CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> product1FlashPhaseStatus;

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

	// Constructor and members initialisation
	Unit()  :
		name(unitName),	description(unitDescription),
		validationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED), dirty(false),
		feed1(new MaterialPort(COBIATEXT("Feed 1"), CAPEOPEN_1_2::CAPE_INLET, name, validationStatus)),
		product1(new MaterialPort(COBIATEXT("Product 1"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		//energy(new EnergyPort(COBIATEXT("Energy"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		//work(new RealParameter(name, COBIATEXT("work"), 0, -(std::pow(2, 64)), std::pow(2, 64),
		//	CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		//temperatureLow(new RealParameter(name, COBIATEXT("temperatureLow"), 273.15, 73.15, 1273.15,
		//	CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		//temperatureHigh(new RealParameter(name, COBIATEXT("temperatureHigh"), 273.15, 73.15, 1273.15,
		//	CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		portCollection(new PortCollection(name)), paramCollection(new ParameterCollection(name)),
		reactionPackage(new ReactionPackage(COBIATEXT("reactionPackage1"), COBIATEXT("reactionPackage1Desc"))) {

		// Add ports to port collection
		portCollection->addPort(feed1);
		portCollection->addPort(product1);
		//portCollection->addPort(energy);

		// Add parameters to parameter collection
		paramCollection = NULL; // No params
		
		//// Add energy port parameters to energy port (a parameter collection)
		//energy->addParameter(work);
		//energy->addParameter(temperatureLow);
		//energy->addParameter(temperatureHigh);

		// Set parameter dimensionality
		//temperatureLow->putDimensionality(4, 1); // K
		//temperatureHigh->putDimensionality(4, 1); // K
		//work->putDimensionality(0, 2); // m
		//work->putDimensionality(1, 1); // kg
		//work->putDimensionality(2, -3); // s
		
		// Set reaction package
		// *** HARDCODING FOR PROOF OF CONCEPT ***
		reactionPackage->addComponents(COBIATEXT("Ethylbenzene"));
		reactionPackage->addComponents(COBIATEXT("Styrene"));
		reactionPackage->addComponents(COBIATEXT("Hydrogen"));

		
		std::vector<CapeInteger> rxn1Stoichiometry({ -1, 1, 1 });
		reactionPackage->addReactions(COBIATEXT("Rxn1"), COBIATEXT("Ethylene Dehydrogenation"),
			CapeReactionType::CAPE_KINETIC, CapeReactionPhase::CAPE_VAPOR, rxn1Stoichiometry,
			0, CapeReactionBasis::CAPE_MOLE_FRACTION, 4240, -90826.3, 125000);

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

		SolverPtr solver = new Solver(this->feed1->getMaterial(), this->product1->getMaterial(),
			this->reactionPackage->getReactions());
		solver->getInitialConditions();
		solver->odeSolver();
		solver->getInitialConditions();
		solver->setProducts();


		// Allow all phases to take part in product flash
		product1Material.SetPresentPhases(product1PhaseIDs, product1FlashPhaseStatus);

		// Flash product(s) at specified T & P
		CAPEOPEN_1_2::CapeThermoEquilibriumRoutine product1Equilibrium(product1Material);

		// Prepare T & P flash specifications for product flash
		// specification format: CapeArrayRealImpl = { propertyIdentifier, basis, phaseLabel [, compoundIdentifier] }
		// basis is undefined when it is not a dependency of the property (e.g. T, P)
		CapeArrayStringImpl flashCondT(3), flashCondP(3);
		flashCondT[0] = COBIATEXT("temperature");
		flashCondT[2] = COBIATEXT("overall");
		flashCondP[0] = COBIATEXT("pressure");
		flashCondP[2] = COBIATEXT("overall");
		product1Equilibrium.CalcEquilibrium(flashCondT, flashCondP, ConstCapeEmptyString());


		//// Energy Balance
		//energy->getCollection()[L"temperatureLow"].putValue(T[0]);

		//// Enthalpy calculations
		//CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;
		//CapeArrayStringImpl phaseLabels, props;
		//props.resize(1);
		//props[0] = L"enthalpy";
		//CapeArrayRealImpl enthalpy, phaseFraction, phaseEnthalpy;

		//// Feed enthalpy
		//CAPEOPEN_1_2::CapeThermoPropertyRoutine feed1Routine(feed1Material);
		//feed1Material.GetPresentPhases(phaseLabels, phaseStatus);
		//feed1Enthalpy.resize(phaseLabels.size() + 1);
		//feed1phaseFraction.resize(phaseLabels.size());
		//phaseEnthalpy.resize(phaseLabels.size());
		//for (size_t i = 0; i < phaseLabels.size(); i++) {
		//	feed1Routine.CalcSinglePhaseProp(props, phaseLabels[i]);
		//	feed1Material.GetSinglePhaseProp(props[0], phaseLabels[i], 
		//		ConstCapeString(COBIATEXT("mole")), enthalpy);
		//	feed1Enthalpy[i] = enthalpy[0];
		//	feed1Material.GetSinglePhaseProp(ConstCapeString(COBIATEXT("phaseFraction")), phaseLabels[i],
		//		ConstCapeString(COBIATEXT("mole")), phaseFraction);
		//	feed1phaseFraction[i] = phaseFraction[0];
		//	phaseEnthalpy[i] = feed1Enthalpy[i] * feed1phaseFraction[i];
		//}
		//// Calculate feed overall enthalpy
		//feed1Enthalpy[phaseLabels.size()] = 0;
		//for (size_t i = 0; i < phaseLabels.size(); i++) {
		//	feed1Enthalpy[phaseLabels.size()] += phaseEnthalpy[i];
		//}
		//phaseLabels.clear();
		//phaseStatus.clear();
		//phaseEnthalpy.clear();

		//// Product enthalpy
		//CAPEOPEN_1_2::CapeThermoPropertyRoutine product1Routine(product1Material);
		//product1Material.GetPresentPhases(phaseLabels, phaseStatus);
		//product1Enthalpy.resize(phaseLabels.size() + 1);
		//product1phaseFraction.resize(phaseLabels.size());
		//phaseEnthalpy.resize(phaseLabels.size());
		//for (size_t i = 0; i < phaseLabels.size(); i++) {
		//	product1Routine.CalcSinglePhaseProp(props, phaseLabels[i]);
		//	product1Material.GetSinglePhaseProp(props[0], phaseLabels[i],
		//		ConstCapeString(COBIATEXT("mole")), enthalpy);
		//	product1Enthalpy[i] = enthalpy[0];
		//	product1Material.GetSinglePhaseProp(ConstCapeString(COBIATEXT("phaseFraction")), phaseLabels[i],
		//		ConstCapeString(COBIATEXT("mole")), phaseFraction);
		//	product1phaseFraction[i] = phaseFraction[0];
		//	phaseEnthalpy[i] = product1Enthalpy[i] * product1phaseFraction[i];
		//}
		//// Calculate product overall enthalpy
		//product1Enthalpy[phaseLabels.size()] = 0;
		//for (size_t i = 0; i < phaseLabels.size(); i++) {
		//	product1Enthalpy[phaseLabels.size()] += phaseEnthalpy[i];
		//}

		//// energy->getCollection()[L"work"].putValue(-(realParameter2->getValue()));
		//// Negative sign becaue energy is set as output in a heater
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean val = true;

		// Validate ICapeParameterSpecification & ICapeRealParameterSpecification
		// TODO: Iterate over parameter collection instead

		//if (work->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
		//	val = work->Validate(message) && work->Validate(work->getValue(), message);
		//}

		//if (temperatureLow->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
		//	val = temperatureLow->Validate(message) && temperatureLow->Validate(temperatureLow->getValue(), message);
		//}

		//if (temperatureHigh->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
		//	val = temperatureHigh->Validate(message) && temperatureHigh->Validate(temperatureHigh->getValue(), message);
		//}

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
		// Storing a reference in case any COSE implementation is needed in future
		this->simulationContext = context;
		// For now an error is thrown
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
		CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> collection(portCollection);
		for (CAPEOPEN_1_2::CapeUnitPort p : collection) {
			p.Disconnect();
		}
		// In case a reference to the simulation context is stored, it too must be released at Terminate
		simulationContext.clear();

	}
	CAPEOPEN_1_2::CapeEditResult Edit(CapeWindowId parent) {
		// TODO
		// For now, an error is thrown
		throw cape_open_error(COBIAERR_NotImplemented);
	}

	//CAPEOPEN_1_2::ICapePersist

	void Save(/*in*/ CAPEOPEN_1_2::CapePersistWriter writer,/*in*/ CapeBoolean clearDirty) {
		writer.Add(ConstCapeString(COBIATEXT("name")), name);
		writer.Add(ConstCapeString(COBIATEXT("description")), description);
		// TODO: Iterate over parameter collection
		// and import parameter name using identification interface
		//writer.Add(ConstCapeString(COBIATEXT("work")), work->getValue());
		//writer.Add(ConstCapeString(COBIATEXT("temperatureLow")), temperatureLow->getValue());
		//writer.Add(ConstCapeString(COBIATEXT("temperatureHigh")), temperatureHigh->getValue());


		if (clearDirty) {
			dirty = false;
		}
	}

	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
		// TODO: Iterate over parameter collection
		// and import parameter name using identification interface
		//work->putValue(reader.GetReal(ConstCapeString(COBIATEXT("work"))));
		//temperatureLow->putValue(reader.GetReal(ConstCapeString(COBIATEXT("temperatureLow"))));
		//temperatureHigh->putValue(reader.GetReal(ConstCapeString(COBIATEXT("temperatureHigh"))));

	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};