#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "EnergyPort.h"
#include "PortCollection.h"
#include "RealParameter.h"
#include "ParameterCollection.h"
#include "ReactionPackage.h"
#include "Helpers.h"
#include "PFR_Adiabatic_Solver.h"

#define projectVersion COBIATEXT("0.4.0")
#define unitDescription COBIATEXT("PFR with integrated reaction package")

#ifdef _DEBUG
	#ifdef _WIN64
		#define unitName COBIATEXT("CO PMC x64 Debug")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704A9
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xa9
	#else
		#define unitName COBIATEXT("CO PMC x86 Debug")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704A8
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xa8
	#endif
#endif

#ifndef _DEBUG
	#ifdef _WIN64
		#define unitName COBIATEXT("CO PMC x64")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704FF
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xff
	#else
		#define unitName COBIATEXT("CO PMC x86")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A704AA
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x04,0xaa
	#endif
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
	EnergyPortPtr energy;
	RealParameterPtr reactorLength, reactorVolume, conversion;
	RealParameterPtr work, temperatureLow, temperatureHigh;

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
		energy(new EnergyPort(COBIATEXT("Energy"), CAPEOPEN_1_2::CAPE_OUTLET, name, validationStatus)),
		portCollection(new PortCollection(name)), paramCollection(new ParameterCollection(name)),
		reactionPackage(new ReactionPackage(COBIATEXT("reactionPackage1"), COBIATEXT("reactionPackage1Desc"))),
		reactorLength(new RealParameter(name, COBIATEXT("Length"), 6, 0, std::pow(2, 64),
			CAPEOPEN_1_2::CAPE_INPUT, validationStatus, dirty)),
		reactorVolume(new RealParameter(name, COBIATEXT("Volume"), 0.77, 0, std::pow(2, 64),
			CAPEOPEN_1_2::CAPE_INPUT, validationStatus, dirty)),
		conversion(new RealParameter(name, COBIATEXT("Conversion"), 0, 0, 1,
			CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		work(new RealParameter(name, COBIATEXT("work"), 0, -(std::pow(2, 64)), std::pow(2, 64),
			CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		temperatureLow(new RealParameter(name, COBIATEXT("temperatureLow"), 0, 0, 1273.15,
			CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)),
		temperatureHigh(new RealParameter(name, COBIATEXT("temperatureHigh"), 0, 0, 1273.15,
			CAPEOPEN_1_2::CAPE_OUTPUT, validationStatus, dirty)) {

		// Add ports to port collection
		portCollection->addPort(feed1);
		portCollection->addPort(product1);
		portCollection->addPort(energy);

		// Add parameters to parameter collection
		paramCollection->addParameter(reactorLength);
		paramCollection->addParameter(reactorVolume);
		paramCollection->addParameter(conversion);
		
		// Add energy port parameters to energy port (a parameter collection)
		energy->addParameter(work);
		energy->addParameter(temperatureLow);
		energy->addParameter(temperatureHigh);


		// Set parameter dimensionality
		reactorLength->putDimensionality(0, 1);		// m
		reactorVolume->putDimensionality(0, 3);		// m
		reactorVolume->putDimensionality(8, 0);		// -
		work->putDimensionality(0, 2);				// m
		work->putDimensionality(1, 1);				// kg
		work->putDimensionality(2, -3);				// s
		temperatureLow->putDimensionality(4, 1);	// K
		temperatureHigh->putDimensionality(4, 1);	// K

		
		// Set reaction package
		// *** HARDCODING FOR PROOF OF CONCEPT ***
		reactionPackage->addComponents(COBIATEXT("Ethylbenzene"));
		reactionPackage->addComponents(COBIATEXT("Styrene"));
		reactionPackage->addComponents(COBIATEXT("Hydrogen"));

		std::vector<CapeInteger> rxn1Stoichiometry({ -1, 1, 1 });
		reactionPackage->addReactions(COBIATEXT("Rxn1"), COBIATEXT("Ethylene Dehydrogenation"),
			CapeReactionType::CAPE_KINETIC, CapeReactionPhase::CAPE_VAPOR, rxn1Stoichiometry,
			0, CapeReactionBasis::CAPE_MOLE_FRACTION, 4240, 90826.3, 125000);
	}

	// Deconstructor
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

		// CAPE-OPEN unit operations may not have side effects on material objects connected to feeds.
		// Product material object id copied from feed material object allowing to perfrom calculations
		// before setting its properties to override feed peoperties' values
		this->product1->getMaterial().CopyFromMaterial(this->feed1->getMaterial());

		// Initiate Solver
		SolverPtr solver = new Solver(this->product1, this->reactionPackage, this->paramCollection);
		
		// Solving Steps
		solver->getInitialConditions();
		solver->odeSolver();
		solver->setProduct();
		solver->flashProduct(product1PhaseIDs, product1FlashPhaseStatus);

		// Set unit output params
		conversion->putValue(solver->calConversion());

		// Set energy output params - TODO: Set temperatureLow and temperatureHigh
		energy->getCollection()[COBIATEXT("work")].putValue(solver->dH());
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean val = true;

		// Validate ICapeParameterSpecification & ICapeRealParameterSpecification
		// TODO: Iterate over parameter collection instead
		if (reactorLength->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = reactorLength->Validate(message) && reactorLength->Validate(reactorLength->getValue(), message);
		}
		if (reactorVolume->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = reactorVolume->Validate(message) && reactorVolume->Validate(reactorVolume->getValue(), message);
		}
		if (conversion->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = conversion->Validate(message) && conversion->Validate(conversion->getValue(), message);
		}
		if (work->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = work->Validate(message) && work->Validate(work->getValue(), message);
		}
		if (temperatureLow->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = temperatureLow->Validate(message) && temperatureLow->Validate(temperatureLow->getValue(), message);
		}
		if (temperatureHigh->getValStatus() != CAPEOPEN_1_2::CAPE_VALID) {
			val = temperatureHigh->Validate(message) && temperatureHigh->Validate(temperatureHigh->getValue(), message);
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
		writer.Add(ConstCapeString(COBIATEXT("reactorLength")), reactorLength->getValue());
		writer.Add(ConstCapeString(COBIATEXT("reactorVolume")), reactorVolume->getValue());
		writer.Add(ConstCapeString(COBIATEXT("conversion")), conversion->getValue());
		writer.Add(ConstCapeString(COBIATEXT("work")), work->getValue());
		writer.Add(ConstCapeString(COBIATEXT("temperatureLow")), temperatureLow->getValue());
		writer.Add(ConstCapeString(COBIATEXT("temperatureHigh")), temperatureHigh->getValue());

		if (clearDirty) {
			dirty = false;
		}
	}

	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
		// TODO: Iterate over parameter collection
		// and import parameter name using identification interface
		reactorLength->putValue(reader.GetReal(ConstCapeString(COBIATEXT("reactorLength"))));
		reactorVolume->putValue(reader.GetReal(ConstCapeString(COBIATEXT("reactorVolume"))));
		conversion->putValue(reader.GetReal(ConstCapeString(COBIATEXT("conversion"))));
		work->putValue(reader.GetReal(ConstCapeString(COBIATEXT("work"))));
		temperatureLow->putValue(reader.GetReal(ConstCapeString(COBIATEXT("temperatureLow"))));
		temperatureHigh->putValue(reader.GetReal(ConstCapeString(COBIATEXT("temperatureHigh"))));
	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};