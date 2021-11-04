#pragma once
#include <COBIA.h>
#include "Reaction.h"
#include "Helpers.h"

using namespace COBIA;

class ReactionPackage :
	public CapeOpenObject<ReactionPackage>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<ReactionPackage> {

	// Members
	CapeStringImpl packageName, packageDesc;
	CapeArrayStringImpl compIDs;

	std::vector<ReactionPtr> reactions;

public:

	CapeInteger numComponents, numReactions;

	const CapeStringImpl getDescriptionForErrorSource() {
		return packageName + COBIATEXT(" ReactionPackage");
	}

	ReactionPackage(const COBIACHAR* _packageName, const COBIACHAR* _packageDesc) : 
		packageName(_packageName), packageDesc(_packageDesc) {
		numComponents = 0;
		numReactions = 0;
	}

	~ReactionPackage() {
	}

	//CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->packageName;
	}
	
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = this->packageDesc;
	}
	
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	void addComponents(/*in*/ CapeStringImpl _CompIDs) {
		compIDs.emplace_back(_CompIDs);
		numComponents += 1;
	}

	void addReactions(/*in*/ CapeStringImpl _reactionName, /*in*/ CapeStringImpl _reactionDesc,
		/*in*/ CapeReactionType _type, /*in*/ CapeReactionPhase _phase,
		/*in*/ CapeArrayIntegerImpl _stoichiometry, /*in*/ CapeInteger _basisComponentIndex,
		/*in*/ CapeReactionBasis _basis, /*in*/ CapeReal _forwarRateConstant,
		/*in*/ CapeReal _forward_ArrheniusEnergy, /*in*/ CapeReal _heatOfReaction) {
		reactions.emplace_back(
			new Reaction(_reactionName, _reactionDesc, _type, _phase, _stoichiometry,
				_basisComponentIndex, _basis, _forwarRateConstant, _forward_ArrheniusEnergy,
				_heatOfReaction));
		numReactions += 1;
	}

	std::vector<ReactionPtr> getReactions() {
		return reactions;
	}
};

using ReactionPackagePtr = CapeOpenObjectSmartPointer<ReactionPackage>;
