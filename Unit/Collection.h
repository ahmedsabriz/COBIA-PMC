#pragma once
#include <COBIA.h>

using namespace COBIA;

template <typename CollectionItemInterface, typename CollectionItem> class Collection :
	public CapeOpenObject<Collection<CollectionItemInterface,CollectionItem>>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Collection<CollectionItemInterface, CollectionItem>>,
	public CAPEOPEN_1_2::CapeCollectionAdapter<CollectionItemInterface,Collection<CollectionItemInterface, CollectionItem>> {

	/*
	Implementation of Collection Common Interface
	"The ICapeCollection interface provides a means of collecting together
	lists of CAPE-OPEN items/entities (eg. parameters, ports, …)."
	*/

	// Members
	CapeStringImpl& unitName;
	std::vector<CollectionItem> items;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Collection of ") + unitName;
	}

	Collection(CapeStringImpl& _unitName) : unitName(_unitName) {
	}

	~Collection() {
	}

	// Method to adding ports to the collection
	void addItem(CollectionItem item) {
		items.emplace_back(item);
	}

	// Method to removing ports from the collection - Not Implemented

	// Lookup by index and return a CollectionItem
	CollectionItem getItemImpl(/*in*/ CapeInteger index) {
		if ((index < 0) || (index >= (CapeInteger)items.size())) {
			throw cape_open_error(COBIAERR_NoSuchItem);
		}
		return items[index];
	}
	
	std::vector<CollectionItem> iterateOverItems() {
		return items;
	}

	// CAPEOPEN_1_2::ICapeIdentification
	void getComponentName(/*out*/ CapeString name) {
		name = COBIATEXT("Collection");
	}
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Collection");
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	

	// CAPEOPEN_1_2::ICapeCollection<CAPEOPEN_1_2::ICapeUnitPort>

	// Lookup by index
	CollectionItemInterface Item(/*in*/ CapeInteger index) {
		if ((index < 0) || (index >= (CapeInteger)items.size())) {
			throw cape_open_error(COBIAERR_NoSuchItem);
		}
		return items[index];
	}

	// Lookup by name
	CollectionItemInterface Item(/*in*/ CapeString name) {
		CapeStringImpl itemName;
		for (CollectionItem& item : items) {
			CAPEOPEN_1_2::CapeIdentification identification(item);
			identification.getComponentName(itemName);
			if (itemName == name) {
				return item;
			}
		}
		// If not found (no return)
		throw cape_open_error(COBIAERR_NoSuchItem);
	}
	CapeInteger getCount() {
		return (CapeInteger)items.size();
	}
	
};

template <typename CollectionItemInterface, typename CollectionItem> using CollectionPtr = CapeOpenObjectSmartPointer<Collection<CollectionItemInterface,CollectionItem>>;
