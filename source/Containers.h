#pragma once

#include <cstdlib>
#include <set>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/composite_key.hpp> 
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "Host.h"
#include "Event.h"

#ifdef _DEBUG
#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE
#endif

struct age_tag{};
struct household_tag{};
struct eligible_age_household_tag{};
struct eligible_household_tag{};
struct age_household_tag{};
struct age_neighborhood_tag{};
struct carriage_tag{};

typedef boost::multi_index::multi_index_container<
	boost::shared_ptr<Host>,
	boost::multi_index::indexed_by< 
	    boost::multi_index::hashed_unique< // 0 - ID index
		    boost::multi_index::const_mem_fun<Host, int, &Host::getID>>, 
		boost::multi_index::ordered_non_unique< // 1 - Age index
		    boost::multi_index::tag<age_tag>,
			boost::multi_index::const_mem_fun<Host, int, &Host::getAgeInY>>, 
		boost::multi_index::hashed_non_unique < // 2 - Household index
		    boost::multi_index::tag<household_tag>,
			boost::multi_index::const_mem_fun < Host, int, &Host::getHousehold >> ,
	    boost::multi_index::ordered_non_unique< // 3 - Eligible by age & household
		    boost::multi_index::tag<eligible_age_household_tag>,
			boost::multi_index::composite_key <
			    Host,
				boost::multi_index::const_mem_fun<Host, int, &Host::getAgeInY>,
				boost::multi_index::const_mem_fun<Host, bool, &Host::isEligible>,
				boost::multi_index::const_mem_fun<Host, int, &Host::getHousehold>>>,
		boost::multi_index::ordered_non_unique < // 4 - Eligible by household (all single adults)
		    boost::multi_index::tag<eligible_household_tag>,
			boost::multi_index::composite_key <
		        Host,
				boost::multi_index::const_mem_fun<Host, bool, &Host::isEligible>,
				boost::multi_index::const_mem_fun<Host, int, &Host::getHousehold >> >,
		boost::multi_index::ordered_non_unique < // 5 - Neighborhood & age
		    boost::multi_index::tag<age_neighborhood_tag>,
			boost::multi_index::composite_key <
		        Host,
				boost::multi_index::const_mem_fun<Host, int, &Host::getNeighborhood>,
				boost::multi_index::const_mem_fun<Host, int, &Host::getAgeInY >> >,
		boost::multi_index::ordered_non_unique < // 6 - Household & age
		    boost::multi_index::tag<age_household_tag>,
			boost::multi_index::composite_key <
	            Host,
				boost::multi_index::const_mem_fun<Host, int, &Host::getHousehold>,
				boost::multi_index::const_mem_fun<Host, int, &Host::getAgeInY >> >,
		boost::multi_index::ordered_non_unique <  // 7 - Carriage status
	        boost::multi_index::tag<carriage_tag>,
			boost::multi_index::const_mem_fun < Host, bool, &Host::isInfected >> >> HostContainer;

typedef HostContainer::nth_index<0>::type HostsByID;
typedef HostContainer::nth_index<1>::type HostsByAge;
typedef HostContainer::nth_index<2>::type HostsByHousehold;
typedef HostContainer::nth_index<3>::type HostsEligibleByAgeHousehold;
typedef HostContainer::nth_index<4>::type HostsEligibleByHousehold;
typedef HostContainer::nth_index<5>::type HostsByAgeNeighborhood;
typedef HostContainer::nth_index<6>::type HostsByAgeHousehold;
typedef HostContainer::nth_index<7>::type HostsByCarriage;

typedef std::set<int> HouseholdContainer;
