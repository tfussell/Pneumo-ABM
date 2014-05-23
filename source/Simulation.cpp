#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <boost/random.hpp>

#include "Simulation.h"
#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Containers.h"

#define HFLU_INDEX (INIT_NUM_STYPES-1)

Simulation::Simulation(int trt, int sid, SimPars * spPtr) : numInfecteds() {
    simID = sid;
    treatment = trt;
    simParsPtr = spPtr;
    t = 0;
    demOutputStrobe = 0;
    eventCtr = 0;
    idCtr = 1;
    hholdCtr = 1;
    demComplete = 0.0;
    rng.seed((unsigned int)sid);

    // Set up multinomial to determine # of households in each size category
    // and # of people in each age class
    std::array<int, INIT_NUM_AGE_CATS> init_age_dist = {0};
    rmultinom(simParsPtr->get_demPMF_row(INIT_AGE_INDEX), N0, INIT_NUM_AGE_CATS, init_age_dist.data(), rng);

    // Create hosts
    double thisAge;
    double DOB;
    int randNeighborhood;
    for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
        while(init_age_dist[a] > 0) {
            thisAge = 365.0*((double)a + r01(rng));
            DOB = t - thisAge;
            randNeighborhood = (int)(floor((double)NUM_NEIGHBORHOODS*r01(rng)));
            Host * newHostPtr = new Host(t, DOB, idCtr, hholdCtr, randNeighborhood, currentEvents, simParsPtr, rng);
            allHosts.insert(boost::shared_ptr<Host>(newHostPtr));
            allHouseholds.insert(hholdCtr);
            idCtr++;
            hholdCtr++;
            init_age_dist[a]--;
        } // end while (= this age category exhausted)
    } // all age categories done

    // Assign to households
    HostsByAge::iterator it = allHosts.get<age_tag>().end();
    while(it != allHosts.get<age_tag>().begin()) {
        it--;
        if((*it)->isEligible()) { // If host is an adult and not yet paired
            double pairs = r01(rng);
            if(pairs < PROB_PAIR / 2.0) {
                int thisID = (*it)->getID();
                pairHost(thisID);
            }
        }
        else if((*it)->isAdult() == false) { // If host is a kid, pick hosts randomly until find adult and join household
            bool a = false;
            int rID = 0;
            HostsByID::iterator ita;
            while(a == false) {
                rID = static_cast<int>(ceil((double)N0*r01(rng)));
                ita = allHosts.find(rID);
                a = (*ita)->isAdult();
            }
            int newHH = (*ita)->getHousehold();
            int newNH = (*ita)->getNeighborhood();
            HostsByID::iterator itt = allHosts.find((*it)->getID());
            allHouseholds.erase((*itt)->getHousehold());
            allHosts.modify(itt, [=](boost::shared_ptr<Host> h) { h->setHousehold(newHH); });
            allHosts.modify(itt, [=](boost::shared_ptr<Host> h) { h->setNeighborhood(newNH); });
        }
    } // end while
    initOutput();
}

Simulation::~Simulation() {
    closeOutput();
}


// Stream functions
void Simulation::closeOutput() {
    ageDistStream.close();
    hhDistStream.close();
    demTimesStream.close();

    coinfectionHistStream.close();
    coinfectionHFHistStream.close();
    epidTimesStream.close();
    totCarriageStream.close();
}

void Simulation::initOutput() {
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        ageDistFile = makeBigName("age_dist_neighborhood", n);
        ageDistStream.open("../../outputs/" + ageDistFile, std::ios::out);
        ageDistStream.close();
    }
    hhDistFile = makeName("hh_dist");
    demTimesFile = makeName("dem_times");
    epidTimesFile = makeName("epid_times");
    coinfectionHistFile = makeName("coinfection_dist_pneumo");
    coinfectionHFHistFile = makeName("coinfection_dist_hflu-pneumo");
    totCarriageFile = makeName("totCarriage");

    hhDistStream.open("../../outputs/" + hhDistFile, std::ios::out);
    demTimesStream.open("../../outputs/" + demTimesFile, std::ios::out);
    epidTimesStream.open("../../outputs/" + epidTimesFile, std::ios::out);
    coinfectionHistStream.open("../../outputs/" + coinfectionHistFile, std::ios::out);
    coinfectionHFHistStream.open("../../outputs/" + coinfectionHFHistFile, std::ios::out);
    totCarriageStream.open("../../outputs/" + totCarriageFile, std::ios::out);

    for(int s = 0; s < INIT_NUM_STYPES; s++) {
        for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
            infectionStream.open("../../outputs/" + makeBiggerName("infections", s, "neighborhood", n), std::ios::out);
            infectionStream.close();
            infectedStream.open("../../outputs/" + makeBiggerName("infecteds", s, "neighborhood", n), std::ios::out);
            infectedStream.close();
        }
    }
}

void Simulation::writeDemOutput() {

    // I. Output age distribution to file
    HostsByAge::iterator it = allHosts.get<age_tag>().end();
    it--;
    int maxAge = INIT_NUM_AGE_CATS - 1;

    // For 0 to max age, count ages 
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        ageDistFile = makeBigName("age_dist_neighborhood", n);
        ageDistStream.open("../../outputs/" + ageDistFile, std::ios::app);
        for(int a = 0; a < maxAge + 1; a++) {
            ageDistStream << allHosts.get<age_neighborhood_tag>().count(boost::make_tuple(n, a)) << "\t";
        }
        ageDistStream << std::endl;
        ageDistStream.close();
    }

    // Update dem times
    demTimesStream << demOutputStrobe << std::endl;

    // II. Output household distribution to file
    HostsByHousehold& sorted_index2 = allHosts.get<household_tag>();
    int hhSize = 0;
    int maxSize = 1;
    int households[HHOLD_SIZE_BUFFER]; // holds counts for sizes 1...HHOLD_SIZE_BUFFER+1
    for(int s = 0; s < HHOLD_SIZE_BUFFER; s++) {
        households[s] = 0;
    }
    for(auto household : allHouseholds) {
        hhSize = static_cast<int>(sorted_index2.count(household));
        if(hhSize > maxSize) {
            maxSize = hhSize;
        }
        if(hhSize > floor(0.8*(double)HHOLD_SIZE_BUFFER)) {
            std::cout << "Set MAX_HHOLD_SIZE to be larger. Encountered household with " << hhSize << " members." << std::endl;
            assert(hhSize < HHOLD_SIZE_BUFFER);
        }
        households[hhSize - 1]++;
    }
    for(int s = 0; s < HHOLD_SIZE_BUFFER; s++) {
        hhDistStream << households[s] << "\t";
    }
    hhDistStream << std::endl;
}

void Simulation::writeTheta() {
    thetaFile = makeName("theta");
    thetaStream.open("../../outputs/" + thetaFile, std::ios::out);
    for(auto host : allHosts.get<age_tag>()) { // For each host...
        if(host->getAgeInY() < 6) {
            thetaStream << t - host->getDOB() << "\t" << host->getSummedTheta() << std::endl;
        }
    }
    thetaStream.close();
}

void Simulation::writeEpidOutput() {
    // I. Output infections by age to file
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) { // Fix when adapt to more than one neighborhood
        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            infectionStream.open("../../outputs/" + makeBiggerName("infections", s, "neighborhood", n), std::ios::app);
            for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                infectionStream << numInfecteds[a][s][n] << "\t";
            }
            infectionStream << std::endl;
            infectionStream.close();
        }
    }

    // II. Write time file
    epidTimesStream << epidOutputStrobe << std::endl;

    // III. Write infecteds by age to file for each serotype
    int actualInfecteds[INIT_NUM_AGE_CATS][INIT_NUM_STYPES][NUM_NEIGHBORHOODS] = {0};

    // IV. Write coinfections & total carriage to files
    int numCarryingPneumo = 0;
    int coinf[HFLU_INDEX];
    int coinfHflu[HFLU_INDEX];
    for(int c = 0; c < HFLU_INDEX; c++) {
        coinf[c] = 0;
        coinfHflu[c] = 0;
    }
    int thisC;
    auto pit = allHosts.get<carriage_tag>().equal_range(true);
    for(auto fit = pit.first; fit != pit.second; fit++) {
        if((*fit)->isInfectedPneumo()) {
            numCarryingPneumo++;
        }
        if((*fit)->isInfectedPneumo() && (*fit)->getAgeInY() < COCOL_AGE_LIMIT) {
            bool hasHflu = (*fit)->isInfectedHflu();
            thisC = -1;
            for(int z = 0; z < HFLU_INDEX; z++) {
                if((*fit)->isInfectedZ(z)) {
                    thisC++;
                    if(hasHflu) {
                        coinfHflu[z]++;
                    }
                }
            }
            coinf[thisC]++;
        }

        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            if((*fit)->isInfectedZ(s)) {
                actualInfecteds[(*fit)->getAgeInY()][s][(*fit)->getNeighborhood()]++;
            }
        } // end for each serotype
    } // end for all infected hosts
    for(int c = 0; c < HFLU_INDEX; c++) {
        coinfectionHistStream << coinf[c] << "\t";
        coinfectionHFHistStream << coinfHflu[c] << "\t";
    }
    coinfectionHistStream << std::endl;
    coinfectionHFHistStream << std::endl;
    totCarriageStream << numCarryingPneumo << std::endl;

    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            infectedStream.open("../../outputs/" + makeBiggerName("infecteds", s, "neighborhood", n), std::ios::app);
            for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                infectedStream << actualInfecteds[a][s][n] << "\t";
            }
            infectedStream << std::endl;
            infectedStream.close();
        }
    }
}


// Member function definitions

void Simulation::runDemSim() {
    double percentDone = 0.0;
    std::cout << "  Starting demographic component with seed=" << simID << "." << std::endl;
    writeDemOutput();
    demOutputStrobe += STROBE_DEM;
    EventQueue::iterator eventIter = currentEvents.begin();
    while((*eventIter).time < DEM_SIM_LENGTH  && currentEvents.size() != 0 && allHosts.size() > 0)   {
        Event thisEvent = *eventIter;
        t = thisEvent.time;
        while(demOutputStrobe < t) {
            writeDemOutput();
            demOutputStrobe += STROBE_DEM;
        }
        while(percentDone / 100.0 <= t / DEM_SIM_LENGTH) {
            std::cout << "\t" << percentDone << "% of this component complete." << std::endl;
            percentDone += PROGRESS_INTERVAL;
        }
        executeEvent(thisEvent);
        eventCtr++;
        currentEvents.erase(eventIter);
        eventIter = currentEvents.begin();
    } // while time < DEM_SIM_LENGTH

    // Set time to end of simulation and reset strobe to calibrate with epid strobing
    t = DEM_SIM_LENGTH;
    demComplete = t;
    std::cout << "\t100% of this component complete." << std::endl;
}

std::array<double, INIT_NUM_STYPES> Simulation::runTestEpidSim() {
    if(allHosts.size() == 0) {
        std::cerr << "No hosts remaining for epidemiological simulation. Cancelling." << std::endl;
        assert(false);
    }
    std::cout << "  Entering test simulation at t=" << demComplete << "." << std::endl;
    double percentDone = 0.0;

    // Initialize host population with infections
    demOutputStrobe = t;
    epidOutputStrobe = t;
    seedInfections();
    EventQueue::iterator eventIter = currentEvents.begin();
    double nextTimeStep = t + EPID_DELTA_T;

    // holds prevalences at strobing periods
    std::array<double, INIT_NUM_STYPES> serotype_prevalence_rates = {0};

    double prevYear = DEM_SIM_LENGTH + TEST_EPID_SIM_LENGTH - (NUM_TEST_SAMPLES * 365.0); // first simulation time to start sampling
    int prevSamples = 0;

    while(t < TEST_EPID_SIM_LENGTH + demComplete)
    {
        // Calculate new infections for every host and add events to stack
        calcSI();
        eventIter = currentEvents.begin();
        while((*eventIter).time < nextTimeStep) {
            while(demOutputStrobe < t) {
                writeDemOutput();
                demOutputStrobe += STROBE_DEM;
            }
            while(epidOutputStrobe < t) {
                writeEpidOutput();
                epidOutputStrobe += STROBE_EPID;
            }
            if(prevYear < t) {
                auto period_prevalence_rates = calculateSerotypePrevalenceRates();
                for(int i = 0; i < INIT_NUM_STYPES; i++)
                {
                    serotype_prevalence_rates[i] += period_prevalence_rates[i] / NUM_TEST_SAMPLES;
                }
                //std::cout << "\tOutputting prevalence sample #" << prevSamples + 1 << "; prevalence of pneumo under 5 is " << std::accumulate(prevalences[prevSamples].begin(), prevalences[prevSamples].end(), 0.0) << std::endl;
                prevSamples++;
                prevYear += 365.0;
            }
            while(percentDone / 100.0 <= (t - demComplete) / TEST_EPID_SIM_LENGTH) {
                std::cout << "\t" << percentDone << "% of this test component complete." << std::endl;
                percentDone += PROGRESS_INTERVAL;
            }

            Event thisEvent = *eventIter;
            t = thisEvent.time;
            executeEvent(thisEvent);
            eventCtr++;
            currentEvents.erase(eventIter);
            eventIter = currentEvents.begin();
        }
        t = nextTimeStep;
        nextTimeStep += EPID_DELTA_T;
    }

    std::cout << "\t100% of this test component complete." << std::endl;

    return serotype_prevalence_rates;
}

void Simulation::runEpidSim() {
    if(allHosts.size() == 0) {
        std::cerr << "No hosts remaining for epidemiological simulation. Cancelling." << std::endl;
        assert(false);
    }
    std::cout << "  Entering epidemiological simulation at t=" << demComplete << "." << std::endl;
    double percentDone = 0.0;

    // Initialize host population with infections
    demOutputStrobe = t;
    epidOutputStrobe = t;
    seedInfections();
    EventQueue::iterator eventIter = currentEvents.begin();
    double nextTimeStep = t + EPID_DELTA_T;

    while(t < EPID_SIM_LENGTH + demComplete)
    {
        // Calculate new infections for every host and add events to stack
        calcSI();
        eventIter = currentEvents.begin();
        while((*eventIter).time < nextTimeStep) {
            while(demOutputStrobe < t) {
                writeDemOutput();
                demOutputStrobe += STROBE_DEM;
            }
            while(epidOutputStrobe < t) {
                writeEpidOutput();
                epidOutputStrobe += STROBE_EPID;
            }
            while(percentDone / 100.0 < (t - demComplete) / EPID_SIM_LENGTH) {
                std::cout << "\t" << percentDone << "% of this component complete." << std::endl;
                percentDone += PROGRESS_INTERVAL;
            }
            Event thisEvent = *eventIter;
            t = thisEvent.time;
            executeEvent(thisEvent);
            eventCtr++;
            currentEvents.erase(eventIter);
            eventIter = currentEvents.begin();
        }
        t = nextTimeStep;
        nextTimeStep += EPID_DELTA_T;
    }
    std::cout << "\t100% of this component complete." << std::endl;
    writeTheta();
    std::cout << "  At end of simulation, " << allHosts.size() << " hosts and " << allHouseholds.size() << " households; cumulatively, " << idCtr - 1 << " hosts and " << hholdCtr - 1 << " households. " << eventCtr << " total events." << std::endl;
}


// PRIVATE SIMULATION FUNCTIONS
void Simulation::executeEvent(Event & te) {
    switch(te.type) {
    case Event::Type::Death:
        killHost(te.hostID);
        break;
    case Event::Type::Fledge:
        fledgeHost(te.hostID);
        break;
    case Event::Type::Pair:
        pairHost(te.hostID);
        break;
    case Event::Type::Birth:
        birthHost(te.hostID);
        break;
    case Event::Type::Birthday:
        ageHost(te.hostID);
        break;
    case Event::Type::Infection:
        infectHost(te.hostID, te.s);
        break;
    case Event::Type::Recovery:
        recoverHost(te.hostID, te.s);
        break;
    case Event::Type::Vaccination:
        vaccinateHost(te.hostID);
        break;
    default:
        throw std::runtime_error("invalid event");
    }
}

void Simulation::seedInfections() {
    // Currently assuming: 
    // .... only single colonizations occur for each serotype
    // .... not considering co-colonizations when calculating duration of infection at this point
    for(auto host : allHosts.get<age_tag>())
    {
        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            double thisRNG = r01(rng);
            if(thisRNG < simParsPtr->get_serotypePar_ij(INIT_INFECTEDS_INDEX, s)) {
                int id = host->getID();
                infectHost(id, s);
            }
        }
    }
}

void Simulation::killHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    // If only member of household, remove household from allHouseholds
    int hid = (*it)->getHousehold();
    auto hsize = allHosts.get<household_tag>().count(hid);
    if(hsize == 1) {
        allHouseholds.erase(hid);
    }

    // Update numInfecteds
    if(demComplete != 0) {
        if((*it)->isInfected()) {
            int age = (*it)->getAgeInY();
            int nhood = (*it)->getNeighborhood();
            for(int s = 0; s < INIT_NUM_STYPES; s++) {
                numInfecteds[age][s][nhood] -= (*it)->isInfectedZ(s);
            }
        }
    }
    allHosts.erase(it);
}

void Simulation::ageHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    if(demComplete != 0) {
        if((*it)->isInfected()) {
            int age = (*it)->getAgeInY();
            int nhood = (*it)->getNeighborhood();
            for(int s = 0; s < INIT_NUM_STYPES; s++) {
                numInfecteds[age][s][nhood] -= (*it)->isInfectedZ(s);
                numInfecteds[age + 1][s][nhood] += (*it)->isInfectedZ(s);
            } // end for each serotype
        } // end for if host infected
    } // end for if in epid part of simulation
    allHosts.modify(it, [](boost::shared_ptr<Host> h) { h->incrementAge(); });
}

void Simulation::pairHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    if((*it)->isPaired() == false) {

        // Identify absorbing household (= initiating host's household )
        int hhold1 = (*it)->getHousehold();

        // Find eligible partner 
        int thisAge = (*it)->getAgeInY();
        int hhold2 = hhold1; // ensure do not pair w/in own household
        int count = 0;
        int partnerAge = 0;
        int attempts = 0;
        int rInd = 0;
        int id2 = 0;
        int famCount = 1;
        bool foundPartner = false;
        while(count == 0 && attempts < DATE_BUFFER) { // First attempt to find partner in preferred age range
            partnerAge = calcPartnerAge(thisAge);
            count = static_cast<int>(allHosts.get<eligible_age_household_tag>().count(boost::make_tuple(partnerAge, true)));
            famCount = static_cast<int>(allHosts.get<eligible_age_household_tag>().count(boost::make_tuple(partnerAge, true, hhold1)));
            if(count == 0 || famCount == count) {
                attempts++;
            }
            else {
                while(hhold1 == hhold2 && attempts < DATE_BUFFER) {
                    rInd = static_cast<int>(ceil((double)count * r01(rng)));
                    auto pit = allHosts.get<eligible_age_household_tag>().equal_range(boost::make_tuple(partnerAge, true));
                    auto partnerIt = pit.first;
                    for(int i = 1; i < rInd; i++) {
                        partnerIt++;
                    }
                    hhold2 = (*partnerIt)->getHousehold();
                    if(hhold1 != hhold2) { // Found partner; will update marital status
                        id2 = (*partnerIt)->getID();
                        foundPartner = true;
                    }
                    attempts++;
                } // end while hhold1==hhold2 && attempts < MAX_DATES
            } // end if (count == 0 ){} else{}
        } // end while count==0 && attempts < MAX_DATES

        // If no one found, search outside original range
        if(foundPartner == false) {
            count = static_cast<int>(allHosts.get<eligible_household_tag>().count(boost::make_tuple(true)));
            int famCount = static_cast<int>(allHosts.get<eligible_household_tag>().count(boost::make_tuple(true, hhold1)));
            attempts = 0;
            hhold2 = hhold1;

            // Check if there are any not in household
            if(count - famCount > 0) {
                while(hhold1 == hhold2) {
                    rInd = static_cast<int>(ceil((double)count * r01(rng)));
                    auto pit = allHosts.get<eligible_household_tag>().equal_range(boost::make_tuple(true));
                    auto partnerIt = pit.first;
                    for(int i = 1; i < rInd; i++) {
                        partnerIt++;
                    }
                    hhold2 = (*partnerIt)->getHousehold();
                    if(hhold1 != hhold2) { // Found partner; will update marital status
                        id2 = (*partnerIt)->getID();
                        foundPartner = true;
                    }
                } // end while (hhold1==hhold2)
            } // end if count > 0
        } // end if count == 0

        if(foundPartner == true) {
            partner2Hosts(id, id2);
        } // end if partner found (count > 0)
    } // end if host eligible
}

void Simulation::partner2Hosts(int hid1, int hid2) {
    // Update both partners' partner IDs and get household IDs
    HostsByID::iterator it1 = allHosts.find(hid1); // Original host
    int hhold1 = (*it1)->getHousehold();
    int nhood1 = (*it1)->getNeighborhood();
    allHosts.modify(it1, [=](boost::shared_ptr<Host> h) { h->setPartner(hid2); }); // Partner host
    HostsByID::iterator it2 = allHosts.find(hid2); // Partner host
    int hhold2 = (*it2)->getHousehold();
    allHosts.modify(it2, [=](boost::shared_ptr<Host> h) { h->setPartner(hid1); });

    // Check if initiating host still lives at home. If so, fledge and give new household & neighborhood.
    if((*it1)->hasFledged() == false) {
        allHosts.modify(it1, [](boost::shared_ptr<Host> h) { h->setFledge(true); });
        int famSize = static_cast<int>(allHosts.get<household_tag>().count(hhold1));
        if(famSize != 1) {
            allHosts.modify(it1, [=](boost::shared_ptr<Host> h) { h->setHousehold(hholdCtr); });
            int randNeighborhood = (int)(floor((double)NUM_NEIGHBORHOODS*r01(rng)));
            allHosts.modify(it1, [=](boost::shared_ptr<Host> h) { h->setNeighborhood(randNeighborhood); });
            if((*it1)->isInfected()) {
                for(int s = 0; s < INIT_NUM_STYPES; s++) {
                    int infections = (*it1)->isInfectedZ(s);
                    if(infections > 0) {
                        numInfecteds[(*it1)->getAgeInY()][s][nhood1] -= infections;
                        numInfecteds[(*it1)->getAgeInY()][s][randNeighborhood] += infections;
                    }
                }
            }
            allHouseholds.insert(hholdCtr);
            hholdCtr++;
            hhold1 = (*it1)->getHousehold();
            nhood1 = (*it1)->getNeighborhood();
        }
    }

    // Update household of partner and partner's kids, if any
    int numFamily = static_cast<int>(allHosts.get<household_tag>().count(hhold2));
    bool fledged = (*it2)->hasFledged();

    // Remove household2 from master list of households, IF host lives alone or has fledged
    if(numFamily == 1) {
        allHouseholds.erase(hhold2);
    }

    if(numFamily > 1 && fledged == true) { // New partner living independently with own family
        // Get IDs of all people whose households need updating
        std::vector<int> familyIDs(numFamily, 0);
        int indID = 0;
        int f = 0;
        auto pit = allHosts.get<household_tag>().equal_range(hhold2);
        for(auto fit = pit.first; fit != pit.second; fit++) {
            indID = (*fit)->getID();
            familyIDs[f] = indID;
            f++;
        }
        // Now update households and neighborhoods
        HostsByID::iterator it;
        for(int f = 0; f < numFamily; f++) {
            indID = familyIDs[f];
            it = allHosts.find(indID);
            int thisAge = (*it)->getAgeInY();
            int origNhood = (*it)->getNeighborhood();
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setHousehold(hhold1); });
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setNeighborhood(nhood1); });
            if((*it)->isInfected()) {
                for(int s = 0; s < INIT_NUM_STYPES; s++) {
                    int infections = (*it)->isInfectedZ(s);
                    numInfecteds[thisAge][s][origNhood] -= infections;
                    numInfecteds[thisAge][s][nhood1] += infections;
                }
            }
        }
        allHouseholds.erase(hhold2);
    }
    else {  // New partner still living in family of origin, so just update partner
        int origNhood = (*it2)->getNeighborhood();
        int age = (*it2)->getAgeInY();
        allHosts.modify(it2, [=](boost::shared_ptr<Host> h) { h->setHousehold(hhold1); });
        allHosts.modify(it2, [=](boost::shared_ptr<Host> h) { h->setNeighborhood(nhood1); });
        if((*it2)->isInfected()) {
            for(int s = 0; s < INIT_NUM_STYPES; s++) {
                int infections = (*it2)->isInfectedZ(s);
                numInfecteds[age][s][origNhood] -= infections;
                numInfecteds[age][s][nhood1] += infections;
            }
        }
        allHosts.modify(it2, [=](boost::shared_ptr<Host> h) { h->setFledge(true); });
    }
}

int Simulation::calcPartnerAge(int a) {
    double partnerAge = (double)MATURITY_AGE - 1.0;
    double rand_offset = 0.0;
    while(partnerAge < MATURITY_AGE) {
        rand_offset = r01(rng) * 2.0 * STD_AGE_PAIR - STD_AGE_PAIR;
        partnerAge = a + rand_offset;
    }
    return static_cast<int>(floor(partnerAge));
}

void Simulation::fledgeHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    if((*it)->isPaired() == false) {
        // First consider case of orphan living alone - no new household needed
        int currentHH = (*it)->getHousehold();
        int currentNhood = (*it)->getNeighborhood();
        int hhSize = static_cast<int>(allHosts.get<household_tag>().count(currentHH));
        if(hhSize == 1) {
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setFledge(true); });
        }
        else {
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setHousehold(hholdCtr); });
            int randNeighborhood = (int)(floor((double)NUM_NEIGHBORHOODS*r01(rng)));
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setNeighborhood(randNeighborhood); });
            int age = (*it)->getAgeInY();
            int infections;
            if((*it)->isInfected()) {
                for(int s = 0; s < INIT_NUM_STYPES; s++) {
                    infections = (*it)->isInfectedZ(s);
                    if(infections > 0) {
                        numInfecteds[age][s][currentNhood] -= infections;
                        numInfecteds[age][s][randNeighborhood] += infections;
                    }
                }
            }
            allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setFledge(true); });
            allHouseholds.insert(hholdCtr);
            hholdCtr++;
        } // end if come from household > 1
    }
}

void Simulation::birthHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    int hhold = (*it)->getHousehold();
    int nhood = (*it)->getNeighborhood();
    double DOB = t;
    Host * newHostPtr;
    newHostPtr = new Host(t, DOB, idCtr, hhold, nhood, currentEvents, simParsPtr, rng);
    allHosts.insert(boost::shared_ptr<Host>(newHostPtr));
    idCtr++;
}

void Simulation::infectHost(int id, int s) { // Could rewrite to take itr instead of id
    HostsByID::iterator it = allHosts.find(id);
    if(s < HFLU_INDEX || (*it)->isInfectedHflu() == false) {
        numInfecteds[(*it)->getAgeInY()][s][(*it)->getNeighborhood()]++;
    }
    if(s < HFLU_INDEX || (*it)->isInfectedHflu() == false) {
        allHosts.modify(it, [&](boost::shared_ptr<Host> h) { h->becomeInfected(s, t, currentEvents, rng); });
        allHosts.modify(it, [](boost::shared_ptr<Host> h) { h->setInf(true); });
    }
}

void Simulation::recoverHost(int id, int s) {
    HostsByID::iterator it = allHosts.find(id);

    // Host clears just one infection. All infections represented in numInfecteds[][].
    if(s != HFLU_INDEX || (s == HFLU_INDEX && (*it)->isInfectedZ(s) == 1)) {
        numInfecteds[(*it)->getAgeInY()][s][(*it)->getNeighborhood()]--;
    }
    allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->recover(s, t); });
    if((*it)->totStrains() == 0) {
        allHosts.modify(it, [=](boost::shared_ptr<Host> h) { h->setInf(false); });
    }
}

void Simulation::vaccinateHost(int id) {
    HostsByID::iterator it = allHosts.find(id);
    allHosts.modify(it, [](boost::shared_ptr<Host> h) { h->getVaccinated(); });
}

std::array<double, INIT_NUM_STYPES> Simulation::calculateSerotypePrevalenceRates()
{
    // Get population sizes of kids <5
    int N_total = 0;
    int I_total_pneumo = 0;
    std::array<int, INIT_NUM_STYPES> serotype_counts = {0};
    int I_total_hflu = 0;

    // Now count how many kids in each age group are infected -- note that not sufficient to use numInfecteds, which counts 'effective' number of infecteds 
    // (i.e., the number of infections) for purposes of calculating the force of infection
    int hostAge = 0;
    for(auto host : allHosts.get<age_tag>())
    {
        hostAge = host->getAgeInY();

        if(hostAge < 5)
        {
            N_total++;

            std::vector<int> concurrent_serotypes;

            for(int i = 0; i < INIT_NUM_STYPES; i++)
            {
                if(host->isInfectedZ(i))
                {
                    concurrent_serotypes.push_back(i);
                }
            }

            while(concurrent_serotypes.size() > 1)
            {
                auto random_serotype_index = static_cast<std::size_t>(r01(rng) * concurrent_serotypes.size());
                std::swap(concurrent_serotypes[random_serotype_index], concurrent_serotypes.back());
                concurrent_serotypes.erase(concurrent_serotypes.begin() + concurrent_serotypes.size() - 1);
            }

            if(!concurrent_serotypes.empty())
            {
                serotype_counts[concurrent_serotypes.front()]++;
                I_total_pneumo++;
            }

            I_total_hflu += host->isInfectedHflu();
        }
    } // end for each host

    std::cout << "\tThere are " << N_total << " kids <5 y old; " << I_total_pneumo << " (" << 100.0*(double)I_total_pneumo / (double)N_total << "%) carry pneumo and "
        << I_total_hflu << " (" << 100.0*(double)I_total_hflu / (double)N_total << "%) carry Hflu" << std::endl;

    std::array<double, INIT_NUM_STYPES> serotype_prevalence_rates;

    for(int i = 0; i < INIT_NUM_STYPES; i++)
    {
        serotype_prevalence_rates[i] = serotype_counts[i] / (double)N_total;
    }

    return serotype_prevalence_rates;
}


void Simulation::calcSI() {

#if defined( NO_HHOLDS ) && defined( NO_AGE_ASSORT ) // RANDOM MIXING: NO HOUSEHOLDS AND NO AGE-ASSORTATIVITY
    int N_total = 0;
    int nhood = 0;
    HostsByAge& sorted_index = allHosts.get<age_tag>();
    for(int n = 0; n < INIT_NUM_AGE_CATS; n++) {
        N_total += static_cast<int>(sorted_index.count(n));
    }
    int I_total[INIT_NUM_STYPES];
    for(int s = 0; s < INIT_NUM_STYPES; s++) {
        I_total[s] = 0;
    }
    for(int s = 0; s < INIT_NUM_STYPES; s++) {
        for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
            I_total[s] += numInfecteds[a][s][nhood];
        }
    }
    double prInf;
    double rTrans;
    double infectionTime;
    for(auto it : allHosts.get<age_tag>()) {
        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            rTrans = simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(I_total[s] - it->isInfectedZ(s)) / (double)(N_total - 1); // beta*I/N (excluding self)
            prInf = it->getSusc(s) * (rTrans + simParsPtr->get_serotypePar_ij(IMMIGRATION_INDEX, s));
            prInf = 1.0 - exp(-prInf * EPID_DELTA_T);
            if(r01(rng) < prInf) {
                infectionTime = (double)r01(rng) * (double)EPID_DELTA_T + (double)t;
                if(infectionTime <= t) {
                    std::cout << "\tAdding epsilon to event time." << std::endl;
                    infectionTime += pow(10, APPROX_NOW);
                }
                if(infectionTime < it->getDOD()) {
                    addEvent(infectionTime, Event::Type::Infection, it->getID(), s);
                }
            }
        } // end for each serotype
    } // end for each host
#endif 

#if defined( NO_HHOLDS ) && !defined( NO_AGE_ASSORT ) // AGE-ASSORTATIVE MIXING BUT NO HOUSEHOLDS
    // Get population ages and normalize age-assortative matrix, alpha
    int N_total = 0;
    int N_age[INIT_NUM_AGE_CATS];
    double alpha[INIT_NUM_AGE_CATS][INIT_NUM_AGE_CATS];
    HostsByAge& sorted_index = allHosts.get<age>();
    for(int n = 0; n < INIT_NUM_AGE_CATS; n++) {
        N_age[n] = sorted_index.count(n);
        N_total += N_age[n];
        for(int n2 = 0; n2 < INIT_NUM_AGE_CATS; n2++) {
            alpha[n][n2] = 0.0;
        }
    }
    double runningSum;
    for(int i = 0; i < INIT_NUM_AGE_CATS; i++) { // for donors of age j
        if(N_age[i] > 0) {
            runningSum = 0.0;
            for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
                if((N_age[j] > 1) || ((N_age[j] == 1) && (i != j))) { // only consider transmission b/w two age groups if both present and not just self
                    runningSum += simParsPtr->get_waifw_ij(i, j);
                }
            }
            if(runningSum > 0) {
                for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
                    if((N_age[j] > 1) || ((N_age[j] == 1) && (i != j))) {
                        alpha[i][j] = simParsPtr->get_waifw_ij(i, j) / runningSum;
                    }
                }
            }
        }
    }

    // Calculate force for every host
    int hostID;
    int hostAge;
    int numInfectionsZ;
    int nhood;
    double susc_z;
    double prInf;
    double rTrans;
    double infectionTime;
    for(auto host : allHosts.get<age>()) {
        hostID = host->getID();
        hostAge = host->getAgeInY();
        nhood = host->getNeighborhood();
        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            susc_z = host->getSusc(s);
            numInfectionsZ = host->isInfectedZ(s);
            rTrans = 0.0;
            for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                if(a != hostAge && N_age[a] > 0) {
                    rTrans += simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(numInfecteds[a][s][nhood]) / (double)N_age[a] * alpha[hostAge][a];
                }
                else if(N_age[a] > 1) { // no force from self if only member of cohort
                    rTrans += simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(numInfecteds[a][s][nhood] - numInfectionsZ) / (double)(N_age[a] - 1) * alpha[hostAge][a];
                }
            } // end for each age	
            assert(rTrans >= 0.0);
            prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij(IMMIGRATION_INDEX, s);
            prInf = 1.0 - exp(-prInf * EPID_DELTA_T);
            if(r01(rng) < prInf) {
                infectionTime = r01(rng) * EPID_DELTA_T + t;
                if(infectionTime <= t) {
                    std::cout << "\tAdding epsilon to event time." << std::endl;
                    infectionTime += pow(10, APPROX_NOW);
                }
                if(infectionTime < host->getDOD()) {
                    addEvent(infectionTime, Event::Type::Infection, hostID, s);
                }
            }
        } // end for each serotype
    } // end for each host
#endif


#if !defined( NO_HHOLDS ) && defined( NO_AGE_ASSORT ) // HOUSEHOLDS BUT RANDOM MIXING
    // Count total size of each neighborhood
    int neighborhoodSizes[NUM_NEIGHBORHOODS];
    int thisSize;
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        thisSize = 0;
        for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
            thisSize += allHosts.get<an>().count(boost::make_tuple(n, a));
        }
        neighborhoodSizes[n] = thisSize;
    } // end for each neighborhood

    // Calculate force for every host
    int hostID;
    int otherID;
    int otherHhold;
    int hostHhold;
    int hostNhood;
    int hholdSize;
    int numInfectionsZ;
    int numHholdInfecteds;
    double numNhoodInfecteds;
    double susc_z;
    double prInf;
    double rTrans;
    double infectionTime;
    double n_weight;
    for(auto host : allHosts.get<age>()) {
        hostID = host->getID();
        hostHhold = host->getHousehold();
        hostNhood = host->getNeighborhood();
        hholdSize = allHosts.get<household>().count(hostHhold) - 1; // do not count self in household 

        for(int s = 0; s < INIT_NUM_STYPES; s++) {
            susc_z = (*it)->getSusc(s);
            numInfectionsZ = (*it)->isInfectedZ(s);
            numHholdInfecteds = 0;
            rTrans = 0.0;

            if(hholdSize > 0) {
                auto pit = allHosts.get<household>().equal_range(hostHhold);
                for(auto fit = pit.first; fit != pit.second; fit++) {
                    numHholdInfecteds += (*fit)->isInfectedZ(s); // adding self in numerator (easier than checking family members' IDs)
                } // end for each host in household
                rTrans += RHO_H * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(numHholdInfecteds - numInfectionsZ) / (double)hholdSize; // not counting self in denominator
            } // end for household size > 0

            // Do non-household contacts, weighted for each contacted neighborhood
            for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
                n_weight = simParsPtr->get_normalized_neighbor(hostNhood, n);
                if(n_weight > 0) {
                    numNhoodInfecteds = 0;
                    for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                        numNhoodInfecteds += (double)numInfecteds[a][s][n];
                    } // end for this age
                    if(n != hostNhood) {
                        rTrans += (1.0 - RHO_H) * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)n_weight*(numNhoodInfecteds) / (double)(neighborhoodSizes[n]);
                    }
                    else {
                        rTrans += (1.0 - RHO_H) * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)n_weight*(numNhoodInfecteds - numHholdInfecteds) / (double)(neighborhoodSizes[n] - 1 - hholdSize);
                    }
                } // end for if n_weight > 0 
            } // end for each neighborhood

            prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij(IMMIGRATION_INDEX, s);
            prInf = 1.0 - exp(-prInf * EPID_DELTA_T);
            if(r01(rng) < prInf) {
                infectionTime = r01(rng) * EPID_DELTA_T + t;
                if(infectionTime <= t) {
                    cout << "\tAdding epsilon to event time." << endl;
                    infectionTime += pow(10, APPROX_NOW);
                }
                if(infectionTime < (*it)->getDOD()) {
                    addEvent(infectionTime, INFECTION_EVENT, hostID, s);
                }
            } // end for if r01(rng)<prInf
        } // end for each serotype
    } // end for each host
#endif


#if !defined( NO_HHOLDS ) && !defined( NO_AGE_ASSORT ) // HOUSEHOLDS AND AGE-ASSORTATIVE MIXING
    int N_total[NUM_NEIGHBORHOODS];
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        N_total[n] = 0;
    }
    int N_age[INIT_NUM_AGE_CATS][NUM_NEIGHBORHOODS];
    double alpha[INIT_NUM_AGE_CATS][INIT_NUM_AGE_CATS][NUM_NEIGHBORHOODS];
    HostsByAge& sorted_index = allHosts.get<age>();
    for(int h = 0; h < NUM_NEIGHBORHOODS; h++) {
        for(int n = 0; n < INIT_NUM_AGE_CATS; n++) {
            N_age[n][h] = allHosts.get<an>().count(boost::make_tuple(h, n));
            N_total[h] += N_age[n][h];
            for(int n2 = 0; n2 < INIT_NUM_AGE_CATS; n2++) {
                alpha[n][n2][h] = 0.0;
            }
        } // end for age n
    } // end for neighborhood h

    // Calculate alpha_nh for each neighborhood
    double runningSum;
    for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
        for(int i = 0; i < INIT_NUM_AGE_CATS; i++) { // for donors of age j
            if(N_age[i][n] > 0) {
                runningSum = 0.0;
                for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
                    if((N_age[j][n] > 1) || ((N_age[j][n] == 1) && (i != j))) { // only consider transmission b/w two age groups if both present
                        runningSum += simParsPtr->get_waifw_ij(i, j);
                    }
                }
                if(runningSum > 0) {
                    for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
                        if((N_age[j][n] > 1) || ((N_age[j][n] == 1) && (i != j))) {
                            alpha[i][j][n] = simParsPtr->get_waifw_ij(i, j) / runningSum;
                        }
                    }
                }
            }
        } // end for donors of age j
    } // end for this neighborhood n

    // Calculate force for every host
    int hostID;
    int otherID;
    int hostHhold;
    int otherHhold;
    int hostNhood;
    int hholdSize;
    int hostAge;
    int numInfectionsZ;
    int numHholdInfecteds;
    double numNhoodInfecteds;
    double n_weight;
    double alpha_hh[INIT_NUM_AGE_CATS]; // alpha_ij for this host i
    double susc_z;
    double prInf;
    double rTrans;
    double infectionTime;
    double sumAlphas;
    int hholdAges[INIT_NUM_AGE_CATS];
    for(auto host : allHosts.get<age>()) {
        hostID = host->getID();
        hostHhold = host->getHousehold();
        hostNhood = host->getNeighborhood();
        hostAge = host->getAgeInY();
        hholdSize = allHosts.get<household>().count(hostHhold) - 1; // do not count self in household
        sumAlphas = 0.0;
        for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
            hholdAges[a] = 0;
            alpha_hh[a] = 0;
        }

        // Update alpha_hh vector for household members
        if(hholdSize > 0) {
            for(int a = 0; a < INIT_NUM_AGE_CATS; a++) { // ...for each age 
                auto pit = allHosts.get<ah>().equal_range(boost::make_tuple(hostHhold, a));
                for(auto fit = pit.first; fit != pit.second; fit++) {
                    if((*fit)->getID() != hostID) {
                        hholdAges[a]++;
                    }
                } // end for each household member of this age
                if(hholdAges[a] > 0) { // if there's a non-self household member of this age
                    sumAlphas += simParsPtr->get_waifw_ij(hostAge, a);
                }
            } // end for each age
            if(sumAlphas > 0) {
                for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                    if(hholdAges[a] > 0) {
                        alpha_hh[a] = simParsPtr->get_waifw_ij(hostAge, a) / sumAlphas;
                    }
                }
            }
        } // end for if hholdSize > 0

        for(int s = 0; s < INIT_NUM_STYPES; s++) { // for each serotype
            susc_z = (*it)->getSusc(s);
            rTrans = 0.0;
            numInfectionsZ = (*it)->isInfectedZ(s);
            for(int a = 0; a < INIT_NUM_AGE_CATS; a++) {
                numHholdInfecteds = 0;
                if(hholdSize > 0) {
                    auto pit = allHosts.get<ah>().equal_range(boost::make_tuple(hostHhold, a));
                    for(auto fit = pit.first; fit != pit.second; fit++) {
                        numHholdInfecteds += (*fit)->isInfectedZ(s); // note that this number includes the host
                    } // end for all hosts in household of this age
                    if(hholdAges[a] > 0 && a == hostAge) { // Calculate household transmission for this age
                        rTrans += RHO_H * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(numHholdInfecteds - numInfectionsZ) / (double)hholdAges[a] * alpha_hh[a];
                    }
                    else if(hholdAges[a] > 0 && a != hostAge) {
                        rTrans += RHO_H * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * (double)(numHholdInfecteds) / (double)hholdAges[a] * alpha_hh[a];
                    }
                } // end if hholdSize > 0

                for(int n = 0; n < NUM_NEIGHBORHOODS; n++) {
                    n_weight = simParsPtr->get_normalized_neighbor(hostNhood, n);
                    if(n_weight > 0) {
                        numNhoodInfecteds = numInfecteds[a][s][n];
                        if(n == hostNhood && a == hostAge && (N_age[a][n] > hholdAges[a] + 1)) { // same neighborhood and same age as host
                            rTrans += (1.0 - RHO_H) * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * n_weight*(double)(numNhoodInfecteds - numHholdInfecteds) / (double)(N_age[a][n] - hholdAges[a] - 1) * alpha[hostAge][a][n];
                        }
                        else if(n == hostNhood && a != hostAge && (N_age[a][n] > hholdAges[a])) { // same neighborhood but different age
                            rTrans += (1.0 - RHO_H) * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * n_weight*(double)(numNhoodInfecteds - numHholdInfecteds) / (double)(N_age[a][n] - hholdAges[a]) * alpha[hostAge][a][n];
                        }
                        else if(n != hostNhood && N_age[a][n] > 0) { // different neighborhood, but need people there
                            rTrans += (1.0 - RHO_H) * simParsPtr->get_serotypePar_ij(BETA_INDEX, s) * n_weight*(double)(numNhoodInfecteds) / (double)(N_age[a][n]) * alpha[hostAge][a][n];
                        }
                    } // end if n_weight > 0
                } // end for each neighborhood
            } // end for every age

            prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij(IMMIGRATION_INDEX, s);
            prInf = 1.0 - exp(-prInf * EPID_DELTA_T);
            if(r01(rng) < prInf) {
                infectionTime = r01(rng) * EPID_DELTA_T + t;
                if(infectionTime < (*it)->getDOD()) {
                    addEvent(infectionTime, INFECTION_EVENT, hostID, s);
                }
            }
        } // end for each serotype
    } // end for each host
#endif
}

std::string Simulation::d2str(double d) {
    std::stringstream t;
    t << d;
    return t.str();
}

std::string Simulation::makeName(std::string suffix) {
    std::string thisCtr = d2str(simID);
    std::string thisTr = d2str(treatment);
    std::string thisName = "tr_" + thisTr + "_sim_" + thisCtr + "_" + suffix;
    return thisName;
}

std::string Simulation::makeBigName(std::string suffix, int index) {
    std::string thisName = "tr_" + d2str(treatment) + "_sim_" + d2str(simID) + "_" + suffix + "_" + d2str(index);
    return thisName;
}

std::string Simulation::makeBiggerName(std::string suffix1, int index1, std::string suffix2, int index2) {
    std::string thisName = "tr_" + d2str(treatment) + "_sim_" + d2str(simID) + "_" + suffix1 + "_" + d2str(index1) + "_" + suffix2 + "_" + d2str(index2);
    return thisName;
}

void Simulation::addEvent(double et, Event::Type event_type, int hid, int s) {
    while(currentEvents.find(Event(et)) != currentEvents.end()) {
        et += pow(10, APPROX_NOW);
        std::cout << "\tSimulation is adjusting event time to prevent collision (host id " << hid << ", event id " << int(event_type) << ", strain " << s << ", event time " << et << ").\n";
    }
    Event thisEvent(et, event_type, hid, s);
    currentEvents.insert(thisEvent);
}
