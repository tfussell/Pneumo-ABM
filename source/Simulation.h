#pragma once

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/random.hpp>

#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Containers.h"
#include "SimPars.h"

class Simulation
{
public:
    Simulation(int trt, int sid, SimPars * spPtr);
    ~Simulation();

    // MEMBER FUNCTION PROTOTYPES
    void runDemSim();
    void runEpidSim();
    std::array<double, NUM_STYPES> runTestEpidSim();

private:
    // SIMULATION OBJECTS
    boost::mt19937 rng;
    double t;
    double demOutputStrobe;
    double epidOutputStrobe;
    double demComplete;
    int idCtr;
    int hholdCtr;
    int simID;
    int eventCtr; // for memory management and debugging
    HostContainer allHosts; // shared_ptr to Hosts
    HouseholdContainer allHouseholds; // set of household ids
    // counts number infected with each serotype by age
    std::array<std::array<std::array<int, NUM_NEIGHBORHOODS>, NUM_STYPES>, INIT_NUM_AGE_CATS> numInfecteds;
    EventQueue currentEvents; // current events queue
    SimPars * simParsPtr;
    int treatment;

    // SIMULATION FUNCTIONS
    void ageHost(int id);
    int calcPartnerAge(int a);
    void executeEvent(Event & te);
    void seedInfections();
    void killHost(int id);
    void pairHost(int id);
    void partner2Hosts(int id1, int id2);
    void fledgeHost(int id);
    void birthHost(int id);
    void infectHost(int id, int s);
    void recoverHost(int id, int s);
    void vaccinateHost(int id, const std::string &vaccine);
    void calcSI();
    std::string d2str(double d);
    std::string makeName(std::string s);
    std::string makeBigName(std::string, int);
    std::string makeBiggerName(std::string, int, std::string, int);
    void addEvent(double et, Event::Type event_type, int hid, int s);
    std::array<double, NUM_STYPES> calculateSerotypePrevalenceRates();

    // STREAM MANAGEMENT
    void writeDemOutput();
    void writeEpidOutput();
    void initOutput();
    void closeOutput();
    void writeTheta();

    std::ofstream ageDistStream;
    std::ofstream hhDistStream;
    std::ofstream demTimesStream;
    std::ofstream infectionStream;
    std::ofstream epidTimesStream;
    std::ofstream infectedStream;
    std::ofstream coinfectionHistStream;
    std::ofstream coinfectionHFHistStream;
    std::ofstream thetaStream;
    std::ofstream totCarriageStream;

    std::string ageDistFile;
    std::string hhDistFile;
    std::string demTimesFile;
    std::string epidTimesFile;
    std::string coinfectionHistFile;
    std::string coinfectionHFHistFile;
    std::string thetaFile;
    std::string totCarriageFile;
};
