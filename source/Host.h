#pragma once

#include "Event.h"
#include "Infection.h"
#include "Parameters.h"
#include "SimPars.h"
#include <boost/random.hpp>

class Host
{
public:
    Host(double, double, int, int, int, EventQueue &, SimPars *, boost::mt19937 &);
    ~Host();

    // SET FUNCTION PROTOTYPES
    void setHousehold(int);
    void setNeighborhood(int);
    void setGroup(int);
    void setDOB(double);
    void incrementAge();
    void setPartner(int);
    void setFledge(bool);
    void setInf(bool);
    void getVaccinated(const std::string &vaccine);

    // GET FUNCTION PROTOTYPES
    int getAgeInY() const;
    int getID() const;
    int getGroup() const;
    int getHousehold() const;
    int getPartner() const;
    bool isPaired() const;
    bool isAdult() const;
    bool isEligible() const;
    bool hasFledged() const;
    int isInfectedZ(int) const;
    bool isInfected() const;
    bool isInfectedPneumo() const;
    bool isInfectedHflu() const;
    int totStrains() const;
    double getSusc(int) const;
    double getDOD() const;
    int getNeighborhood() const;
    int getSummedTheta() const;
    double getDOB() const;

    // MEMBER FUNCTION PROTOTYPES
    void calcLifeHist(double, EventQueue &, double, boost::mt19937 &);
    double calcDeath(boost::mt19937 &);
    void addEvent(double, Event::Type type, int, EventQueue &);
    void addEvent(double, Event::Type type, int, int, EventQueue &);
    double calcFledge(boost::mt19937 &);
    double calcPairAge(boost::mt19937 &);
    int calcNumBirths(boost::mt19937 &);
    double calcBirthAge(boost::mt19937 &);
    void becomeInfected(int, double, EventQueue &, boost::mt19937 &);
    void recover(int, double);
    double calcRecovery(int, double, double, boost::mt19937 &);
    double calcRecovery(int, double, Infection &, boost::mt19937 &);
    void calcSusc(double);

private:
    double DOB;
    double DOD;
    int age;
    int id;
    int group;
    int household;
    int neighborhood;
    int partner;
    bool fledge;
    InfectionMap carriage;
    int immune[NUM_STYPES];
    int carriageSummary[NUM_STYPES];
    double susc[NUM_STYPES];
    bool inf;
    SimPars * simParsPtr;
	std::string activeVaccine;
};
