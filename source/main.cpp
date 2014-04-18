#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Simulation.h"
#include "SimPars.h"

#define HFLU_INDEX (INIT_NUM_STYPES-1)

// Function prototypes
void adjustTreatment(int treatmentNumber, double treatment, int simNumber);
void initializeBeta(int treatmentNumber, double beta, int simNumber);
void printAssumptions();
double adjustBeta(double preve, double weight, double sb, int treatmentNumber, int simNumber);
std::string d2str(double d);
std::string makeName(int treatmentIdx, int simIdx, std::string suffix);
void printTotalTime(time_t t1, time_t t2);

void match_prevalence(int treatmentNumber, int simNumber, double treatment, double startingBeta)
{
    std::cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << std::endl;
    int matchAttempts = 0;
    double error = 0.0;
    double usedBeta = startingBeta;
    double prevError = 10.0;
    double oldPrevError = 1.0;
    double thisBeta = 0.0;
    double weight = INIT_WEIGHT;
    while(matchAttempts < MAX_MATCH_ATTEMPTS && abs(prevError) > PREV_ERROR_THOLD) {
        SimPars thesePars(treatmentNumber, simNumber);
        SimPars * spPtr = &thesePars;
        Simulation thisSim(treatmentNumber, simNumber, spPtr);
        std::cout << "  Attempt #" << matchAttempts + 1 << std::endl;
        if(prevError == 10.0) {
            usedBeta = startingBeta;
        }
        thisSim.runDemSim();
        prevError = thisSim.runTestEpidSim() - TARGET_PREV;
        std::cout << "Prevalence error=" << prevError << "; weight=" << weight << std::endl;
        if(abs(prevError) > PREV_ERROR_THOLD) {
            if(prevError*oldPrevError < 0) { // if error changed signs and overshot, reduce weight
                weight *= COOL_DOWN;
            }
            else if(abs(TARGET_PREV - prevError) / abs(TARGET_PREV - oldPrevError) > TEMP_THOLD) { // if climbing too slowly, increase weight
                weight *= WARM_UP;
            }
            thisBeta = adjustBeta(prevError, weight, startingBeta, treatmentNumber, simNumber);
            usedBeta = thisBeta;
            oldPrevError = prevError;
        }
        matchAttempts++;
    }
    error = prevError;

    if(prevError < PREV_ERROR_THOLD) {
        std::cout << "\tBeginning simulation #" << simNumber << std::endl;
        time_t tic;
        tic = time(NULL);
        SimPars thesePars(treatmentNumber, simNumber);
        SimPars * spPtr = &thesePars;
        Simulation thisSim(treatmentNumber, simNumber, spPtr);
        thisSim.runDemSim();
        thisSim.runEpidSim();
        time_t toc;
        toc = time(NULL);
        printTotalTime(tic, toc);
    }
    else {
        std::cout << "Acceptable prevalence not found." << std::endl;
    }

    // Output used beta
    std::ofstream thisBetaStream;
    std::ofstream errorStream;
    std::string betaFile = makeName(treatmentNumber, simNumber, "beta_used");
    std::string errorFile = makeName(treatmentNumber, simNumber, "prevalence_errors");
    thisBetaStream.open("../../outputs/" + betaFile, std::ios::out);
    errorStream.open("../../outputs/" + errorFile, std::ios::out);
    thisBetaStream << treatment << "\t" << usedBeta;
    errorStream << error;
    thisBetaStream.close();
    errorStream.close();
}

void run_simulation(int simNumber, int treatmentNumber)
{
    std::cout << "\tBeginning simulation #" << simNumber << std::endl;
    time_t tic;
    tic = time(NULL);
    SimPars thesePars(treatmentNumber, simNumber);
    SimPars * spPtr = &thesePars;
    Simulation thisSim(treatmentNumber, simNumber, spPtr);
    thisSim.runDemSim();
    thisSim.runEpidSim();
    time_t toc;
    toc = time(NULL);
    printTotalTime(tic, toc);
}

int main()
{
    int treatmentNumber = 1;
    double treatment = 1;
    int simNumber = 1;
    std::cout << "Treatment number " << treatmentNumber << ", treatment value " << treatment << ", simulation number " << simNumber << std::endl;
    printAssumptions();

    double startingBeta;
    double betaTable[100][100]; // buffer

    // Get beta values
    std::ifstream thisFile;
    thisFile.open("../../inputs/Betas_used.txt", std::ios::in);
    if(!thisFile) {
        std::cerr << "Error reading Betas_used.txt." << std::endl;
        exit(1);
    }
    double thisVal;
    int b = 0;
    while(!thisFile.eof()) { // Should instead grab beta closest to treatment value in Betas_used
        thisFile >> thisVal;
        betaTable[b][1] = thisVal;
        b++;
    }

    // Get corresponding treatment value (sigma)
    std::ifstream thisFileT;
    thisFileT.open("../../inputs/Treatments.txt", std::ios::in);
    if(!thisFileT) {
        std::cerr << "Error reading Treatments.txt." << std::endl;
        exit(1);
    }
    double thisValT;
    int bT = 0;
    double diff = 10000.0;
    int bestTreatment = 0;
    while(!thisFileT.eof()) {
        thisFileT >> thisValT;
        if(abs(treatment - thisValT) < diff) { // Finds index of previous treatment most similar to this one
            bestTreatment = bT;
            diff = abs(treatment - thisValT);
        }
        betaTable[bT][0] = thisValT;
        bT++;
    }

    startingBeta = betaTable[bestTreatment][1];
    std::cout << "Best match to treatment value " << treatment << " is beta = " << startingBeta << " (associated treatment value is " << betaTable[bestTreatment][0] << ")" << std::endl;
    initializeBeta(treatmentNumber, startingBeta, simNumber);
    adjustTreatment(treatmentNumber, treatment, simNumber);

    //match_prevalence(treatmentNumber, simNumber, treatment, startingBeta);
    run_simulation(simNumber, treatmentNumber);
}

double adjustBeta(double preve, double w, double sb, int treatmentNumber, int simNumber) {
    double newBetas[INIT_NUM_STYPES];
    for(int b = 0; b < INIT_NUM_STYPES; b++) {
        newBetas[b] = 0.0;
    }

    // Either write starting betas or adjust old betas
    if(preve == 10.0) {
        for(int b = 0; b < HFLU_INDEX; b++) {
            newBetas[b] = sb;
        }
        newBetas[HFLU_INDEX] = HFLU_BETA;
    }
    else {
        std::ifstream thisFile3;
        std::string filename = makeName(treatmentNumber, simNumber, "BETA");
        thisFile3.open("../../inputs" + filename, std::ios::in);
        if(!thisFile3) {
            std::cerr << "Error reading " << filename << std::endl;
            exit(1);
        }
        double thisVal;
        int b = 0;
        while(!thisFile3.eof()) {
            thisFile3 >> thisVal;
            std::cout << "beta=" << thisVal << ";";
            if(b < HFLU_INDEX) {
                newBetas[b] = thisVal*(1.0 - w*preve);
            }
            else {
                newBetas[b] = HFLU_BETA;
            }
            std::cout << "newBetas[" << b << "]=" << newBetas[b] << std::endl;
            b++;
        }
        thisFile3.close();
    }

    std::ofstream betaStream;
    std::string filename2 = makeName(treatmentNumber, simNumber, "BETA");
    betaStream.open("../../outputs/" + filename2, std::ios::out);
    for(int b = 0; b < INIT_NUM_STYPES; b++) {
        betaStream << newBetas[b] << "\t";
    }
    betaStream.close();

    double returnBeta = newBetas[0];
    return returnBeta;
}

void initializeBeta(int treatmentNumber, double beta, int simNumber) {
    std::ofstream betaStream;
    std::string filename2 = makeName(treatmentNumber, simNumber, "BETA");
    betaStream.open("../../outputs/" + filename2, std::ios::out);
    for(int b = 0; b < HFLU_INDEX; b++) {
        betaStream << beta << "\t";
    }
    betaStream << HFLU_BETA;
    betaStream.close();
}

void adjustTreatment(int treatmentNumber, double treatment, int simNumber) {
    double thisXI[INIT_NUM_STYPES][INIT_NUM_STYPES];
    for(int i = 0; i < INIT_NUM_STYPES; i++) {
        for(int j = 0; j < INIT_NUM_STYPES; j++) {
            thisXI[i][j] = 0;
        }
    }

    for(int i = 0; i < HFLU_INDEX; i++) {
        thisXI[i][i] = treatment;
    }
    thisXI[HFLU_INDEX][HFLU_INDEX] = HFLU_SIGMA;

    std::ofstream xiStream;
    std::string XIFile = makeName(treatmentNumber, simNumber, "XI");
    xiStream.open("../../outputs/" + XIFile, std::ios::out);
    for(int i = 0; i < INIT_NUM_STYPES; i++) {
        for(int j = 0; j < INIT_NUM_STYPES; j++) {
            xiStream << thisXI[i][j];
            if((j + 1) % INIT_NUM_STYPES == 0) {
                xiStream << std::endl;
            }
            else {
                xiStream << "\t";
            }
        }
    }
    xiStream.close();
}

std::string d2str(double d) {
    std::stringstream t;
    t << d;
    return t.str();
}

std::string makeName(int treatmentIdx, int simIdx, std::string suffix) {
    std::string thisCtr = d2str(simIdx);
    std::string thisTr = d2str(treatmentIdx);
    std::string thisName = "tr_" + thisTr + "_sim_" + thisCtr + "_" + suffix;
    return thisName;
}

void printAssumptions() {
    std::cout << "Model assumptions:" << std::endl;
#ifdef NO_HHOLDS
    std::cout << "\t--no household structure" << std::endl;
#else 
    cout << "\t--households present" << endl;
#endif

#ifdef NO_AGE_ASSORT
    std::cout << "\t--no age-assortative mixing" << std::endl;
#else
    cout << "\t--age-assortative mixing" << endl;
#endif

#ifdef MATCH_PREVALENCE
    cout << "\t--matching prevalence" << endl;
#else
    std::cout << "\t--not matching prevalence" << std::endl;
#endif
}

void printTotalTime(time_t t1, time_t t2) {
    auto time_difference = static_cast<int>(t2 - t1);
    if(time_difference < 60) {
        std::cout << "  Total simulation time: " << time_difference << " seconds" << std::endl;
    }
    else if(time_difference < 3600) {
        int nmin = static_cast<int>(floor(time_difference / 60.0));
        std::cout << "  Total simulation time: " << nmin << " min " << time_difference - nmin * 60 << " seconds " << std::endl;
    }
    else {
        int nhour = static_cast<int>(floor(time_difference / 3600.0));
        int nmin = static_cast<int>(floor((time_difference - nhour * 3600) / 60));
        int nsec = static_cast<int>(time_difference - nhour * 3600 - nmin * 60);
        std::cout << "  Total simulation time: " << nhour << " h " << nmin << " m " << nsec << " s " << std::endl;
    }
}

