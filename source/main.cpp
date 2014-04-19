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
void adjustBetas(const std::array<double, INIT_NUM_STYPES> &prevalence_error);
std::string d2str(double d);
std::string makeName(int treatmentIdx, int simIdx, std::string suffix);
void printTotalTime(time_t t1, time_t t2);

void match_prevalence(int treatmentNumber, int simNumber, double treatment, double startingBeta)
{
    std::cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << std::endl;

    int matchAttempts = 0;
    double sum_errors = 10.0;

    std::array<double, INIT_NUM_STYPES> target_prevalence, prevalence_error;
    target_prevalence.fill(0.01);
    prevalence_error.fill(0);

    while(matchAttempts < MAX_MATCH_ATTEMPTS && abs(sum_errors) > PREV_ERROR_THOLD) {
        SimPars thesePars(treatmentNumber, simNumber);
        Simulation thisSim(treatmentNumber, simNumber, &thesePars);
        std::cout << "  Attempt #" << matchAttempts + 1 << std::endl;

        thisSim.runDemSim();
        auto current_prevalence = thisSim.runTestEpidSim();

        for(int i = 0; i < INIT_NUM_STYPES; i++)
        {
            prevalence_error[i] = current_prevalence[i] - target_prevalence[i];
            if(i == INIT_NUM_STYPES - 1)
            {
                prevalence_error[i] = 0;
            }
        }

        sum_errors = std::accumulate(prevalence_error.begin(), prevalence_error.end(), 0.0, [](double prev, double v) { return prev + abs(v); });
        std::cout << "Prevalence error=" << sum_errors << std::endl;

        if(sum_errors > PREV_ERROR_THOLD) {
            adjustBetas(prevalence_error);
        }
        matchAttempts++;
    }

    if(sum_errors < PREV_ERROR_THOLD) {
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

    /*
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
    */
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
        throw std::runtime_error("");
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
        throw std::runtime_error("");
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

    match_prevalence(treatmentNumber, simNumber, treatment, startingBeta);
    //run_simulation(simNumber, treatmentNumber);
}

void adjustBetas(const std::array<double, INIT_NUM_STYPES> &prevalence_error) {
    std::array<double, INIT_NUM_STYPES> betas;

    {
        std::ifstream betas_in;
        std::string filename = makeName(1, 1, "BETA");
        betas_in.open("../../outputs/" + filename, std::ios::in);
        if(!betas_in) {
            std::cerr << "Error reading " << filename << std::endl;
            throw std::runtime_error("");
        }
        double thisVal;
        int b = 0;
        while(!betas_in.eof()) {
            betas_in >> thisVal;
            if(b < HFLU_INDEX) {
                if(prevalence_error[b] < 0)
                {
                    betas[b] = std::min(1.0, thisVal + 0.01);
                }
                else
                {
                    betas[b] = std::max(0.0, thisVal - 0.01);
                }
            }
            else {
                betas[b] = HFLU_BETA;
            }
            std::cout << "beta[" << b << "]: " << thisVal << "->" << betas[b] << std::endl;
            b++;
        }
    }

    std::ofstream betas_out;
    std::string filename2 = makeName(1, 1, "BETA");
    betas_out.open("../../outputs/" + filename2, std::ios::out);
    for(int b = 0; b < INIT_NUM_STYPES; b++) {
        betas_out << betas[b];
        if(b != INIT_NUM_STYPES - 1)
        {
            betas_out << "\t";
        }
    }
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
    std::cout << "\t--matching prevalence" << std::endl;
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

