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
void printAssumptions();
std::string d2str(double d);
std::string makeName(int treatmentIdx, int simIdx, std::string suffix);
void printTotalTime(time_t t1, time_t t2);

double log_factorial(double n)
{
    return n * log(n) - n + (log(n * (1 + 4 * n * (1 + 2 * n)))) / 6 + log(3.141592653589793238462) / 2;
}

double calculate_likelihood(const std::array<double, INIT_NUM_STYPES + 1> &expected, const std::array<int, INIT_NUM_STYPES + 1> &observed)
{
    int n = std::accumulate(observed.begin(), observed.end(), 0);
    double log_likelihood = log_factorial(n);

    for(int i = 0; i < INIT_NUM_STYPES + 1; i++)
    {
        if(i != HFLU_INDEX)
        {
            log_likelihood -= log_factorial(observed[i]);
            log_likelihood += observed[i] * log(expected[i]);
        }
    }

    return log_likelihood;
}

void adjust_ranks(const std::array<double, INIT_NUM_STYPES> &errors, std::array<double, INIT_NUM_STYPES> &ranks)
{
    for(int i = 0; i < INIT_NUM_STYPES - 1; i++)
    {
        if(errors[i] > 0)
        {
            ranks[i] = std::min(ranks[i] + 0.5, (double)(INIT_NUM_STYPES - 1));
        }
        else
        {
            ranks[i] = std::max(ranks[i] - 0.5, 1.0);
        }
        std::cout << "newRanks[" << i << "]=" << ranks[i] << std::endl;
    }
}

void adjustBeta(double preve, double w, std::array<double, INIT_NUM_STYPES> &betas)
{
    for(int b = 0; b < INIT_NUM_STYPES; b++)
    {
        if(b != HFLU_INDEX)
        {
            betas[b] *= (1.0 - w*preve);
        }
        else
        {
            betas[b] = HFLU_BETA;
        }
        std::cout << "newBetas[" << b << "]=" << betas[b] << std::endl;
    }
}

void match_prevalence(int treatmentNumber, int simNumber, double treatment, double startingBeta)
{
    std::cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << std::endl;

    std::array<double, INIT_NUM_STYPES> betas;
    betas.fill(startingBeta);

    std::array<double, INIT_NUM_STYPES> serotype_ranks;
    for(int i = 0; i < INIT_NUM_STYPES; i++)
    {
        serotype_ranks[i] = i + 1;
    }

    double best_likelihood = std::numeric_limits<double>::lowest();
    std::array<double, INIT_NUM_STYPES> best_betas = betas;
    std::array<double, INIT_NUM_STYPES> best_serotype_ranks = serotype_ranks;

    double prevError = 10.0;
    double oldPrevError = 1.0;
    double weight = INIT_WEIGHT;
    std::array<int, INIT_NUM_STYPES + 1> observed_prevalence = {283, 237, 184, 117, 90, 85, 84, 70, 56, 54, 53, 51, 49, 49, 43, 38, 34, 34, 29, 25, 23, 21, 19, 18, 15, 0, 10079};

    for(int attempt = 0; attempt < 10; attempt++)
    {
        SimPars thesePars(treatmentNumber, simNumber);
        thesePars.set_serotype_ranks(serotype_ranks);
        thesePars.set_betas(betas);
        Simulation thisSim(treatmentNumber, simNumber, &thesePars);

        thisSim.runDemSim();
        auto serotype_counts = thisSim.runTestEpidSim();

        std::array<double, INIT_NUM_STYPES + 1> expected_prevalence = {0};
        std::array<double, INIT_NUM_STYPES> prevalence_error = {0};

        double expected_population = std::accumulate(serotype_counts.begin(), serotype_counts.end(), 0.0);
        double infected_prevalence = 0;
        int observed_population = std::accumulate(observed_prevalence.begin(), observed_prevalence.begin() + INIT_NUM_STYPES + 1, 0);
        double target_prevalence = std::accumulate(observed_prevalence.begin(), observed_prevalence.begin() + INIT_NUM_STYPES, 0.0) / observed_population;

        for(int i = 0; i < INIT_NUM_STYPES; i++)
        {
            infected_prevalence += serotype_counts[i] / (double)expected_population;
            prevalence_error[i] = expected_prevalence[i] - observed_prevalence[i] / (double)observed_population;

            if(i == HFLU_INDEX)
            {
                prevalence_error[i] = 0;
            }

            expected_prevalence[i] = (serotype_counts[i] + 0.01) / (expected_population + INIT_NUM_STYPES * 0.01);
        }

        expected_prevalence.back() = 1 - infected_prevalence;

        double likelihood = calculate_likelihood(expected_prevalence, observed_prevalence);

        if(likelihood > best_likelihood)
        {
            best_likelihood = likelihood;
            best_betas = betas;
            best_serotype_ranks = serotype_ranks;
        }
        else
        {
            betas = best_betas;
            serotype_ranks = best_serotype_ranks;
        }

        std::cout << "Likelihood=" << likelihood << "(best=" << best_likelihood << ")" << std::endl;

        if(attempt % 2 == 0)
        {
            std::cout << "Fitting overall prevalence" << std::endl;

            prevError = infected_prevalence - target_prevalence;
            std::cout << "Observed prevalence=" << target_prevalence << "; expected prevalence=" << infected_prevalence << "; error=" << prevError << "; weight = " << weight << std::endl;

            if(abs(prevError) > PREV_ERROR_THOLD)
            {
                // if error changed signs and overshot, reduce weight
                if(prevError * oldPrevError < 0)
                {
                    weight *= COOL_DOWN;
                }
                // if climbing too slowly, increase weight
                else if(abs(target_prevalence - prevError) / abs(target_prevalence - oldPrevError) > TEMP_THOLD)
                {
                    weight *= WARM_UP;
                }

                adjustBeta(prevError, weight, betas);

                oldPrevError = prevError;
            }
        }
        else
        {
            std::cout << "Fitting per-serotype prevalence" << std::endl;
            adjust_ranks(prevalence_error, serotype_ranks);
        }
    }
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
    adjustTreatment(treatmentNumber, treatment, simNumber);

    match_prevalence(treatmentNumber, simNumber, treatment, 0.1);
    //run_simulation(simNumber, treatmentNumber);
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

