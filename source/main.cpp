#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <thread>

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

double calculate_likelihood(const std::array<double, INIT_NUM_STYPES> &expected, const std::array<int, INIT_NUM_STYPES + 1> &observed)
{
    int n = std::accumulate(observed.begin(), observed.end(), 0);
    double log_likelihood = log_factorial(n);

    double proportion_uninfected = 1 + (0.0001 * (INIT_NUM_STYPES - 1));
    proportion_uninfected -= std::accumulate(expected.begin(), expected.end(), 0.0);

    for(int i = 0; i < INIT_NUM_STYPES + 1; i++)
    {
        if(i == HFLU_INDEX)
        {
            continue;
        }

        log_likelihood -= log_factorial(observed[i]);
        log_likelihood += observed[i] * log(i == INIT_NUM_STYPES ? proportion_uninfected : (expected[i] + 0.0001));
    }

    return log_likelihood;
}

void adjust_ranks(const std::array<double, INIT_NUM_STYPES> &errors, std::array<double, INIT_NUM_STYPES> &ranks)
{
    for(int i = 0; i < INIT_NUM_STYPES - 1; i++)
    {
        if(errors[i] > 0)
        {
            ranks[i] = std::min(ranks[i] + 0.1, (double)(INIT_NUM_STYPES - 1));
        }
        else
        {
            ranks[i] = std::max(ranks[i] - 0.1, 1.0);
        }
    }
    std::cout << std::endl;
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
    }
}

void print_state(double beta, std::array<double, INIT_NUM_STYPES> &ranks, 
    const std::array<double, INIT_NUM_STYPES> &observed_prevalence, 
    double observed_total_prevalence,
    const std::array<double, INIT_NUM_STYPES> &expected_prevalence, 
    double expected_total_prevalence,
    const std::array<double, INIT_NUM_STYPES> &prevalence_error,
    double total_prevalence_error,
    double likelihood,
    double best_likelihood,
    double weight)
{
    std::cout << "Likelihood=" << likelihood << "(best=" << best_likelihood << ")" << std::endl;
    std::cout << "Observed prevalence=" << observed_total_prevalence << "; expected prevalence=" << expected_total_prevalence << "; error=" << total_prevalence_error << "; weight = " << weight << std::endl;
    std::cout << "Contact rate (beta) = " << beta << std::endl;
    std::cout << "[" << std::endl;
    for(int z = 0; z < INIT_NUM_STYPES - 1; z++)
    {
        std::cout << "\t" << z << " : [" << ranks[z] << ", " << observed_prevalence[z] << ", " << expected_prevalence[z] << ", " << prevalence_error[z] << "]," << std::endl;
    }
    std::cout << "]" << std::endl;
}

void match_prevalence(int treatmentNumber, int simNumber, double treatment, double startingBeta)
{
    std::cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << std::endl;

    std::array<double, INIT_NUM_STYPES> betas;
    betas.fill(startingBeta);
    betas[HFLU_INDEX] = HFLU_BETA;

    std::array<double, INIT_NUM_STYPES> serotype_ranks;
    for(int i = 0; i < INIT_NUM_STYPES; i++)
    {
        serotype_ranks[i] = 1 + i;
    }

    serotype_ranks = {{1, 1, 2.3, 5.5, 7.1, 8.3, 7.9, 9.3, 10.3, 10.7, 10.9, 11.3, 11.7, 11.7, 13.1, 13.9, 14.5, 15.5, 16.7, 17.7, 19.9, 19.7, 21.9, 24.9, 25}};

    double best_likelihood = std::numeric_limits<double>::lowest();
    std::array<double, INIT_NUM_STYPES> best_betas = betas;
    std::array<double, INIT_NUM_STYPES> best_serotype_ranks = serotype_ranks;

    double total_prevalence_error = 10.0;
    double previous_total_prevalence_error = 1.0;
    double weight = INIT_WEIGHT;

    std::array<int, INIT_NUM_STYPES + 1> observed_counts = {283, 237, 184, 117, 90, 85, 84, 70, 56, 54, 53, 51, 49, 49, 43, 38, 34, 34, 29, 25, 23, 21, 19, 18, 15, 0, 1079};
    int observed_population = std::accumulate(observed_counts.begin(), observed_counts.end(), 0);
    std::array<double, INIT_NUM_STYPES>  observed_prevalence;
    double observed_total_prevalence = 0;
    for(int z = 0; z < INIT_NUM_STYPES; z++)
    {
        observed_prevalence[z] = observed_counts[z] / (double)observed_population;
        observed_total_prevalence += observed_prevalence[z];
    }

    bool fitting_beta = true;

    for(int attempt = 0; attempt < 1; attempt++)
    {
        SimPars thesePars(treatmentNumber, simNumber);
        thesePars.set_serotype_ranks(serotype_ranks);
        thesePars.set_betas(betas);
        Simulation thisSim(treatmentNumber, simNumber, &thesePars);

        thisSim.runDemSim();
        auto expected_prevalence = thisSim.runTestEpidSim();

        std::array<double, INIT_NUM_STYPES> prevalence_error;

        for(int i = 0; i < INIT_NUM_STYPES; i++)
        {
            prevalence_error[i] = expected_prevalence[i] - observed_prevalence[i];

            if(i == HFLU_INDEX)
            {
                prevalence_error[i] = 0;
            }
        }

        double likelihood = calculate_likelihood(expected_prevalence, observed_counts);

        if(likelihood > best_likelihood)
        {
            best_likelihood = likelihood;
            best_betas = betas;
            best_serotype_ranks = serotype_ranks;
        }

        double expected_total_prevalence = std::accumulate(expected_prevalence.begin(), expected_prevalence.end(), 0.0);
        total_prevalence_error = expected_total_prevalence - observed_total_prevalence;

        print_state(betas[0], serotype_ranks, observed_prevalence, 
            observed_total_prevalence, expected_prevalence, 
            expected_total_prevalence, prevalence_error, 
            total_prevalence_error, likelihood, best_likelihood, weight);

        if(abs(total_prevalence_error) < PREV_ERROR_THOLD)
        {
            fitting_beta = false;
        }
        else
        {
            fitting_beta = true;
        }

        if(fitting_beta)
        {
            std::cout << "Fitting beta" << std::endl;

            // if error changed signs and overshot, reduce weight
            if(total_prevalence_error * previous_total_prevalence_error < 0)
            {
                weight *= COOL_DOWN;
            }
            // if climbing too slowly, increase weight
            else if(abs(expected_total_prevalence - total_prevalence_error) 
                / abs(expected_total_prevalence - previous_total_prevalence_error) > TEMP_THOLD)
            {
                weight *= WARM_UP;
            }

            adjustBeta(total_prevalence_error, weight, betas);

            previous_total_prevalence_error = total_prevalence_error;
        }
        else
        {
            std::cout << "Fitting ranks" << std::endl;
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

int main(int argc, const char *argv[])
{
    int treatmentNumber = 1;
    double treatment = 1;
    int simNumber = argc > 1 ? std::stoi(argv[1]) : 1;
    std::cout << "Treatment number " << treatmentNumber << ", treatment value " << treatment << ", simulation number " << simNumber << std::endl;
    //printAssumptions();

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
    //std::cout << "Best match to treatment value " << treatment << " is beta = " << startingBeta << " (associated treatment value is " << betaTable[bestTreatment][0] << ")" << std::endl;
    adjustTreatment(treatmentNumber, treatment, simNumber);

    boost::random::mt19937_64 rng;
    rng.seed(simNumber);
    boost::random::uniform_int_distribution<int> dist;
    double beta = 1;
    for(int j = 0; j < 12; j++)
    {
        match_prevalence(treatmentNumber, dist(rng), treatment, beta);
    }
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

