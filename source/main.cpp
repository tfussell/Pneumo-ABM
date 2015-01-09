#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <thread>
#include <unordered_map>

#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Simulation.h"
#include "SimPars.h"

#define HFLU_INDEX (NUM_STYPES-1)

// Function prototypes
void adjustTreatment(int treatmentNumber, double treatment, int simNumber);
void printAssumptions();
std::string d2str(double d);
std::string makeName(int treatmentIdx, int simIdx, std::string suffix);
void printTotalTime(time_t t1, time_t t2);

double log_factorial(double n)
{
    if(n == 0)
    {
        return log(1.0);
    }

    if(n < 10)
    {
        return lgamma(n + 1);
    }

    return n * log(n) - n + (log(n * (1 + 4 * n * (1 + 2 * n)))) / 6 + log(3.141592653589793238462) / 2;
}

double calculate_likelihood(const std::array<double, NUM_STYPES> &expected, const std::array<double, NUM_STYPES + 1> &observed)
{
    double n = std::accumulate(observed.begin(), observed.end(), 0.0);
    double log_likelihood = log_factorial(n);

    double proportion_uninfected = 1 + (0.0001 * (NUM_STYPES - 1));

    for(auto a : expected)
    {
        proportion_uninfected -= a + 0.0001;
    }

    for(int i = 0; i < NUM_STYPES + 1; i++)
    {
        if(i == HFLU_INDEX)
        {
            continue;
        }

        log_likelihood -= log_factorial(observed[i]);
        log_likelihood += observed[i] * log(i == NUM_STYPES ? proportion_uninfected : (expected[i] + 0.0001));
    }

    return log_likelihood;
}

void adjust_ranks(const std::array<double, NUM_STYPES> &errors, std::array<double, NUM_STYPES> &ranks)
{
    const double delta = 1;
    //auto sum_error = std::accumulate(errors.begin(), errors.end(), 0.0, [](double a, double b) { return a + std::abs(b); });
    std::array<double, NUM_STYPES> scaled_errors;
    std::copy(errors.begin(), errors.end(), scaled_errors.begin());
    std::for_each(scaled_errors.begin(), scaled_errors.end(), [=](double &d) { d *= delta; });

    for(int i = 0; i < NUM_STYPES - 1; i++)
    {
        ranks[i] += scaled_errors[i];

        double over = 0;
        if(ranks[i] > (double)(NUM_STYPES - 1))
        {
            over = ranks[i] - (double)(NUM_STYPES - 1);
            ranks[i] = (double)(NUM_STYPES - 1);
        }
        else if(ranks[i] < 1.0)
        {
            over = 1.0 - ranks[i];
            ranks[i] = 1.0;
        }

        if(over == 0)
        {
            continue;
        }

        for(int j = 0; j < NUM_STYPES - 1; j++)
        {
            if(j == i)
            {
                continue;
            }

            ranks[j] = std::min(std::max(ranks[j] + over, 1.0), (double)(NUM_STYPES - 1));
        }
    }
}

void adjustBeta(double preve, double w, std::array<double, NUM_STYPES> &betas)
{
    for(int b = 0; b < NUM_STYPES; b++)
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

void match_prevalence(const std::array<double, NUM_STYPES> &external_ranks, int treatmentNumber, int simNumber, double treatment, double startingBeta)
{
    //std::cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << std::endl;

    std::array<double, NUM_STYPES> betas;
    betas.fill(startingBeta);
    betas[HFLU_INDEX] = HFLU_BETA;

    std::array<double, NUM_STYPES> serotype_ranks;
    for(int i = 0; i < NUM_STYPES; i++)
    {
        serotype_ranks[i] = 1 + i;
    }

   // serotype_ranks = external_ranks;

    double best_likelihood = std::numeric_limits<double>::lowest();
    std::array<double, NUM_STYPES> best_betas = betas;
    std::array<double, NUM_STYPES> best_serotype_ranks = serotype_ranks;

    double total_prevalence_error = 10.0;
    double previous_total_prevalence_error = 1.0;
    double weight = 0.8;

    std::array<double, NUM_STYPES + 1> observed_counts = {{283.0, 237.0, 184.0, 117.0, 90.0, 85.0, 84.0, 70.0, 56.0, 54.0, 53.0, 51.0, 49.0,
        49.0, 43.0, 38.0, 34.0, 34.0, 29.0, 25.0, 23.0, 21.0, 19.0, 18.0, 15.0, 13.0, 13.0, 12.0, 8.0, 7.0, 6.0, 4.0,
        4.0, 4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0}};

    observed_counts[NUM_STYPES] = 972;

    double observed_population = std::accumulate(observed_counts.begin(), observed_counts.end(), 0.0);
    std::array<double, NUM_STYPES>  observed_prevalence;
    double observed_total_prevalence = 0;
    for(int z = 0; z < NUM_STYPES; z++)
    {
        observed_prevalence[z] = observed_counts[z] / (double)observed_population;
        observed_total_prevalence += observed_prevalence[z];
    }

    bool fitting_beta = true;

	auto in_bounds = [](const std::array<double, NUM_STYPES> &expected_prevalence)
	{
		static const std::array<double, NUM_STYPES - 1> min_stop = { {
				0.088872449380254,
				0.073536416538417,
				0.056011852441207,
				0.034188101414215,
				0.025557987103319,
				0.023975128446894,
				0.023659234193544,
				0.019263454231972,
				0.014928520162317,
				0.014315434427329,
				0.014009559187710,
				0.013399206774066,
				0.012790809465225,
				0.012790809465225,
				0.010978605025572,
				0.009485614546237,
				0.008304699165959,
				0.008304699165959,
				0.006849012305781,
				0.005704588602655,
				0.005140538886973,
				0.004582886283969,
				0.004032566044146,
				0.003760514787317,
				0.002959044722980,
				0.002439490134574,
				0.002439490134574,
				0.002185148997272,
				0.001216898229325,
				0.000991529257695,
				0.000775697150860,
				0.000383884541820,
				0.000383884541820,
				0.000383884541820,
				0.000217895296624,
				0.000217895296624,
				0.000217895296624,
				0.000217895296624,
				0.000217895296624,
				0.000217895296624,
				0.000217895296624,
				0.000085296338325,
				0.000085296338325,
				0.000085296338325,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000008914681385,
				0.000000000000000,
				0.000000000000000,
				0.000000000000000,
				0.000000000000000 } };

		static const std::array<double, NUM_STYPES - 1> max_stop = { {
				0.111256066355603,
				0.094232887945256,
				0.074474850534474,
				0.049170387167968,
				0.038810363450686,
				0.036876749335286,
				0.036489359861464,
				0.031039632574562,
				0.025530276287365,
				0.024737177717045,
				0.024339977652025,
				0.023544215824669,
				0.022746551379733,
				0.022746551379733,
				0.020340943746101,
				0.018319633758304,
				0.016689587309438,
				0.016689587309438,
				0.014632395091039,
				0.012967439515854,
				0.012127222488374,
				0.011280974318645,
				0.010427850017154,
				0.009998387337851,
				0.008696404901838,
				0.007814884055676,
				0.007814884055676,
				0.007369195582954,
				0.005542834682162,
				0.005071772367006,
				0.004592693147404,
				0.003602232109533,
				0.003602232109533,
				0.003602232109533,
				0.003083934870743,
				0.003083934870743,
				0.003083934870743,
				0.003083934870743,
				0.003083934870743,
				0.003083934870743,
				0.003083934870743,
				0.002541565200340,
				0.002541565200340,
				0.002541565200340,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001960267967797,
				0.001298058009173,
				0.001298058009173,
				0.001298058009173,
				0.001298058009173 } };

		for (int z = 0; z < NUM_STYPES - 1; z++)
		{
			if (expected_prevalence[z] < min_stop[z] || expected_prevalence[z] > max_stop[z])
			{
				return false;
			}
		}

		return true;
	};

	while (true)
	{
        SimPars thesePars(treatmentNumber, simNumber);
        thesePars.set_serotype_ranks(serotype_ranks);
        thesePars.set_betas(betas);
        Simulation thisSim(treatmentNumber, simNumber, &thesePars);

        thisSim.runDemSim();
        auto expected_prevalence = thisSim.runTestEpidSim();

        std::array<double, NUM_STYPES> prevalence_error;

        for(int i = 0; i < NUM_STYPES; i++)
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

        double sum_absolute_errors = std::accumulate(prevalence_error.begin(), prevalence_error.end(), 0.0, [](double t, double d) { return t + std::abs(d); });

        std::cout << "Sum of absolute errors=" << sum_absolute_errors << std::endl;
        std::cout << likelihood << "(best=" << best_likelihood << ")" << std::endl;
        std::cout << "Observed prevalence=" << observed_total_prevalence << "; expected prevalence=" << expected_total_prevalence << "; error=" << total_prevalence_error << "; weight = " << weight << std::endl;
        std::cout << "Contact rate (beta) = " << betas[0] << std::endl;
        std::cout << "[" << std::endl;
        for(int z = 0; z < NUM_STYPES - 1; z++)
        {
            std::cout << "\t" << SerotypeNames[z] << "(" << z << ")" << " : [" << serotype_ranks[z] << ", " << observed_prevalence[z] << ", " << expected_prevalence[z] << ", " << prevalence_error[z] << "]," << std::endl;
        }
        std::cout << "]" << std::endl;

        if(abs(total_prevalence_error) < PREV_ERROR_THOLD)
        {
            fitting_beta = false;
        }
        else
        {
            fitting_beta = true;
        }

		if (in_bounds(expected_prevalence)) break;

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
    //std::cout << "\tBeginning simulation #" << simNumber << std::endl;
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
    double treatment = 0.3;
    int simNumber = 1;
    
    std::array<double, NUM_STYPES> external_ranks;

    for(int i = 1; i <= std::min(NUM_STYPES, argc - 1); i++)
    {
        external_ranks[i - 1] = std::stod(argv[i]);
    }

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
    adjustTreatment(treatmentNumber, treatment, simNumber);

    double beta = 0.0583273/*0.0580585*/;
    match_prevalence(external_ranks, treatmentNumber, simNumber, treatment, beta);
}

void adjustTreatment(int treatmentNumber, double treatment, int simNumber) {
    double thisXI[NUM_STYPES][NUM_STYPES];
    for(int i = 0; i < NUM_STYPES; i++) {
        for(int j = 0; j < NUM_STYPES; j++) {
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
    for(int i = 0; i < NUM_STYPES; i++) {
        for(int j = 0; j < NUM_STYPES; j++) {
            xiStream << thisXI[i][j];
            if((j + 1) % NUM_STYPES == 0) {
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

