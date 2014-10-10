#include <cstdlib>

#include "SimPars.h"
#include "Parameters.h"

#define HFLU_INDEX (NUM_STYPES-1)

SimPars::SimPars(int tID, int sID) {
    treatmentID = tID;
    simID = sID;
    initializeDemInput();
    initializeEpidInput();

}

SimPars::~SimPars() {
}

// PUBLIC MEMBER FUNCTIONS
const age_index& SimPars::get_demPMF_row(int row) const {
    return demPMFs[row];
}

const double SimPars::get_serotypePar_ij(int i, int j) const {
    return serotypePars[i][j];
}

const double SimPars::get_XI_ij(int i, int j) const {
    return XI[i][j];
}

const double SimPars::get_waifw_ij(int i, int j) const {
    return normWAIFW[i][j];
}

const kids_index& SimPars::get_parity() const {
    return parity_pmf;
}

const double SimPars::get_neighbor(int i, int j) const {
    return neighborhoods[i][j];
}

const double SimPars::get_normalized_neighbor(int i, int j) const {
    return NNN[i][j];
}

const double SimPars::get_Hflu_prob(int s) const {
    return Hflu_probs[s];
}

const double SimPars::get_reductions(int s) const {
    return reductions[s];
}



// PRIVATE MEMBER FUNCTIONS
void SimPars::initializeDemInput() {

    // First initialize all matrices to zero
    for(int i = 0; i < NUM_SOCIODEM_FILES; i++) {
        for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
            demPMFs[i][j] = 0.0;
        }
    }
    for(int i = 0; i < INIT_NUM_AGE_CATS; i++) {
        for(int j = 0; j < INIT_NUM_AGE_CATS; j++) {
            WAIFW[i][j] = 0.0;
            normWAIFW[i][j] = 0.0;
        }
    }
    for(int i = 0; i < PARITY_BUFFER; i++) {
        parity_pmf[i] = 0.0;
    }
    for(int n1 = 0; n1 < NUM_NEIGHBORHOODS; n1++) {
        for(int n2 = 0; n2 < NUM_NEIGHBORHOODS; n2++) {
            neighborhoods[n1][n2] = 0;
            NNN[n1][n2] = 0;
        }
    }

    // Sociodemographic data
    for(int i = 0; i < NUM_SOCIODEM_FILES; i++) {
        std::ifstream thisFile;
        thisFile.open("../../inputs/" + std::string(SOCIODEM_FILENAMES[i]), std::ios::in);
        if(!thisFile) {
            std::cerr << "Error reading " << SOCIODEM_FILENAMES[i] << " file." << std::endl;
            throw std::runtime_error("");
        }
        double currentMax = 0.0;
        int maxIndex = 0;
        double cumulativeSum = 0.0;
        double thisVal;
        int ii = 0;
        while(!thisFile.eof()) {
            thisFile >> thisVal;
            if(ii < INIT_NUM_AGE_CATS) {
                demPMFs[i][ii] = thisVal;
                cumulativeSum += thisVal;
                if(thisVal > currentMax) {
                    currentMax = thisVal;
                    maxIndex = ii;
                }
                ii++;
            }
        }
        if(cumulativeSum != 1.0) {
            //std::cout << "Cumulative sum of " << SOCIODEM_FILENAMES[i] << " is " << cumulativeSum << ", not 1. ";
            double diff = 1.0 - cumulativeSum;
            if(abs(diff) < ERR_EPSILON) {
                if(diff > 0) {
                    //std::cout << "Adding difference (" << diff << ") to highest rate." << std::endl;
                    demPMFs[i][maxIndex] += diff;
                }
                else {
                    //std::cout << "Adding difference (" << diff << ") to highest rate." << std::endl;
                    demPMFs[i][maxIndex] += diff;
                }
            }
            else {
                //std::cout << "Difference (" << diff << ") is too great (max acceptable is " << ERR_EPSILON << "). Fix the input file." << std::endl;
                throw std::runtime_error("");
            }
        }
    }


    // WAIFW
    std::ifstream thisFile3;
    thisFile3.open("../../inputs/WAIFW.txt", std::ios::in);
    if(!thisFile3) {
        std::cerr << "Error reading " << WAIFW_FILENAME << "." << std::endl;
        throw std::runtime_error("");
    }
    double thisVal;
    int ii = 0;
    int i = 0;
    int j = 0;
    while(!thisFile3.eof()) {
        thisFile3 >> thisVal;
        if(ii < INIT_NUM_AGE_CATS * INIT_NUM_AGE_CATS) {
            if(thisVal == 0) {
                thisVal = WAIFW_NONZERO;
            }
            WAIFW[i][j] = thisVal;
            j++;
            if(j == INIT_NUM_AGE_CATS) {
                j = 0;
                i++;
            }
            ii++;
        }
    }

    // Normalize WAIFW
    for(int j = 0; j < INIT_NUM_AGE_CATS; j++) { // for each column (FROM)
        double colSum = 0.0;
        for(int i = 0; i < INIT_NUM_AGE_CATS; i++) {
            colSum += WAIFW[i][j];
        }
        for(int i = 0; i < INIT_NUM_AGE_CATS; i++) {
            normWAIFW[i][j] = WAIFW[i][j] / colSum;
        }
    }

    // PARITY_PMF
    std::ifstream thisFile2;
    thisFile2.open("../../inputs/" + std::string(PARITY_FILENAME), std::ios::in);
    if(!thisFile2) {
        std::cerr << "Error reading " << PARITY_FILENAME << "." << std::endl;
        throw std::runtime_error("");
    }
    thisVal = 0;
    i = 0;
    double currentMax = 0.0;
    int maxIndex = 0;
    double cumulativeSum = 0.0;
    while(!thisFile2.eof()) {
        thisFile2 >> thisVal;
        if(i < PARITY_BUFFER) {
            parity_pmf[i] = thisVal;
            cumulativeSum += thisVal;
            if(thisVal > currentMax) {
                currentMax = thisVal;
                maxIndex = i;
            }
            i++;
        }
    }
    double diff = 1.0 - cumulativeSum;
    if(abs(diff) >= 0 && abs(diff) < ERR_EPSILON) {
        //std::cout << "Cumulative sum of " << PARITY_FILENAME << " is " << cumulativeSum << ", not 1." << std::endl;
        if(diff < 0) {
            //std::cout << "Adding difference (" << diff << ") to highest rate." << std::endl;
            parity_pmf[maxIndex] += diff;
        }
        else {
            //std::cout << "Adding difference (" << diff << ") to highest rate." << std::endl;
            parity_pmf[maxIndex] += diff;
        }
    }
    else {
        //std::cerr << "Difference (" << diff << ") exceeds acceptable threshold (" << ERR_EPSILON << ")." << std::endl;
        throw std::runtime_error("");
    }

    // NEIGHBORHOODS
    std::ifstream thisFile4;
    thisFile4.open("../../inputs/" + std::string(NEIGHBORHOOD_FILENAME), std::ios::in);
    if(!thisFile4) {
        std::cerr << "Error reading " << NEIGHBORHOOD_FILENAME << "." << std::endl;
        throw std::runtime_error("");
    }
    double thisVal2;
    int ii2 = 0;
    int i2 = 0;
    int j2 = 0;
    while(!thisFile4.eof() && ii2 < NUM_NEIGHBORHOODS * NUM_NEIGHBORHOODS) {
        thisFile4 >> thisVal2;
        if(ii2 < NUM_NEIGHBORHOODS * NUM_NEIGHBORHOODS) {
            neighborhoods[i2][j2] = thisVal2;
            j2++;
            if(j2 == NUM_NEIGHBORHOODS) {
                j2 = 0;
                i2++;
            }
            ii2++;
        }
    }

    // Normalize neighborhoods
    for(int j = 0; j < NUM_NEIGHBORHOODS; j++) { // for each column (FROM)
        double colSum = 0.0;
        for(int i = 0; i < NUM_NEIGHBORHOODS; i++) {
            colSum += neighborhoods[i][j];
        }
        for(int i = 0; i < NUM_NEIGHBORHOODS; i++) {
            NNN[i][j] = neighborhoods[i][j] / colSum;
        }
    }

}

void SimPars::set_betas(const std::array<double, NUM_STYPES> &betas)
{
    for(int s = 0; s < NUM_STYPES; s++) {
        serotypePars[BETA_INDEX][s] = betas[s];
    }
}

void SimPars::set_serotype_ranks(const std::array<double, NUM_STYPES> &ranks)
{
    double maxDuration = 220.0;
    double minDuration = BASE_DURATION;
    double rangeDuration = maxDuration - minDuration;

    for(int s = 0; s < NUM_STYPES; s++) {
        if(MAX_REDUCTION > 0) {
            if(s < HFLU_INDEX) {
                double x = (ranks[s] - 1) / (NUM_STYPES - 2);
                serotypePars[MEAN_DURATION_INDEX][s] = maxDuration - x * rangeDuration;
                reductions[s] = MAX_REDUCTION * (1 - x);
            }
            else {
                serotypePars[MEAN_DURATION_INDEX][s] = 1;
                reductions[s] = 0.0;
            }
        }
        else {
            reductions[s] = 0.0;
        }
    }
}

void SimPars::initializeEpidInput() {

    // Epidemiological data
    // Read in files
    for(int f = 0; f < NUM_EPID_FILES; f++) {
        std::ifstream thisFile;
        thisFile.open("../../inputs/" + std::string(EPID_FILENAMES[f]), std::ios::in);
        if(!thisFile) {
            std::cerr << "Error reading " << EPID_FILENAMES[f] << " file." << std::endl;
            throw std::runtime_error("");
        }
        double thisVal;
        int ii = 0;
        while(!thisFile.eof()) {
            thisFile >> thisVal;
            if(ii < NUM_STYPES) {
                serotypePars[f][ii] = thisVal;
                ii++;
            }
        }
    }

    // XI
    std::ifstream thisFile;
    std::string XIFile = makeName(treatmentID, 1, "XI");
    thisFile.open("../../outputs/" + XIFile, std::ios::in);
    if(!thisFile) {
        std::cerr << "Error reading " << XIFile << "." << std::endl;
        throw std::runtime_error("");
    }
    double thisVal;
    int ii = 0;
    int i = 0;
    int j = 0;
    while(!thisFile.eof()) {
        thisFile >> thisVal;
        if(ii < NUM_STYPES * NUM_STYPES) {
            XI[i][j] = thisVal;
            j++;
            if(j == NUM_STYPES) {
                j = 0;
                i++;
            }
            ii++;
        }
    }

    // HFLU_PROBS
    std::ifstream hfStream;
    hfStream.open("../../inputs/" + std::string(HFPROB_FILENAME), std::ios::in);
    if(!hfStream) {
        std::cerr << "Error reading " << HFPROB_FILENAME << "." << std::endl;
        throw std::runtime_error("");
    }
    i = 0;
    while(!hfStream.eof() && (i < NUM_STYPES - 1)) {
        hfStream >> thisVal;
        Hflu_probs[i] = thisVal;
        i++;
    }

    // Set up competition terms
    double maxDuration = serotypePars[MEAN_DURATION_INDEX][0];
    double minDuration = BASE_DURATION;
    double rangeDuration = maxDuration - minDuration;
    double thisDuration;
    for(int s = 0; s < NUM_STYPES; s++) {
        if(MAX_REDUCTION > 0) {
            if(s < HFLU_INDEX) {
                thisDuration = serotypePars[MEAN_DURATION_INDEX][s];
                reductions[s] = MAX_REDUCTION*(thisDuration - minDuration) / rangeDuration;
            }
            else {
                reductions[s] = 0.0;
            }
        }
        else {
            reductions[s] = 0.0;
        }
    }
}

std::string SimPars::d2str(double d) {
    std::stringstream t;
    t << d;
    return t.str();
}

std::string SimPars::makeName(int treatmentIdx, int simIdx, std::string suffix) {
    std::string thisCtr = d2str(simIdx);
    std::string thisTr = d2str(treatmentIdx);
    std::string thisName = "tr_" + thisTr + "_sim_" + thisCtr + "_" + suffix;
    return thisName;
}
