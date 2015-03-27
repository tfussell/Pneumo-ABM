#pragma once

#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

// All temporal units are days
// Haemophilus influenzae is the last serotype with index NUM_STYPES-1

// MODEL OPTIONS
#define MATCH_PREVALENCE // turn off to run with transmission rates in Betas_used.txt and Treatments.txt
#define NO_HHOLDS // if defined, contact rates are independent of household status
#define NO_AGE_ASSORT // if defined, contact rates are independent of host age
#define SIM_PCV // if on, introduces vaccine

// SIMULATION PARAMETERS
// ...input:
const double DEM_SIM_LENGTH = (double)300 * 365.0; // period of demographic burn-in (no epid dynamics); should be a multiple of strobing interval
const double EPID_SIM_LENGTH = (double)150.0*365.0; // duration of epid dynamics; total sim time is EPID_SIM_LENGTH + DEM_SIM_LENGTH
const double EPID_DELTA_T = (double)1.0; // time step for calculating force of colonization (in days); keep low (~1 day) to avoid error
const double APPROX_NOW = -10.0; // pow(10,APPROX_NOW) used for adjusting times to prevent event collisions

// ...fitting of transmission rate (beta) (parameters below are ignored unless MATCH_PREVALENCE defined)
const double PREV_ERROR_THOLD = 0.01; // allowed error to fit target prevalence (0.01 = target prevalence +/- 1%) 
const int NUM_TEST_SAMPLES = 20; // number of epid strobes used to determine if target criterion has been met (prevalence is average of samples) 
const double COOL_DOWN = 0.8; // decrease in jump size if vacillating around target prevalence prevalence
const double WARM_UP = 1.1; // increase in jump size if undershooting target prevalence
const double TEMP_THOLD = 0.7; // if prevalence is less that TEMP_THOLD fraction of target, jump size (temperature of search) increases (WARM_UP)
const double TEST_EPID_SIM_LENGTH = (double)150 * 365.0; // duration of epid dynamics for simulations to fit prevalence
const double INIT_WEIGHT = 0.8; // initial coefficient for jump size (jump size = INIT_WEIGHT * difference in prevalence)

// ...output:
const bool CHECK_INPUT = true; // performs non-comprehensive logic checks on input
const double STROBE_EPID = (double)365.0; // interval (in time) between epid observations - Matlab post-processing assumes 365
const double STROBE_DEM = (double)365.0; // interval (in time) between demographic observations - see above
const double PROGRESS_INTERVAL = 10.0; // % interval at which to report progress to screen for each component
const int COCOL_AGE_LIMIT = 5; // in *years*; Haemophilus-pneumo and pneumo-pneumo co-colonization stats printed for hosts <COCOL_AGE_LIMIT

// SOCIODEMOGRAPHIC PARAMETERS
const int N0 = 150000; // initial population size
const double MATURITY_AGE = (double)15.0;
const int TSTEPS_AGE = 365; // EPID_DELTA_T per age
const int INIT_NUM_AGE_CATS = 111; // if NO_AGE_ASSORT *not* defined, number of age categories (assume categories are YEARS), older ages borrow rates from NUM_AGES
const char * const WAIFW_FILENAME = {"WAIFW.txt"}; // unless NO_AGE_ASSORT defined, age-associated contact matrix (will be normalized)
const int NUM_NEIGHBORHOODS = 1; // do not change--multiple neighborhoods deprecated
const int NUM_SOCIODEM_FILES = 5; // used below
const char * const SOCIODEM_FILENAMES[NUM_SOCIODEM_FILES] = {
    "LSPAN_PMF.txt", // probability mass function of lifespans, indexed by age (in y)
    "FLEDGE_PMF.txt",  // pmf of leaving home of origin (if not already paired), indexed by age (in y)
    "PAIR_PMF.txt",  // pmf of initiating a pairing, by age (in y)
    "BIRTH_AGE_PMF.txt",  // pmf of age at reproduction (in y)
    "INIT_AGE_PMF.txt" // initial age distribution (in y)
};
enum { LSPAN_INDEX, FLEDGE_INDEX, PAIR_INDEX, BIRTH_AGE_INDEX, INIT_AGE_INDEX }; // MUST match file order above
const double INIT_FRAC_PAIRED = 0.40; // fraction of mature adults that are partnered 
const double PROB_PAIR = 0.9; // fraction of the population that should be paired at some point in life
const double STD_AGE_PAIR = 3.0; // half of acceptable range in ages between partners
const char * const PARITY_FILENAME = {"PARITY_PMF.txt"}; // probability of having [i] kids
const char * const NEIGHBORHOOD_FILENAME = {"NEIGHBORHOODS.txt"}; // adjacency matrix for neighborhoods
const char * const HFPROB_FILENAME = {"HFLU_PROBS.txt"}; // probability of serotype being immediately cleared in presence of Hflu
const int HHOLD_SIZE_BUFFER = 30; // buffer for maximum household size
const int PARITY_BUFFER = 50; // max number of offspring permitted
const int DATE_BUFFER = 100; // max number of searches for partner w/in ideal age range; can usually ignore
const double ERR_EPSILON = 0.00007; // acceptable remainder in sum of PMFs; if less than this amount, will add to highest rate class

// EPIDEMIOLOGICAL PARAMETERS
// ...initialization
const int NUM_STYPES = 43 + 1; // = initial number of pneumo serotypes + Haemophilus influenzae
const int NUM_EPID_FILES = 3;
const char * const EPID_FILENAMES[NUM_EPID_FILES] = {
    "INIT_INFECTEDS.txt", // initial fraction of population colonized with each serotype and H. influenzae
    "IMMIG_RATES.txt", // external immigration rate for each serotype and H. influenzae
    "M_DURATION_INFECTION.txt", // mean strain-specific intrinsic duration of infection
};
enum { INIT_INFECTEDS_INDEX, IMMIGRATION_INDEX, MEAN_DURATION_INDEX, BETA_INDEX }; // MUST match file order above		

// ...transmission
const double BASE_DURATION = 25.0; // minimum mean intrinsic duration of carriage
const int COINFECTION_BUFFER = 100; // maximum possible number of coinfections
const double WAIFW_NONZERO = 0.000001; // unless NO_AGE_ASSORT defined, replacement for WAIFW entry ij = 0 (allows w/in household contacts)
const double RHO_H = 0.4; // unless NO_HHOLDS defined, fraction of transmission within household
const double HFLU_BETA = 0; // transmission rate of H. influenzae

// ...immunity
const char * const XI_FILENAME = {"XI.txt"}; // if VARY_SSI not defined, strain-specific cross-immunity (from col j to row i)
const double REC_EPS = 0.25; // coefficient (epsilon) for reduction in duration from current & past carriage
const double MAX_REDUCTION = 0.25; // FOI reduced by (1-MAX_REDUCTION) for res. stype 0; interpolated linearly to zero for others; if 0, no competition from resident
const double HFLU_SIGMA = 0.3; // reduction in susceptibility to H. flu if host has carried it before
const double RSCC_HFLU = 1.0; // reduction in susceptibility to H. flu if carrying H. flu (ensures single-strain dynamics w/o co-colonization)
const double EPSILON = 0.0001; // future time (in days) to 'instantaneous' recovery 

// ...vaccination (ignored unless SIM_PCV defined)
const double INIT_VACCINE_EFFICACY = 0.73; // percent reduction in susceptibility to serotypes in vaccine (susceptibility is max(1-vaccine efficacy,XI))
const double VACCINE_AGE = 30 * 6; // in days--host age at which vaccinated
const double VACCINE_COVERAGE_AGE = 30 * 13; // in days--host age at which coverage is evaluated--this is different from VACCINE_AGE because it's based off real data which isn't necessarily collected for patients at VACCINE_AGE

const std::vector<int> SampleTimes =
{
	1981,
	1982,
	1983,
	1984,
	1985,
	1986,
	1987,
	1988,
	1989,
	1990,
	1991,
	1992,
	1993,
	1994,
	1995,
	1996,
	1997,
	1998,
	1999,
	2000,
    2001,
    2002,
    2003,
	2004,
	2005,
	2007,
	2009,
	2010,
	2011,
	2012,
	2014,
	2016,
	2018,
	2020,
	2050
};

const std::unordered_map<std::string, std::vector<std::string>> VaccineTypes = 
{
    {"PCV7", {"4", "6A", "6B", "9V", "14", "18C", "19F", "23F"}},
    {"PCV10", {"4", "6B", "9V", "14", "18C", "19F", "23F", "1", "5", "7F"}},
    {"PCV13", {"4", "6B", "9V", "14", "18C", "19F", "23F", "1", "5", "7F", "3", "6A", "19A"}}
};

const std::array<std::string, NUM_STYPES> SerotypeNames = 
{
	"6A",
	"23F",
	"19F",
	"6B",
	"11A",
	"15B/C",
	"19A",
	"35A/B",
	"14",
	"22F",
	"9A",
	"18C",
	"NT",
	"10",
	"6C",
	"9N",
	"23A",
	"35F",
	"23B",
	"3",
	"34",
	"4",
	"31",
	"15A",
	"38",
	"15F",
	"29",
	"25A",
	"7F",
	"16F",
	"33F",
	"Pool I",
	"17F",
	"21",
	"37",
	"9V",
	"7C",
	"33A",
	"13",
	"18F",
	"36",
	"20",
	"24F",
    "Flu"
};

const std::unordered_map<int, std::pair<double, std::string>> YearlyVaccinationCoverage =
{
	{ 2001, { 0, "PCV7" } },
	{ 2002, { 0.451, "PCV7" } },
	{ 2003, { 0.864, "PCV7" } },
	{ 2004, { 0.860, "PCV7" } },
	{ 2005, { 0.873, "PCV7" } },
	{ 2006, { 0.904, "PCV7" } },
	{ 2007, { 0.949, "PCV7" } },
	{ 2008, { 0.940, "PCV7" } },
	{ 2009, { 0.928, "PCV7" } },
	{ 2010, { 0.941, "PCV13" } },
	{ 2011, { 0.964, "PCV13" } },
	{ 2012, { 0.895, "PCV13" } },
	{ 2013, { 0.956, "PCV13" } },
	{ 2014, { 0.956, "PCV13" } },
	{ 2015, { 0.956, "PCV13" } },
	{ 2016, { 0.956, "PCV13" } },
	{ 2017, { 0.956, "PCV13" } },
	{ 2018, { 0.956, "PCV13" } },
	{ 2019, { 0.956, "PCV13" } },
	{ 2020, { 0.956, "PCV13" } },
	{ 2021, { 0.956, "PCV13" } },
	{ 2022, { 0.956, "PCV13" } },
	{ 2023, { 0.956, "PCV13" } },
	{ 2024, { 0.956, "PCV13" } },
	{ 2025, { 0.956, "PCV13" } },
	{ 2026, { 0.956, "PCV13" } },
	{ 2027, { 0.956, "PCV13" } },
	{ 2028, { 0.956, "PCV13" } },
	{ 2029, { 0.956, "PCV13" } },
	{ 2030, { 0.956, "PCV13" } },
	{ 2031, { 0.956, "PCV13" } },
	{ 2032, { 0.956, "PCV13" } },
	{ 2033, { 0.956, "PCV13" } },
	{ 2034, { 0.956, "PCV13" } },
	{ 2035, { 0.956, "PCV13" } },
	{ 2036, { 0.956, "PCV13" } },
	{ 2037, { 0.956, "PCV13" } },
	{ 2038, { 0.956, "PCV13" } },
	{ 2039, { 0.956, "PCV13" } },
	{ 2040, { 0.956, "PCV13" } },
	{ 2041, { 0.956, "PCV13" } },
	{ 2042, { 0.956, "PCV13" } },
	{ 2043, { 0.956, "PCV13" } },
	{ 2044, { 0.956, "PCV13" } },
	{ 2045, { 0.956, "PCV13" } },
	{ 2046, { 0.956, "PCV13" } },
	{ 2047, { 0.956, "PCV13" } },
	{ 2048, { 0.956, "PCV13" } },
	{ 2049, { 0.956, "PCV13" } },
	{ 2050, { 0.956, "PCV13" } },
	{ 2051, { 0.956, "PCV13" } }
};
