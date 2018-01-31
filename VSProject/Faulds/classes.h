#pragma once

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <functional>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <cmath>
#include <Windows.h> //just needed for multi-threading

//features enabled
#define WIN_ENABLED 1
#define TTL_ENABLED 1
#define DF_ENABLED 1
#define MSS_ENABLED 1
#define OPT_ENABLED 1
#define RST_ENABLED 1

//EM Iterations
#define EM_MAX_ITERATIONS 100

//Hershel+ parameters
#define HP_JITTER_LOSS_THRESHOLD 2
#define HP_OWD_MEAN 0.5
#define HP_OWD_LAMBDA  (1.0 / HP_OWD_MEAN)
#define HP_ERLANG_MEAN 0.5
#define HP_ERLANG_NU 4
#define HP_ERLANG_K 2
#define HP_LOSS_PROB 0.038
#define HP_FEATURE_PROB_PI_RST_OPT 0.01
#define HP_FEATURE_PROB_PI 0.1

//Network delay PMF sizes
#define T_GUESS_LIMIT_IN_SECS 1.0
#define T_NUM_BINS 100
#define OWD_GUESS_LIMIT_IN_SECS 1.0
#define OWD_NUM_BINS 1000
//#define T_STEP_SIZE (T_GUESS_LIMIT_IN_SECS / T_NUM_BINS)
//#define OWD_STEP_SIZE (Q_GUESS_LIMIT_IN_SECS / Q_NUM_BINS)
#define T_STEP_SIZE 0.01
#define OWD_STEP_SIZE 0.001

#define OPT_MSS 1 // max segment size
#define OPT_WIN 2 // window scaling
#define OPT_TS 3 // timestamps
#define OPT_SACK 4 // SACK allowed
#define OPT_EOL 5 // end of list (another type of padding?)
#define OPT_NOP 6 // padding

#define COMPACTION_THRESHOLD 0.001

#ifdef _DEBUG
#define NUM_THREADS 1
#else 
#define NUM_THREADS GetActiveProcessorCount(ALL_PROCESSOR_GROUPS)
#endif

using namespace std;

class Signature{
public:
	int subOS;
	unsigned int ip_int;
	unsigned int predicted_ip_int;
	int predicted_winner_index;
	bool was_correct;
	int win;
	int ttl;
	int df;
	char options[25];
	unsigned long long options_int;
	int option_order;
	int mss;
	int rst;
	int rst_ack;
	int rst_win;
	int rst_seq;
	int rst_nonzero;
	int rst_int;
	vector<double> packet_arrival_timestamp;
	double prob;

	void print(){
		printf("ClassID: %d, PredictedClass: %d\n", ip_int, predicted_ip_int);
		for (double d : packet_arrival_timestamp) printf(" %lf ", d);
		printf("\n");
	}
};

static double factorial(int n){
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

static void normalize(vector<double>& v){
	double sum = 0;
	for (int i = 0; i < v.size(); i++) sum += v[i];
	if (sum > 0) for (int i = 0; i < v.size(); i++) v[i] /= sum;
}

static void getRSTbits(int rst_int, int& rst, int& ra, int& rw, int& rn) {
	rn = rst_int & 1;
	rst_int >>= 1;
	rw = rst_int & 3;
	rst_int >>= 2;
	ra = rst_int & 1;
	rst_int >>= 1;
	rst = rst_int & 1;
}

static void getOptionString(int option_int, char* str, int size) {
	if (option_int == 0) {
		strcpy_s(str, size, "X");
		return;
	}

	//add options to string
	strcpy_s(str, size, "\0");
	int count = 0;
	while (option_int > 0) {
		int opt = option_int & 0x7; //AND by 3 bits to get option
		if (opt == OPT_MSS) strcat_s(str, size, "M");
		else if (opt == OPT_WIN) strcat_s(str, size, "W");
		else if (opt == OPT_SACK) strcat_s(str, size, "S");
		else if (opt == OPT_TS) strcat_s(str, size, "T");
		else if (opt == OPT_NOP) strcat_s(str, size, "N");
		else if (opt == OPT_EOL) strcat_s(str, size, "E");

		option_int = option_int >> 3;
		count++;
	}

	//reverse the string to get the correct order	
	for (int i = 0, j = count - 1; i < j; i++, j--) {
		swap(str[i], str[j]);
	}

	//add a null since we are expert programmers who never forget this...
	strcat_s(str, size, "\0");
}

static char ipbuffer[18];
static char* my_ntoa(UINT ip) {
	unsigned char *bytes = (unsigned char *)&ip;
	_snprintf_s(ipbuffer, sizeof(ipbuffer), "%d.%d.%d.%d", bytes[0], bytes[1], bytes[2], bytes[3]);

	return ipbuffer;
}

//extern unordered_map<int, int> win_map;
//extern unordered_map<int, int> mss_map;
//extern unordered_map<int, int> rst_map;
extern unordered_map<int, int> ttl_map;
extern unordered_map<long long int, int> opt_map;
extern unordered_map<int, unordered_set<int>> option_subsets;
extern unordered_map<unsigned int, string> IP_to_OS_map;


class HershelPlusEM{
public:
	//current pmfs
	static vector<double>* t_pmf;
	static vector<double>* owd_pmf;
	static vector<vector<double>>* loss_pmfs;
	static vector<vector<double>>* n_choose_k;
	static vector<vector<double>>* windowsize_pmf;
	static vector<vector<double>>* ttl_pmf;
	static vector<vector<double>>* df_pmf;
	static vector<vector<double>>* mss_pmf;
	static vector<vector<double>>* opt_pmf;
	static vector<vector<double>>* rst_pmf;

	//new pmfs to be built
	vector<double>* new_t_pmf;
	vector<double>* new_owd_pmf;

	double totalprob;
	int y_nPkts;
	int x_nPkts;
	double subOS_prob;
	double prob_so_far;
	vector<int> accumulator;
	vector<double>* x_timestamps;
	vector<double>* y_timestamps;
	vector<int> owd_indexes;	

	//debug vars
	Signature& sample_sig;
	int dbnum;
	FILE* fdebug;

	HershelPlusEM(Signature& x, vector<Signature>& class_sigs, int osnum, double alpha_j, vector<double>* new_t_pmf, vector<double>* new_owd_pmf, vector<bool>* opt_flags);

private:
	void ExamineCombination(vector<int>& gamma);
	void ProduceLossPatterns(int remaining, int start, int accumulatorSize);
	int options_match_int64(unsigned long long tar_opts, unsigned long long sig_opts);
	double MatchConstantFeatures(Signature& target, Signature& dbsig, int osnum, vector<bool>* opt_flags);
};

class HershelPlus {
public:
	//current pmfs
	static vector<double> owd_pmf;
	static vector<double> t_pmf;

	double totalprob;
	int y_nPkts;
	int x_nPkts;
	double subOS_prob;
	vector<int> accumulator;
	vector<double> x_timestamps;
	vector<double> y_timestamps;

	HershelPlus(Signature& observation, vector<Signature>& class_sigs);

private:
	void ExamineCombination(vector<int>& gamma);
	void ExamineCombinationOptimized(vector<int>& gamma);
	void ProduceLossPatterns(int remaining, int start, int accumulatorSize);
	int options_match_int64(unsigned long long tar_opts, unsigned long long sig_opts);
	double MatchConstantFeatures(Signature& target, Signature& dbsig);
};

void Faulds(map<u_int, vector<Signature>> original_sigs, vector<Signature> disturbed, int max_num_packets, bool read_iteration_data);
void HardEM_MT(map<u_int, vector<Signature>> original_sigs, vector<Signature> disturbed);
