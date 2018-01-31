/*
Implementation of Hershel+ for use with Faulds. 
Stores static vectors for quick lookup/building of feature PMFs
Author: Zain Shamsi
*/


#include "classes.h"

//will be set by Faulds
vector<double>* HershelPlusEM::t_pmf = NULL;
vector<double>* HershelPlusEM::owd_pmf = NULL;
vector<vector<double>>* HershelPlusEM::loss_pmfs = NULL;
vector<vector<double>>* HershelPlusEM::n_choose_k = NULL;
vector<vector<double>>* HershelPlusEM::windowsize_pmf = NULL;
vector<vector<double>>* HershelPlusEM::ttl_pmf = NULL;
vector<vector<double>>* HershelPlusEM::df_pmf = NULL;
vector<vector<double>>* HershelPlusEM::mss_pmf = NULL;
vector<vector<double>>* HershelPlusEM::opt_pmf = NULL;
vector<vector<double>>* HershelPlusEM::rst_pmf = NULL;

void HershelPlusEM::ExamineCombination(vector<int>& gamma){
	double gamma_prob = 0;
	double t_limit = T_GUESS_LIMIT_IN_SECS;

	//but T has to be less than the smallest distance between any x and y
	for (int m = 0; m < x_nPkts; m++){
		double d = (*x_timestamps)[m] - (*y_timestamps)[gamma[m]];
		if (d < t_limit) t_limit = d;
	}

	//for (double t_guess = 0; t_guess <= t_limit; t_guess += T_STEP_SIZE){
	for (int t_index = 0; t_index < T_NUM_BINS; t_index++){
		double probOWD = 1;

		//int t_index = nearbyint(t_guess / T_STEP_SIZE); //t_guess / T_STEP_SIZE; 
		double t_guess = T_STEP_SIZE * t_index;
		if (t_guess > t_limit) break;

		double probT = (*t_pmf)[t_index];

		//for each received packet, determine the OWD and calculate \prod p(OWD)		
		boolean skip_update = false;
		for (int m = 0; m < x_nPkts; m++){
			double OWD = ((*x_timestamps)[m] - (*y_timestamps)[gamma[m]] - t_guess);
			if (OWD < 0){
				probOWD = 0;
				break;
			}
			if (t_guess < OWD) { //heuristic to skip update if T < OWD
				skip_update = true;
			}
			else {
				int owd_bin = OWD / OWD_STEP_SIZE;
				double p = 0;
				if (owd_bin > OWD_NUM_BINS - 1){ 					
					//if past the range, set bin to last and let p=0 to not count this case
					owd_bin = OWD_NUM_BINS - 1;
				}
				else p = (*owd_pmf)[owd_bin];
				probOWD *= p;

				//store index to build new pmf
				owd_indexes[m] = owd_bin;
			}
		}

		//get prob for this combination - owd * T * sofar (includes loss, alpha, constants)
		double p_combination = probOWD * probT * prob_so_far;

		gamma_prob += p_combination;

		//update new pmfs
		if (!skip_update) {
			if (p_combination > 0 && x_nPkts > 1) {
				(*new_t_pmf)[t_index] += p_combination;
				for (int m = 0; m < x_nPkts; m++) {
					//if (owd_indexes[m] != OWD_NUM_BINS-1) //skip if last bin

					if (owd_indexes[m] < OWD_NUM_BINS) //check if valid bin range
						(*new_owd_pmf)[owd_indexes[m]] += p_combination;
				}
			}
		}
	}

	subOS_prob += gamma_prob;
}

void HershelPlusEM::ProduceLossPatterns(int remaining, int start, int accumulatorSize){
	if (remaining > 0){
		for (int i = start; i < y_nPkts - remaining + 1; i++){
			bool possible = false;
			double diff = (*x_timestamps)[accumulatorSize] - (*y_timestamps)[i];
			if (diff >= 0){ // must arrive after the signature's timestamp		
				if (accumulatorSize == 0) // first packet being matched: allow it since the RTT can be anything
					possible = true;
				else {// non-first packet: do a test on jitter
					double jitt = (*x_timestamps)[accumulatorSize] - (*x_timestamps)[accumulatorSize - 1] -
						((*y_timestamps)[i] - (*y_timestamps)[accumulator[accumulatorSize - 1]]);
					// TODO: automatically determine the value of jitter for this cutoff 
					// need the jitter PMF and some threshold: CDF sum above the threshold should be less than some small number
					if (abs(jitt) < HP_JITTER_LOSS_THRESHOLD) // only if less than "JITTER_THRESHOLD" OWD
						possible = true;
				}
			}
			if (possible){
				accumulator[accumulatorSize] = i;
				ProduceLossPatterns(remaining - 1, i + 1, accumulatorSize + 1);
			}
		}
	}
	else {
		//if (produced == 0) // emulates always taking the first x[i].nPkts packets of the signature
		ExamineCombination(accumulator);
		//produced ++;
	}
}

int HershelPlusEM::options_match_int64(unsigned long long tar_opts, unsigned long long sig_opts){
	if (tar_opts == sig_opts) return 2;

	int count = 0;

	//keep track of where we found our last option in sig list
	int last_found_position = -1;
	unsigned long long sig_it;

	while (tar_opts > 0){
		int opt = -1;
		while (opt != OPT_MSS && opt != OPT_WIN && opt != OPT_SACK && opt != OPT_TS){ //skip unimportant bits
			if (tar_opts == 0) return 1; //we've matched so far and reached 0 without another option
			//get last 3 bits from target
			opt = tar_opts & 0x7; //AND by 3 bits to get value
			//move to next option
			tar_opts = tar_opts >> 3;
		}

		sig_it = sig_opts; //reset iterator to beginning
		int position = 0;
		bool found = false;
		while (sig_it > 0){ //go through sig_opt list
			int sopt = sig_it & 0x7;
			if (sopt == opt){ //if found
				if (position < last_found_position) return 0; //found it before the last one, doesnt match ordering
				else last_found_position = position;
				found = true;
				break;
			}
			sig_it = sig_it >> 3;
			position++;
		}
		if (!found){
			//not found in list, was probably enabled by user - do nothing here
		}
	}

	//all there
	return 1;
}

double HershelPlusEM::MatchConstantFeatures(Signature& target, Signature& dbsig, int osnum, vector<bool>* opt_flags){	
	//for constant feature match we just multiply by the probability in the pmf
	double total_prob = 1;	
	
	if (WIN_ENABLED) total_prob *= (*windowsize_pmf)[osnum][target.win];
	if (TTL_ENABLED) total_prob *= (*ttl_pmf)[osnum][ttl_map[target.ttl]];
	if (DF_ENABLED) total_prob *= (*df_pmf)[osnum][target.df];
	if (MSS_ENABLED) total_prob *= (*mss_pmf)[osnum][target.mss];
	
	
	if (OPT_ENABLED) {
		int result = options_match_int64(target.options_int, dbsig.options_int);
		if (result == 2 || result == 0) {
			total_prob *= (*opt_pmf)[osnum][opt_map[target.options_int]];
			(*opt_flags)[osnum] = true;
		}
		else {
			total_prob *= (*opt_pmf)[osnum][opt_map[target.option_order]];
			(*opt_flags)[osnum] = false;
		}
	}	

	
	if (RST_ENABLED) total_prob *= (*rst_pmf)[osnum][target.rst_int];

	return total_prob;
}

HershelPlusEM::HershelPlusEM(Signature& x, vector<Signature>& class_sigs, int osnum, double alpha_j, vector<double>* new_t_pmf, vector<double>* new_owd_pmf, vector<bool>* opt_flags)
		:x_nPkts(x.packet_arrival_timestamp.size()), y_nPkts(class_sigs[0].packet_arrival_timestamp.size()), totalprob(0), 
		new_t_pmf(new_t_pmf), new_owd_pmf(new_owd_pmf), sample_sig(x), dbnum(osnum)
{
	//return probability that x was generated from OS class_sigs
	//other variables help build the PMFs for the next iteration

	accumulator.resize(x_nPkts);
	x_timestamps = &(x.packet_arrival_timestamp);
	owd_indexes.resize(x_nPkts, 0);

	int lost_packets = y_nPkts - x_nPkts;
	
	//first match constant features
	//double constantprob = MatchConstantFeatures(x, class_sigs[0], osnum, opt_flags);
	double constantprob = MatchConstantFeatures(x, class_sigs[0], osnum, opt_flags);
	
	////for loss prob take from PMF but divide by (n choose k) to recover the original p^k (1-p)^(n-k) probability	
	double lossprob = (*loss_pmfs)[y_nPkts][lost_packets] / (*n_choose_k)[y_nPkts][lost_packets];

	prob_so_far = alpha_j * lossprob * constantprob;	

	if (lost_packets >= 0){
		double rtoprob = 0;

		//for each subOS, calculate and sum up the probability
		for (int r = 0; r < class_sigs.size(); r++){
			subOS_prob = 0;
			y_timestamps = &(class_sigs[r].packet_arrival_timestamp);
			//fill(accumulator.begin(), accumulator.end(), 0);

			ProduceLossPatterns(x_nPkts, 0, 0);
			rtoprob += subOS_prob;
		}

		totalprob = rtoprob / class_sigs.size();
	}


	//totalprob = alpha_j * constantprob;
}
