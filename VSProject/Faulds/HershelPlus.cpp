/*
Implementation of Hershel+ to be used for HardEM testing
Hershel+ Optimized uses closed form of Hershel+ equation, as opposed to calculating the integral numerically
Author: Zain Shamsi
*/


#include "classes.h"

vector<double> HershelPlus::owd_pmf;
vector<double> HershelPlus::t_pmf;

void HershelPlus::ExamineCombination(vector<int>& gamma) {
	double gamma_prob = 0;
	double t_limit = T_GUESS_LIMIT_IN_SECS;
	double t_step_size = T_GUESS_LIMIT_IN_SECS / T_NUM_BINS;
	double owd_step_size = OWD_GUESS_LIMIT_IN_SECS / OWD_NUM_BINS;

	//but T has to be less than the smallest distance between any x and y
	for (int t = 0; t < x_nPkts; t++) {
		double d = x_timestamps[t] - y_timestamps[gamma[t]];
		if (d < t_limit) t_limit = d;
	}

	for (double t_guess = 0; t_guess <= t_limit; t_guess += t_step_size) {
		double probQ = 1;

		int t_index = t_guess / t_step_size;
		double probT = t_pmf[t_index];

		//for each received packet, determine the OWD (Q) and calculate \prod p(Q)
		for (int t = 0; t < x_nPkts; t++) {
			double Q = (x_timestamps[t] - y_timestamps[gamma[t]] - t_guess);
			if (Q < 0) {
				probQ = 0;
				break;
			}
			else {
				int owd_index = Q / owd_step_size;
				if (owd_index > owd_pmf.size() - 1) owd_index = owd_pmf.size() - 1;
				double p = owd_pmf[owd_index];
				probQ *= p;
			}
		}

		//the prob of this gamma with this T is p(T) * \prod p(Q)
		gamma_prob += probQ * probT;
	}

	subOS_prob += gamma_prob;
}

void HershelPlus::ExamineCombinationOptimized(vector<int>& gamma){
	double gamma_prob = 1;
	double s = DBL_MAX; //keep track of the min e_m seen
				
	//for each received packet, determine the OWD (Q) and calculate \prod p(Q)
	for (int m = 0; m < x_nPkts; m++){
		double em = x_timestamps[m] - y_timestamps[gamma[m]];
		if (em < 0){
			gamma_prob = 0;
		}
		else {
			int owd_index = em / (OWD_GUESS_LIMIT_IN_SECS / OWD_NUM_BINS);
			if (owd_index > owd_pmf.size() - 1) owd_index = owd_pmf.size() - 1;
			gamma_prob *= owd_pmf[owd_index];

			//gamma_prob *= exp(-HP_JITTER_LAMBDA * em);

			if (em < s) s = em;
		}
	}

	if (gamma_prob > 0){
		//have case for denom = 0
		if (x_nPkts == HP_ERLANG_NU / HP_OWD_LAMBDA) gamma_prob *= (s * s) / 2;
		else gamma_prob *= (exp(-HP_ERLANG_NU * s) * (exp(HP_OWD_LAMBDA * s * x_nPkts) * (HP_OWD_LAMBDA * s * x_nPkts - HP_ERLANG_NU * s - 1) + exp(HP_ERLANG_NU * s))) / pow((HP_ERLANG_NU - HP_OWD_LAMBDA * x_nPkts), 2); //result from Wolfram
		//else gamma_prob *= (1 - exp(min_em * (-HP_GAMMA_NU * HP_JITTER_LAMBDA * x_nPkts)) * (1 + min_em * (HP_GAMMA_NU - HP_JITTER_LAMBDA * x_nPkts))) / pow((HP_GAMMA_NU - HP_JITTER_LAMBDA * x_nPkts), 2); // to match paper (16)
			
		//multiply by the constant for completeness
		gamma_prob *= HP_OWD_LAMBDA * pow(HP_ERLANG_NU, 2);
		//gamma_prob *= pow(HP_JITTER_LAMBDA, x_nPkts) * pow(HP_GAMMA_NU, 2); //to match paper (16)

	}

	subOS_prob += gamma_prob;
}

void HershelPlus::ProduceLossPatterns(int remaining, int start, int accumulatorSize){
	if (remaining > 0){
		for (int i = start; i < y_nPkts - remaining + 1; i++){
			bool possible = false;
			double diff = x_timestamps[accumulatorSize] - y_timestamps[i];
			if (diff >= 0){ // must arrive after the signature's timestamp		
				if (accumulatorSize == 0) // first packet being matched: allow it since the RTT can be anything
					possible = true;
				else {// non-first packet: do a test on jitter
					double jitt = x_timestamps[accumulatorSize] - x_timestamps[accumulatorSize - 1] -
						(y_timestamps[i] - y_timestamps[accumulator[accumulatorSize - 1]]);
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
		ExamineCombination(accumulator); //H+
		//ExamineCombinationOptimized(accumulator); //H+ Optimized
	}
}


/*********************CONSTANT FEATURE MATCHING***********************/

//match int64 options, intersection of both tar and sig should be ordered the same
int HershelPlus::options_match_int64(unsigned long long tar_opts, unsigned long long sig_opts){
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

//calculate likeliest class by adding FEATURE_CHANGE_PROBABILITY onto candidates
double HershelPlus::MatchConstantFeatures(Signature& target, Signature& dbsig){
	//start with total prob = 1
	double total_prob = 1;

	if (target.win == dbsig.win) total_prob *= (1 - HP_FEATURE_PROB_PI);
	else total_prob *= HP_FEATURE_PROB_PI;

	if (target.ttl == dbsig.ttl) total_prob *= (1 - HP_FEATURE_PROB_PI);
	else total_prob *= HP_FEATURE_PROB_PI;

	if (target.df == dbsig.df) total_prob *= (1 - HP_FEATURE_PROB_PI);
	else total_prob *= HP_FEATURE_PROB_PI;

	int oval = options_match_int64(target.options_int, dbsig.options_int);
	//if (oval <= 1) total_prob *= HP_FEATURE_PROB_PI_RST_OPT;
	if (oval == 0) total_prob *= HP_FEATURE_PROB_PI_RST_OPT;
	if (oval == 1) total_prob *= HP_FEATURE_PROB_PI;
	else if (oval > 1) total_prob *= (1 - HP_FEATURE_PROB_PI);

	if (target.mss == dbsig.mss) total_prob *= (1 - HP_FEATURE_PROB_PI);
	else total_prob *= HP_FEATURE_PROB_PI;

	if (target.rst == dbsig.rst && target.rst_ack == dbsig.rst_ack	&& target.rst_win == dbsig.rst_win && target.rst_nonzero == dbsig.rst_nonzero) total_prob *= (1 - HP_FEATURE_PROB_PI_RST_OPT);
	else total_prob *= HP_FEATURE_PROB_PI_RST_OPT;			

	return total_prob;
}

/***********************************************************************/

HershelPlus::HershelPlus(Signature& observation, vector<Signature>& class_sigs)
	:x_nPkts(observation.packet_arrival_timestamp.size()), y_nPkts(class_sigs[0].packet_arrival_timestamp.size()), totalprob(1){

	accumulator.resize(x_nPkts);
	x_timestamps = observation.packet_arrival_timestamp;		
	
	//database RTO vector must be >= size of sample, due to loss
	int lost_packets = class_sigs[0].packet_arrival_timestamp.size() - observation.packet_arrival_timestamp.size();

	if (lost_packets >= 0) {

		//start with constant feature prob
		totalprob = MatchConstantFeatures(observation, class_sigs[0]);

		//add loss probability
		totalprob *= pow(HP_LOSS_PROB, lost_packets) * pow(1 - HP_LOSS_PROB, observation.packet_arrival_timestamp.size());

		//for each subOS tj of dbsig
		double rtoProb = 0;
		for (int tj = 0; tj < class_sigs.size(); tj++) {
			//calculate OWD probability using Estimator
			subOS_prob = 0;
			y_timestamps = class_sigs[tj].packet_arrival_timestamp;
			ProduceLossPatterns(observation.packet_arrival_timestamp.size(), 0, 0);
			rtoProb += subOS_prob;
		}

		//get average over all subOS
		rtoProb /= class_sigs.size();

		//add subOS probability
		totalprob *= rtoProb;

	}
}
