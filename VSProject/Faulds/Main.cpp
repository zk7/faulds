//    Copyright © 2017 IRL at Texas A&M University (http://irl.cse.tamu.edu)
//
//    Faulds is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Faulds is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details
//    http://www.gnu.org/licenses/lgpl.txt.
//
//    Publication:
//	  Z. Shamsi, D.B.H. Cline and D. Loguinov, "Faulds: A Non-Parametric Iterative
//	  Classifier for Internet-Wide OS Fingerprinting" ACM CCS, Nov 2017.
//

#include <time.h>
#include <fstream>

#include "classes.h"
#include "SimpleRNG.h"

#define SIMULATE_RATIOS {0.9, 0.05, 0.05}	
#define SIMULATE_T_MEAN 0.5
#define SIMULATE_DELAY_MEAN 0.5
#define SIMULATE_LOSS 0.1
#define SIMULATE_INDEPENDENT_FEATURE_CHANGE 1
#define CPI 0.8 //probability to stay at known value
#define OPT_PI_SUBSET 1
#define SIMULATE_PI {CPI, CPI, CPI, CPI, CPI, CPI}
//#define SIMULATE_PI {0.3, 0.2, 0.5, 0.4, 0.4, CPI}

#ifdef _DEBUG
#define NUM_DISTURBANCE 100
#else 
#define NUM_DISTURBANCE (1 << 18)
#endif

SimpleRNG srng;
//unordered_map<int, int> win_map;
unordered_map<int, int> ttl_map = { {32, 0}, {64, 1}, {128, 2}, {255, 3} };
unordered_map<long long int, int> opt_map;
//unordered_map<int, int> mss_map;
//unordered_map<int, int> rst_map;
unordered_map<int, unordered_set<int>> option_subsets;
unordered_map<unsigned int, string> IP_to_OS_map;


vector<int> getOptionOrder(long long int options, int& orderedopts) {
	//map number and order of options seen
	vector<int> options_seen;
	orderedopts = 0;

	while (options > 0) {
		int opt = options & 0x7; //AND by 3 bits to get option
		while (opt != OPT_MSS && opt != OPT_WIN && opt != OPT_SACK && opt != OPT_TS) { //skip unimportant bits			
			options = options >> 3; //move to next option
			opt = options & 0x7; //get next option
			if (opt == 0) break;
		}

		if (opt > 0) {
			//store ordering of options
			options_seen.push_back(opt);
		}

		options = options >> 3;
	}

	//build option order
	reverse(options_seen.begin(), options_seen.end()); //reverse to maintain original ordering in the longlong int

	for (int i = 0; i < options_seen.size(); i++) {
		orderedopts <<= 3;
		orderedopts = orderedopts | options_seen[i];
	}

	//output
	//printf("\noption order -> ");
	//int x = orderedopts;
	//while (x > 0) {
	//	int o = x & 7;
	//	if (o == OPT_MSS) printf("M");
	//	else if (o == OPT_WIN) printf("W");
	//	else if (o == OPT_SACK) printf("S");
	//	else if (o == OPT_TS) printf("T");
	//	else printf("What?");
	//	x = x >> 3;
	//}

	return options_seen;
}

bool containMSS(long long int options) {
	while (options > 0) {
		int opt = options & 0x7; //AND by 3 bits to get option		
		if (opt == OPT_MSS) return true;

		options = options >> 3; //move to next option		
	}
	return false;
}

unsigned long long generateRandomOptions() {
	//generate random valid options string
	unsigned long long new_opts = 0;
	vector<int> valid_options = { OPT_MSS, OPT_WIN, OPT_TS, OPT_SACK };
	int max = valid_options.size();
	while (max > 0) {
		//pick random option
		int index = srng.GetUint() % max;
		int option = valid_options[index];

		//randomly decide whether to add it
		int add = srng.GetUint() % 2;

		if (add) {
			//add NOPs
			int num_nops = 0;
			if (option == OPT_WIN) num_nops = 1;
			if (option == OPT_SACK || OPT_TS) num_nops = 2;

			//flip a coin to add before or after
			int nops_before = srng.GetUint() % 2;

			//push NOPS and option in order
			if (nops_before) {
				while (num_nops > 0) {
					new_opts <<= 3;
					new_opts = new_opts | OPT_NOP;
					num_nops--;
				}
			}

			new_opts <<= 3;
			new_opts = new_opts | option;

			if (!nops_before) {
				while (num_nops > 0) {
					new_opts <<= 3;
					new_opts = new_opts | OPT_NOP;
					num_nops--;
				}
			}
		}

		//swap in the last one and reduce max to remove this option from the list
		valid_options[index] = valid_options[max - 1];
		max--;
	}
	return new_opts;
}

int getRSTInt(int rst, int rst_ack, int rst_win, int rst_nonzero) {
	//push bits into one int and return
	int retval = 0;
	retval = retval | (rst & 1);
	retval <<= 1;
	retval = retval | (rst_ack & 1);
	retval <<= 2;
	retval = retval | (rst_win & 3);
	retval <<= 1;
	retval = retval | (rst_nonzero & 1);

	return retval;
}

unordered_set<int> generateSubsets(vector<int> options) {

	int num_subsets = pow(2, options.size()) - 1; //dont take the "everything on" case

	//now generate subsets, iterating over all possibilities in range, and looking at bits that are set
	unordered_set<int> subsets;
	for (int i = 0; i < num_subsets; i++) {
		int sub = 0;

		//look at each bit in i
		for (int b = 0; b < options.size(); b++) {
			int x = 1 << b; //bit to look at
			if (i & x) { //if bit is set, add opt to subset
				sub <<= 3;
				sub = sub | options[b];
			}
		}

		//output for testing
		/*printf("\nsubset -> ");
		int x = sub;
		while (x > 0) {
			int o = x & 7;
			if (o == OPT_MSS) printf("M");
			else if (o == OPT_WIN) printf("W");
			else if (o == OPT_SACK) printf("S");
			else if (o == OPT_TS) printf("T");
			else printf("What?");
			x = x >> 3;
		}*/

		subsets.insert(sub);

	}

	return subsets;
}

unordered_map<unsigned int, string> readMapping(char *filename) {
	unordered_map<unsigned int, string> retmap;
	int BUFFER_SIZE = 8192;

	FILE* fin;

	int ret = fopen_s(&fin, filename, "r");
	if (ret != 0) {
		printf("Error opening file %s!\n", filename);
		exit(-1);
	}

	char* buffer = new char[BUFFER_SIZE];
	unsigned int id;
	char OS[4096];

	while (!feof(fin)) {
		//read next line
		fgets(buffer, BUFFER_SIZE, fin);
		char* bufferptr = buffer;


		sscanf_s(bufferptr, "%u,%[^\n]\n", &id, OS, _countof(OS));
		retmap[id] = string(OS);
	}

	printf("Read %d OS name mappings...\n", retmap.size());

	fclose(fin);
	return retmap;
}

map<u_int, vector<Signature>> readSigs(char *filename, int& max_rto_length){
	map<u_int, vector<Signature>> retmap;
	int sig_count = 0;
	int BUFFER_SIZE = 8192;
	max_rto_length = 0;

	FILE* fin;

	int ret = fopen_s(&fin, filename, "r");
	if (ret != 0){
		printf("Error opening file %s!\n", filename);
		exit(-1);
	}
	char* buffer = new char[BUFFER_SIZE];
	int count = 0;
	double timestamp = 0;
	vector<double> rtos;

	while (!feof(fin)){
		Signature sig;

		//read next line
		fgets(buffer, BUFFER_SIZE, fin);
		char* bufferptr = buffer;

		sscanf_s(bufferptr, "%u,%d,%d,%d,%[^,],%lld,%d,%d,%d,%d,%d,%d%n",
			&sig.ip_int, &sig.win, &sig.ttl, &sig.df, &sig.options, _countof(sig.options), &sig.options_int, &sig.mss, &sig.rst, &sig.rst_ack, &sig.rst_win, &sig.rst_seq, &sig.rst_nonzero, &count);

		
		bufferptr += count;

		while (strcmp(bufferptr, "\n") != 0 && strlen(bufferptr) > 0){
			sscanf_s(bufferptr, ",%lf%n", &timestamp, &count);
			sig.packet_arrival_timestamp.push_back(timestamp);

			if (sig.packet_arrival_timestamp.size() > 20) break;

			bufferptr += count;
		}

		sig.was_correct = false;
		sig.predicted_ip_int = 0;
		
		//fix rst_win to ternary
		if (sig.rst_win == sig.win) sig.rst_win = 1;
		else if (sig.rst_win == 0) sig.rst_win = 0;
		else sig.rst_win = 2;

		//convert rst to int
		sig.rst_int = getRSTInt(sig.rst, sig.rst_ack, sig.rst_win, sig.rst_nonzero);

		//track max arrival size
		if (sig.packet_arrival_timestamp.size() > max_rto_length) max_rto_length = sig.packet_arrival_timestamp.size();

		//generate option subsets
		int original_order;
		vector<int> opt = getOptionOrder(sig.options_int, original_order);
		sig.option_order = original_order;
		if (option_subsets.count(original_order) == 0) {
			unordered_set<int> subsets = generateSubsets(opt);
			option_subsets[original_order] = subsets;
		}

		//add feature values to maps, if not already there. We will keep track of these as the set of all possible for EM
		//if (win_map.count(sig.win) == 0) win_map[sig.win] = win_map.size();
		//if (mss_map.count(sig.mss) == 0) mss_map[sig.mss] = mss_map.size();
		//if (rst_map.count(sig.rst_int) == 0) rst_map[sig.rst_int] = rst_map.size();
		//add database opt string, order string and all subsets to pmf, if not already there
		if (opt_map.count(sig.options_int) == 0) opt_map[sig.options_int] = opt_map.size();
		if (opt_map.count(sig.option_order) == 0) opt_map[sig.option_order] = opt_map.size();
		for (int sub : option_subsets[sig.option_order]) {
			if (opt_map.count(sub) == 0) opt_map[sub] = opt_map.size();
		}

		retmap[sig.ip_int].push_back(sig);
		sig_count++;
		if (sig_count % 10000 == 0) printf("Read %d signatures...\r", sig_count);
	}

	fclose(fin);
	printf("\nStored %zd signatures in map\n", retmap.size());
	return retmap;
}

void compactList(map<u_int, vector<Signature>>& siglist, double threshold) {
	double CLOSENESS_THRESHOLD = threshold;
	int remove_count = 0;

	for (auto pair : siglist) {
		
		//build a vector of subOSes seen for this sig
		//start with the first one automatically
		vector<Signature> primary_subOSes;
		primary_subOSes.push_back(pair.second[0]);

		//go through all other subOSes
		for (int s = 1; s < pair.second.size(); s++) {
			bool found_match = false;

			//match this subOS to the primaries
			for (Signature prim : primary_subOSes){

				//match timestamps to threshold
				bool too_far = false;
				for (int r = 0; r < pair.second[s].packet_arrival_timestamp.size(); r++) {
					if (fabs(pair.second[s].packet_arrival_timestamp[r] - prim.packet_arrival_timestamp[r]) > CLOSENESS_THRESHOLD) {
						too_far = true;
						break;
					}
				}

				if (!too_far) {
					found_match = true;
					break;
				}

			}

			if (!found_match) { //if did not match, this becomes a primary
				primary_subOSes.push_back(pair.second[s]);
			}

		}	

		remove_count += pair.second.size() - primary_subOSes.size();

		//replace the subOSes for this signature
		siglist[pair.first] = primary_subOSes;

	}

	printf("Removed %d subOSes from DB\n", remove_count);
}

vector<Signature> readTestData(char *filename, int max_rto_length, double subsample) {
	vector<Signature> retlist;
	int sig_count = 0;
	int skip_count = 0;
	int BUFFER_SIZE = (1 << 26);

	/*
	map<int, double> windows;
	map<int, double> mss;
	map<u_int, double> options;
	map<int, double> reset;
	*/

	FILE* fin;

	int ret = fopen_s(&fin, filename, "r");
	if (ret != 0) {
		printf("Error opening file %s!\n", filename);
		exit(-1);
	}
	char* buffer = new char[BUFFER_SIZE];
	int count = 0;
	double timestamp = 0;
	int msscount = 0;
	vector<double> rtos;

	printf("\nReading Test Data\n");
	while (!feof(fin)) {
		Signature sig;
		bool skip = false;

		//read next line
		fgets(buffer, BUFFER_SIZE, fin);
		char* bufferptr = buffer;

		//sscanf_s(bufferptr, "%[^,]", sig.options, _countof(sig.options)
		//sscanf_s(bufferptr, "%u,%d,%d,%d,%d,%d,%d,%d,%d,%d,%llu%n",
		//	&sig.ip_int, &sig.win, &sig.ttl, &sig.df, &sig.rst, &sig.rst_ack, &sig.rst_win, &sig.rst_seq, &sig.rst_nonzero, &sig.mss, &sig.options_int, &count);

		sscanf_s(bufferptr, "%u,%d,%d,%d,%[^,],%lld,%d,%d,%d,%d,%d%n",
			&sig.ip_int, &sig.win, &sig.ttl, &sig.df, &sig.options, _countof(sig.options), &sig.options_int, &sig.mss, &sig.rst, &sig.rst_ack, &sig.rst_win, &sig.rst_nonzero, &count);

		bufferptr += count;

		while (strcmp(bufferptr, "\n") != 0 && strlen(bufferptr) > 0) {
			sscanf_s(bufferptr, ",%lf%n", &timestamp, &count);
			sig.packet_arrival_timestamp.push_back(timestamp);

			//only care about signatures we can classify
			if (sig.packet_arrival_timestamp.size() > max_rto_length) { skip = true; break; }

			bufferptr += count;
		}

		//fix rst_win to ternary
		if (sig.rst_win == sig.win) sig.rst_win = 1;
		else if (sig.rst_win == 0) sig.rst_win = 0;
		else sig.rst_win = 2;

		//convert rst to int
		sig.rst_int = getRSTInt(sig.rst, sig.rst_ack, sig.rst_win, sig.rst_nonzero);

		//generate option subsets
		int original_order;
		vector<int> opt = getOptionOrder(sig.options_int, original_order);
		sig.option_order = original_order;
		if (option_subsets.count(original_order) == 0) {
			unordered_set<int> subsets = generateSubsets(opt);
			option_subsets[original_order] = subsets;
		}

		//add feature values to maps, if not already there
		//if (win_map.count(sig.win) == 0) win_map[sig.win] = win_map.size();
		//if (mss_map.count(sig.mss) == 0) mss_map[sig.mss] = mss_map.size();
		//if (rst_map.count(sig.rst_int) == 0) rst_map[sig.rst_int] = rst_map.size();
		//add database opt string, order string and all subsets to pmf, if not already there
		if (opt_map.count(sig.options_int) == 0) opt_map[sig.options_int] = opt_map.size();
		if (opt_map.count(sig.option_order) == 0) opt_map[sig.option_order] = opt_map.size();
		for (int sub : option_subsets[sig.option_order]) {
			if (opt_map.count(sub) == 0) opt_map[sub] = opt_map.size();
		}

		sig.predicted_ip_int = 0;
		sig.predicted_winner_index = 0;

		if (sig.mss > 1500) {
			msscount++;
		}

		if (!skip) {
			if (srng.GetUniform() < subsample)
				retlist.push_back(sig);
		}
		else skip_count++;

		sig_count++;
		if (sig_count % 10000 == 0) {
			printf("Read %d signatures, skipped %d...\r", sig_count, skip_count);
			//break;
		}
	}

	/*
	read distribution of values from test data
	FILE* fout;
	fopen_s(&fout, "winvals.txt", "w");
	double sum = 0;
	double sofar = 0;

	for (auto p : windows) sum += p.second;
	for (auto p : windows) {
		double freq = p.second / sum;
		double cdfprob = sofar + freq;
		fprintf(fout, "%d, %lf\n", p.first, cdfprob);
		sofar += freq;
	}
	fclose(fout);

	fopen_s(&fout, "mssvals.txt", "w");
	sum = 0; 
	sofar = 0;
	for (auto p : mss) sum += p.second;
	for (auto p : mss) {
		double freq = p.second / sum;
		double cdfprob = sofar + freq;
		fprintf(fout, "%d, %lf\n", p.first, cdfprob);
		sofar += freq;
	}
	fclose(fout);

	fopen_s(&fout, "optvals.txt", "w");
	sum = 0;
	sofar = 0;
	for (auto p : options) sum += p.second;
	for (auto p : options) {
		double freq = p.second / sum;
		double cdfprob = sofar + freq;
		fprintf(fout, "%u, %lf\n", p.first, cdfprob);
		sofar += freq;
	}
	fclose(fout);

	fopen_s(&fout, "rstvals.txt", "w");
	sum = 0;
	sofar = 0;
	for (auto p : reset) sum += p.second;
	for (auto p : reset) {
		double freq = p.second / sum;
		double cdfprob = sofar + freq;
		fprintf(fout, "%d, %lf\n", p.first, cdfprob);
		sofar += freq;
	}
	fclose(fout);
	*/

	fclose(fin);
	printf("Found %d signatures (%.2lf%%) with MSS > 1500\n", msscount, ((double)msscount / retlist.size()) * 100);
	printf("Stored %zd signatures in list\n", retlist.size());
	printf("Skipped adding %d signatures\n", skip_count);
	return retlist;
}

void injectUnknown(map<u_int, vector<Signature>>& original_sigs, double fraction) {	
	//injects a manually created signature into the testdata
	//int num_inject = NUM_DISTURBANCE * fraction; //times the dataset	
	if (fraction > 1) {
		printf("Injection fraction greater than 1!\n");
		exit(0);
	}
	int num_inject = NUM_DISTURBANCE * (1 / (1 / fraction - 1)); //by percentage of dataset
	

	Signature newSig;
	newSig.ip_int = UINT_MAX;
	//add it to the original sigs vector
	original_sigs[newSig.ip_int] = { newSig };

	original_sigs[newSig.ip_int][0].prob = num_inject;
	original_sigs[newSig.ip_int][0].win = 8888;
	original_sigs[newSig.ip_int][0].ttl = 255;
	original_sigs[newSig.ip_int][0].df = 1;
	strcpy_s(original_sigs[newSig.ip_int][0].options, "TW");
	original_sigs[newSig.ip_int][0].options_int = 26; //011 010
	original_sigs[newSig.ip_int][0].mss = 888;
	original_sigs[newSig.ip_int][0].rst = 1;
	original_sigs[newSig.ip_int][0].rst_ack = 1;
	original_sigs[newSig.ip_int][0].rst_win = 0;
	original_sigs[newSig.ip_int][0].rst_nonzero = 0;
	original_sigs[newSig.ip_int][0].rst_int = getRSTInt(1, 1, 0, 0);

	original_sigs[newSig.ip_int][0].ip_int = UINT_MAX;

	//original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(0);
	//original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(10);
	//original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(20);
	//original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(30);

	original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(0);
	original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(2.0);
	//original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(3.0);
	original_sigs[newSig.ip_int][0].packet_arrival_timestamp.push_back(4.0);
}

void shrinkDatabase(map<u_int, vector<Signature>>& original_sigs, vector<Signature> samples, double fraction) {
	//randomly shrink database by a fraction for testing unknown OSes in the dataset
	printf("Shrinking Database by %.1lf%%!\n", fraction * 100);
	int num = original_sigs.size() * fraction;

	unordered_set<u_int> removed;

	map<u_int, vector<Signature>>::iterator it;
	for (int i = 0; i < num; i++) {
		it = original_sigs.begin();

		int r = srng.GetUint() % original_sigs.size();
		advance(it, r);
		
		original_sigs.erase(it);
		removed.insert(it->first);
	}

	//count removed from samples
	int count = 0;
	for (Signature s : samples) {
		if (removed.count(s.ip_int) > 0) count++;
	}

	printf("Done shrinking. New database size is %d\n", original_sigs.size());
	printf("%d (%.4lf%%) samples are now unknown.\n", count, (double)count / samples.size() * 100);
}

vector<Signature> createDisturbedZipf(map<u_int, vector<Signature>> original_sigs, double t_mean, double delay_mean, double loss_prob, int max_rto_length, bool shuffle){
	/////////////////////GENERATE ALPHAS////////////////////////
	//first build Zipf pmf
	vector<double> ratios;
	double sum = 0;
	int size = original_sigs.size();
	if (original_sigs.count(UINT_MAX) > 0) size--; //if there is an injected, dont count it for ratio
	for (int j = 0; j < size; j++){
		ratios.push_back(1.0 / pow((j + 1), 1.2));
		sum += ratios[j];
	}
	
	//normalize ratios
	for (int i = 0; i < size; i++)	ratios[i] /= sum;
	
	//shuffle ratios if desired
	if (shuffle) {
		random_shuffle(ratios.begin(), ratios.end());
	}

	//overwrite ratios if needed
	if (original_sigs.size() <= 4) ratios = SIMULATE_RATIOS;

	if (ratios.size() != original_sigs.size()) {
		printf("Ratio size mismatch! Database: %zd, Ratios, %zd\n", original_sigs.size(), ratios.size());
		exit(0);
	}

	/////////////GENERATE LOSS PMF///////////////////////
	vector<vector<double>> generating_loss_cdfs(max_rto_length + 1, vector<double>()); //loss generating cdf

	for (int n = 0; n < generating_loss_cdfs.size(); n++){
		//resize each cdf accordingly
		generating_loss_cdfs[n].resize(n+1);
				
		//then set up for losing k packets out of n
		//binomial pmf
		/*for (int k = 0; k <= n; k++){
			double a = factorial(n) / (factorial(k) * factorial(n - k));
			double prob_lose_k_n = a * pow(loss_prob, k) * pow(1 - loss_prob, n - k);
			if (k == 0) generating_loss_cdfs[n][k] = prob_lose_k_n;
			else generating_loss_cdfs[n][k] = prob_lose_k_n + generating_loss_cdfs[n][k - 1];			
		}*/

		//reverse binomial pmf		
		if (n > 0) {
			for (int k = 0, reverse = n - 1; k < n; k++, reverse--) {
				double a = factorial(n) / (factorial(reverse) * factorial(n - reverse));
				double prob_lose_rev_n = a * pow(loss_prob, reverse) * pow(1 - loss_prob, n - reverse);
				if (k == 0) generating_loss_cdfs[n][k] = prob_lose_rev_n;
				else generating_loss_cdfs[n][k] = prob_lose_rev_n + generating_loss_cdfs[n][k - 1];
			}
		}

		//hardcoded to lose all but one packet if loss exists
		//for (int k = 0; k <= n; k++) {
		//	generating_loss_cdfs[n][k] = 1 - loss_prob;
		//	if (k >= n-1) generating_loss_cdfs[n][k] = 1;			
		//}
	}

	//now we read in a list of possible values for fixed features, gathered from the Internet scan data
	//files are formatted as value,cumulative_probability
	vector<int> winvals;
	vector<int> mssvals;
	vector<long long int> optvals;
	vector<int> rstvals;

	FILE* fin;
	int value;
	double prob;
	int ret = fopen_s(&fin, "winvals.txt", "r");
	if (ret != 0) {
		printf("winvals.txt not found!\n");
		exit(0);
	}
	while (!feof(fin)) {
		//read next line
		fscanf_s(fin, "%d,%lf\n", &value, &prob);
		winvals.push_back(value);
	}
	fclose(fin);

	fopen_s(&fin, "mssvals.txt", "r");
	while (!feof(fin)) {
		//read next line
		fscanf_s(fin, "%d,%lf\n", &value, &prob);
		mssvals.push_back(value);
	}
	fclose(fin);

	fopen_s(&fin, "optvals.txt", "r");
	long long int longvalue;
	while (!feof(fin)) {
		//read next line
		fscanf_s(fin, "%lld,%lf\n", &longvalue, &prob);
		optvals.push_back(longvalue);
	}
	fclose(fin);

	fopen_s(&fin, "rstvals.txt", "r");
	while (!feof(fin)) {
		//read next line
		fscanf_s(fin, "%d,%lf\n", &value, &prob);
		rstvals.push_back(value);
	}
	fclose(fin);

	//set up vars and tracking of the real distributions
	int loss_counter = 0;
	int feature_changed_counter = 0;
	bool rst_lost = false;
	int total_packets_lost = 0;
	double total_packets = 0;
	int total_features_changed = 0;
	int discarded = 0;
	double total_features = 0;
	vector<Signature> ret_list;
	vector<double> real_t_dist(T_NUM_BINS, 0);
	vector<double> real_q_dist(OWD_NUM_BINS, 0);
	vector<double> real_win_dist(65536, 0);
	vector<double> real_ttl_dist(4, 0);	
	vector<double> real_df_dist(2, 0);
	vector<double> real_mss_dist(65536, 0);
	vector<vector<double>> real_loss_dists(max_rto_length + 1, vector<double>());
	for (int n = 0; n < real_loss_dists.size(); n++) real_loss_dists[n].resize(n + 1, 0);

	///////////GENERATE SAMPLES////////////////////////////////////
	printf("Generating %d signatures...\n", NUM_DISTURBANCE);
	map<u_int, vector<Signature>>::iterator it = original_sigs.begin();
	for (int j = 0; it != original_sigs.end(); ++it, ++j){
		int totalsamplesj;
		if (it->first == UINT_MAX) {
			totalsamplesj = it->second[0].prob;
			printf("Injecting %d samples!\n", totalsamplesj);
		}
		else totalsamplesj = NUM_DISTURBANCE * ratios[j];

		//cycle through subOSes until reached totalsamples
		for (int s = 0, i = 0; i < totalsamplesj; s = ((s + 1) % it->second.size()), i++){

			Signature sig(it->second[s]);
			//if (it->first == UINT_MAX) {
			//	ret_list.push_back(sig);
			//	continue;
			//}
			sig.packet_arrival_timestamp.clear();
			sig.prob = 0;
			sig.predicted_ip_int = 0;		
			sig.was_correct = false;

			///////////////////NETWORK FEATURES//////////////////////
			
			//double T = srng.GetUniform(0, t_mean * 2); //generate random server processing time + sending delay as one uniform var
			//double T = srng.GetExponential(t_mean);
			//double T = srng.GetErlang(t_mean, 2);
			double T = srng.GetReverseExponential(t_mean, 2);
			//double T = srng.GetPareto(3, t_mean);
			//double T = srng.GetGamma(2.0, 0.25);

			int t_index = T / T_STEP_SIZE;
			if (t_index < real_t_dist.size()) real_t_dist[t_index]++;
			
			for (int m = 0; m < it->second[s].packet_arrival_timestamp.size(); m++){
				double qi;
				qi = srng.GetExponential(delay_mean); //from http://tcp-reassembly-project.googlecode.com/files/jitter_model.pdf 		
				//qi = srng.GetReverseExponential(delay_mean, 2);
				//qi = srng.GetPareto(3, delay_mean);
				//qi = srng.GetUniform(0, delay_mean * 2);
				qi = srng.GetErlang(delay_mean, 2);
				//qi = 0;

				int q_index = qi / OWD_STEP_SIZE;	
				if (q_index < real_q_dist.size()) real_q_dist[q_index]++; //q_index = real_q_dist.size() - 1;
				
				double departure = T + it->second[s].packet_arrival_timestamp[m]; //departure of packet from server
				double arrival = departure + qi; //arrival of packet at client with OWD

				if (sig.packet_arrival_timestamp.empty()){
					//if no packets till now, push this one back as first packet
					sig.packet_arrival_timestamp.push_back(arrival);
				}
				else { //check for FIFO ordering
					double last_arrival = sig.packet_arrival_timestamp.back();
					//if new packet arrived after arrival of last packet, add like normal
					//else, new packet arrived at same time or before arrival of last one, make arrival equal to last_arrival to preserve order
					//this effectively pushes 0 as the RTO
					if (arrival > last_arrival) sig.packet_arrival_timestamp.push_back(arrival);
					else sig.packet_arrival_timestamp.push_back(last_arrival);
				}
			}

			//decide how many packets to lose
			int packets_lost = 0;
			double loss = srng.GetUniform();
			for (int k = 0; k < generating_loss_cdfs[it->second[s].packet_arrival_timestamp.size()].size(); k++){
				if (loss < generating_loss_cdfs[it->second[s].packet_arrival_timestamp.size()][k]){
					packets_lost = k;
					break;
				}
			}
			total_packets_lost += packets_lost; //keep track of total packets lost
			total_packets += it->second[s].packet_arrival_timestamp.size(); //keep track of total packets observed

			//randomly drop that many packets, keeping track of rst loss
			rst_lost = false;
			for (int lost = 0; lost < packets_lost; lost++){
				int random_position = srng.GetUint() % sig.packet_arrival_timestamp.size();
				if (sig.rst && random_position == it->second[s].packet_arrival_timestamp.size() - 1) rst_lost = true;
				sig.packet_arrival_timestamp.erase(sig.packet_arrival_timestamp.begin() + random_position);
			}

			if (packets_lost) loss_counter++;
			
			////////////////////USER FEATURES/////////////////
			//randomly change constant features. Implemented are 4 ways to do this:
			//1. RAND - randomly pick a value for each feature from its full feature space
			//2. HALFWAY - randomly pick a value for each feature from all possible values seen on the Internet (read from files above)
			//3. POPULAR - pick another OS and copy its features. If it doesnt change, pick feature from second popular OS
			//4. UNIQUE - pick a value that is unique for each OS as much as it can be
			int RAND = 1;
			int HALFWAY = 0;
			int POPULAR = 0;
			int pop1 = 1;
			int pop2 = 3;
			int UNIQUE = 0;

			vector<double> real_pi = SIMULATE_PI;
			bool feature_changed = false;


			double change_prob = srng.GetUniform(0.0, 1.0);
			for (int feature = 0; feature < real_pi.size(); feature++) {
				if (SIMULATE_INDEPENDENT_FEATURE_CHANGE) change_prob = srng.GetUniform(0.0, 1.0);
				if (change_prob > real_pi[feature]) {
					feature_changed = true;
					total_features_changed++;

					if (feature == 0) {
						if (RAND) {
							int newwin = sig.win;
							while (newwin == sig.win) newwin = srng.GetUint() % 65536;
							sig.win = newwin;
						}

						if (HALFWAY) {
							int r = 0;
							do {
								r = srng.GetUint() % winvals.size();
							} while (sig.win == winvals[r]);
							sig.win = winvals[r];
						}
						
						if (POPULAR) {
							sig.win = original_sigs.at(pop1)[0].win;
							if (sig.win == it->second[s].win) sig.win = original_sigs[pop2][0].win;
						}
						if (UNIQUE) sig.win = j;
					}
					if (feature == 1) {
						if (RAND || HALFWAY) {
							int newttl = sig.ttl;
							while (newttl == sig.ttl) {
								auto iter = ttl_map.begin();
								advance(iter, srng.GetUint() % ttl_map.size());
								newttl = iter->first;
							}
							sig.ttl = newttl;
						}

						/*
						int index = ttl_map[sig.ttl];
						int new_index = (index + 1) % ttl_map.size();
						auto iter = ttl_map.begin();
						advance(iter, new_index);
						sig.ttl = iter->first;
						*/

						if (POPULAR) {
							sig.ttl = original_sigs[pop1][0].ttl;
							if (sig.ttl == it->second[s].ttl) sig.ttl = original_sigs[pop2][0].ttl;
						}

						if (UNIQUE) {
							auto iter = ttl_map.begin();
							advance(iter, j % ttl_map.size());
							sig.ttl = iter->first;
						}
					}
					if (feature == 2) {
						if (RAND || HALFWAY) sig.df = sig.df ^ 1;

						if (POPULAR) {
							sig.df = original_sigs[pop1][0].df;
							if (sig.df == it->second[s].df) sig.df = original_sigs[pop2][0].df;
						}

						if (UNIQUE) sig.df = j % 2;
					}
					if (feature == 3) {
						if (RAND) {
							int newmss = sig.mss;
							while (newmss == sig.mss) newmss = srng.GetUint() % 65536; //full range possible
							sig.mss = newmss;
						}

						if (HALFWAY) {
							int r = 0;
							do {
								r = srng.GetUint() % mssvals.size();
							} while (sig.mss == mssvals[r]);
							sig.mss = mssvals[r];
						}

						if (POPULAR) {
							sig.mss = original_sigs[pop1][0].mss;
							if (sig.mss = it->second[s].mss) sig.mss = original_sigs[pop2][0].mss;
						}

						if (UNIQUE) sig.mss = j;
					}

					if (feature == 4) {
						if (RAND) {
							//pick a valid subset with prob OPT_PI_SUBSET
							double opt_prob = srng.GetUniform(0.0, 1.0);
							if (opt_prob < OPT_PI_SUBSET && sig.option_order > 0) { //make sure its not empty option set
								//get a random valid subset
								auto iter = option_subsets[sig.option_order].begin();
								advance(iter, srng.GetUint() % option_subsets[sig.option_order].size());

								sig.options_int = *iter;
							}
							else {								
								sig.options_int = generateRandomOptions();
							}
						}

						if (HALFWAY) {
							int r = 0;
							double opt_prob = srng.GetUniform(0.0, 1.0);
							if (opt_prob < OPT_PI_SUBSET) {
								if (option_subsets[sig.option_order].size() > 0) { //only if there is a subset for this option
									do {
										r = srng.GetUint() % optvals.size();
									} while (sig.options_int == optvals[r] || option_subsets[sig.option_order].count(optvals[r]) == 0); //keep going till we pick a valid subset
									sig.options_int = optvals[r];
								}
							}
							else {
								sig.options_int = generateRandomOptions();
							}
						}

						if (POPULAR) {
							sig.options_int = original_sigs[pop1][0].options_int;
							if (sig.options_int == it->second[s].options_int) sig.options_int = original_sigs[pop2][0].options_int;
						}

						if (UNIQUE) {
							if (option_subsets[sig.option_order].size() > 0) {
								auto iter = option_subsets[sig.option_order].begin();							
								advance(iter, j % option_subsets[sig.option_order].size());
								sig.options_int = *iter;
							}
						}

						getOptionOrder(sig.options_int, sig.option_order);
						//if (!containMSS(sig.options_int)) sig.mss = 0; //MSS fix, but causes divergence in extreme cases
					}
					if (feature == 5) {						
						if (RAND) {
							if (sig.rst) {
								//if there was a RST, either lose it or generate random
								if (rst_lost) sig.rst_int = 0;
								else {
									//randomize current rst
									sig.rst_ack = sig.rst_ack ^ 1;
									sig.rst_win = srng.GetUint() % 3;
									sig.rst_nonzero = sig.rst_nonzero ^ 1;

									sig.rst_int = getRSTInt(sig.rst, sig.rst_ack, sig.rst_win, sig.rst_nonzero);
								}
							}
							else {
								//if there wasnt a rst, inject random one
								sig.rst = 1;
								sig.rst_ack = srng.GetUint() % 2;
								sig.rst_win = srng.GetUint() % 3;
								sig.rst_nonzero = srng.GetUint() % 2;

								sig.rst_int = getRSTInt(sig.rst, sig.rst_ack, sig.rst_win, sig.rst_nonzero);
							}
						}

						if (HALFWAY) {
							int r = 0;
							do {
								r = srng.GetUint() % rstvals.size();
							} while (sig.rst_int == rstvals[r]);
							sig.rst_int = rstvals[r];
						}

						if (POPULAR) {
							sig.rst = original_sigs[pop1][0].rst;
							sig.rst_ack = original_sigs[pop1][0].rst_ack;
							sig.rst_win = original_sigs[pop1][0].rst_win;
							sig.rst_nonzero = original_sigs[pop1][0].rst_nonzero;
							sig.rst_int = original_sigs[pop1][0].rst_int;

							if (sig.rst_int == it->second[0].rst_int) {
								sig.rst = original_sigs[pop2][0].rst;
								sig.rst_ack = original_sigs[pop2][0].rst_ack;
								sig.rst_win = original_sigs[pop2][0].rst_win;
								sig.rst_nonzero = original_sigs[pop2][0].rst_nonzero;
								sig.rst_int = original_sigs[pop2][0].rst_int;
							}
						}

						if (UNIQUE) {
							sig.rst = j % 2;
							if (sig.rst) {
								sig.rst_ack = j % 2;
								sig.rst_win = j % 3;
								sig.rst_nonzero = j % 2;
							}
							else {
								sig.rst_ack = 0;
								sig.rst_win = 0;
								sig.rst_nonzero = 0;
							}

							sig.rst_int = getRSTInt(sig.rst, sig.rst_ack, sig.rst_win, sig.rst_nonzero);
						}						
					}
				}			
			}
			total_features += real_pi.size();
			

			//add signature to the final list if it has packets, if all were lost then we cant observe it
			if (!sig.packet_arrival_timestamp.empty()){
				ret_list.push_back(sig);

				//create pmf for real loss
				real_loss_dists[it->second[s].packet_arrival_timestamp.size()][packets_lost]++;

				//add constant values to their maps, if not already there
				//if (win_map.count(sig.win) == 0) win_map[sig.win] = win_map.size();
				//if (mss_map.count(sig.mss) == 0) mss_map[sig.mss] = mss_map.size();
				//if (rst_map.count(sig.rst_int) == 0) rst_map[sig.rst_int] = rst_map.size();

				if (opt_map.count(sig.options_int) == 0) opt_map[sig.options_int] = opt_map.size();
				if (opt_map.count(sig.option_order) == 0) opt_map[sig.option_order] = opt_map.size();

				if (feature_changed) feature_changed_counter++;

				//track real dist of OS 0 for testing against
				if (j == 0) {
					real_win_dist[sig.win]++;
					real_ttl_dist[ttl_map[sig.ttl]]++;
					real_df_dist[sig.df]++;
					real_mss_dist[sig.mss]++;
				}
			}
			else {
				discarded++;
				i--;
			}
		}
	}

	printf("Created %zd signatures, (%.2f%%) had loss, (%.2f%%) had constant feature changed\n", 
		ret_list.size(), ((double)loss_counter / ret_list.size()) * 100, ((double) feature_changed_counter / ret_list.size()) * 100);
	printf("Total packet loss rate = %.6f%%\n", (total_packets_lost / total_packets) * 100);
	printf("Total feature change rate = %.6f%%\n", (total_features_changed / total_features) * 100);
	printf("Discarded samples = %d\n", discarded);

	normalize(real_t_dist);
	normalize(real_q_dist);
	for (int n = 0; n < real_loss_dists.size(); n++) normalize(real_loss_dists[n]);
	normalize(real_win_dist);
	normalize(real_ttl_dist);
	normalize(real_df_dist);
	normalize(real_mss_dist);

	//print original ratios to file
	printf("Outputting real alphas and distributions to files...\n");
	FILE* falpha;
	fopen_s(&falpha, "alphas.txt", "a");
	fprintf(falpha, "\nT mean = %.2f, OWD mean = %.2f, packet loss rate = %.2f%%, feature change rate = %.2f%%\n\n", t_mean, delay_mean, (total_packets_lost / total_packets) * 100, (total_features_changed / total_features) * 100);	
	fprintf(falpha, "0.0000 | ");
	int id = INT_MIN;
	double count = 0;
	for (int i = 0; i < ret_list.size(); i++) {
		if (ret_list[i].ip_int != id) {
			//output count of last OS
			if (i > 0) {
				fprintf(falpha, "%.9lf ", count / ret_list.size());
			}

			//reset counts and id
			count = 1;
			id = ret_list[i].ip_int;
		}
		else count++;
	}
	fprintf(falpha, "%.9lf ", count / ret_list.size());
	fprintf(falpha, "\n--------------\n");
	fclose(falpha);
	
	//output real distributions
	FILE *ft, *fowd, *floss, *fconst;
	fopen_s(&ft, "real_t_dist.txt", "w");
	fopen_s(&fowd, "real_owd_dist.txt", "w");
	fopen_s(&floss, "real_loss_dist.txt", "w");

	for (double d : real_t_dist) fprintf(ft, "%.9g\n", d);
	for (double d : real_q_dist) fprintf(fowd, "%.9g\n", d);
	for (int n = 0; n < real_loss_dists.size(); n++){
		for (double d : real_loss_dists[n]) fprintf(floss, "%.9g ", d);
		fprintf(floss, "\n");
	}


	fopen_s(&fconst, "real_win_dist.txt", "w");
	for (double d : real_win_dist) if (d > 0) fprintf(fconst, "%.9g\n", d);
	fclose(fconst);

	fopen_s(&fconst, "real_ttl_dist.txt", "w");
	for (double d : real_ttl_dist) fprintf(fconst, "%.9g\n", d);
	fclose(fconst);

	fopen_s(&fconst, "real_df_dist.txt", "w");
	for (double d : real_df_dist) fprintf(fconst, "%.9g\n", d);
	fclose(fconst);

	fopen_s(&fconst, "real_mss_dist.txt", "w");
	for (double d : real_mss_dist) fprintf(fconst, "%.9g\n", d);
	fclose(fconst);

	fclose(floss);
	fclose(ft);
	fclose(fowd);

	return ret_list;
}

int main(int argc, char* argv[]){

	time_t start = clock();
	srng.SetState(5, 5); //fix seed for reproducibility
	srand(5);
	
	/***************** Restart previous run *****************************/
	int readIterData = 0; //to read from previous saved iteration data, just set the command line argument to 1
	if (argc > 1) readIterData = atoi(argv[1]);

	/********** Load database signatures *********************************/

	char* database_filename;
	if (argc > 2) database_filename = argv[2];
	else {
		database_filename = "420OS_db.txt";
	}
	int max_rto_length;
	//map<u_int, vector<Signature>> siglist = readSigs("3-signatures.txt", max_rto_length);
	//map<u_int, vector<Signature>> siglist = readSigs("116OS_db.txt", max_rto_length);
	map<u_int, vector<Signature>> siglist = readSigs(database_filename, max_rto_length);

	//read IP to OS name mappings for 420 db.
	if (siglist.size() == 420) IP_to_OS_map = readMapping("420OS_mapping.txt");

	//compact the database by removing similar subOSes from each signature. This is only to help EM run faster
	compactList(siglist, COMPACTION_THRESHOLD);

	/********** Load observations ********************************************/
	//clear alpha output
	FILE* falpha;
	fopen_s(&falpha, "alphas.txt", "w");
	fclose(falpha);

	char* test_filename;
	if (argc > 3) test_filename = argv[3];
	else {
		test_filename = "test_observations.txt";
	}
	vector<Signature> testdata = readTestData(test_filename, max_rto_length, 1.0);

	//override to create simulated data
	//vector<Signature> testdata = createDisturbedZipf(siglist, SIMULATE_T_MEAN, SIMULATE_DELAY_MEAN, SIMULATE_LOSS, max_rto_length, false);
		
	/**************** For testing unknowns ******************************/
	//inject an OS into the database signatures for testing
	//injectUnknown(siglist, .70);	

	//remove the injected OS from DB now that it has been simulated and added to test data
	//siglist.erase(UINT_MAX); 

	//remove a portion of database signatures
	//shrinkDatabase(siglist, disturbed, .50);

	/************** Run Faulds *****************************************/

	//run Faulds given database and testdata
	Faulds(siglist, testdata, max_rto_length, readIterData);
	//HardEM_MT(siglist, disturbed);

	time_t end = clock();
	printf("Time taken: %lf hours\n", ((double)(end - start) / CLOCKS_PER_SEC) / 60 / 60);
}