/*
The Faulds EM algorithm
Author: Zain Shamsi
*/

#include "classes.h"

struct ThreadVars{
	int thread_id;
	HANDLE mutex;
	HANDLE events[2];
	HANDLE iteration_done_event;

	map<u_int, vector<Signature>>* db_signatures;
	vector<Signature>* sample_signatures;
	vector<double>* prob_array;
	vector<vector<double>>* global_new_loss_pmfs;
	vector<double>* global_new_t_pmf;
	vector<double>* global_new_owd_pmf;
	vector<vector<double>>* global_new_windowsize_pmf;
	vector<vector<double>>* global_new_ttl_pmf;
	vector<vector<double>>* global_new_df_pmf;
	vector<vector<double>>* global_new_mss_pmf;
	vector<vector<double>>* global_new_opt_pmf;
	vector<vector<double>>* global_new_rst_pmf;
	vector<double>* rowsum_dist;
	vector<double>* alphas;
	map<unsigned int, int> classifications;
	int hosts_done;
	int correct;
};

void smooth(vector<double>& pmf, double h) {
	//kernel density estimation
	double sum = 0;
	vector<double> new_pmf(pmf.size(), 0);

	//store symmetric kernel densities for h*2
	vector<double> kernel_densities;
	double gfactor = 1.0 / (sqrt(2 * M_PI)); //for gaussian	
	for (int i = -h; i <= h; i++) {
		double u = i / h;
		double k = gfactor * exp(-(u*u) / 2); //gaussian with sigma = 1
		//double k = 0.75 * (1 - (u*u)); //epanechnikov
		kernel_densities.push_back(k);
	}

	double avgfactor = 1.0 / pmf.size() * h;
	for (int x = 0; x < new_pmf.size(); x++) {
		double newval = 0;
		for (int i = -h; i <= h; i++) { //limit range to 
			if (x + i > 0 && x + i < pmf.size()) {
				if (pmf[x + i] > 0) { //if there is a value

					//multiply kernel density at i by value of bin at i and sum
					newval += kernel_densities[i + h] * pmf[x + i];
				}
			}
		}

		//set new val
		new_pmf[x] = newval;
		sum += newval;
	}

	//normalize and set new values
	if (sum > 0) {
		for (int x = 0; x < new_pmf.size(); x++) {
			pmf[x] = new_pmf[x] / sum;
		}
	}

	//spike smoothing
	//double spike_thresh = 0.01;
	//double above_thresh = 0;
	//vector<int> below_thresh;
	//
	//for (int i = 0; i < pmf.size(); i++) {
	//	if (pmf[i] > spike_thresh) above_thresh += pmf[i];
	//	else below_thresh.push_back(i);
	//}
	//
	//for (int i = 0; i < below_thresh.size(); i++)
	//	pmf[below_thresh[i]] = (1 - above_thresh) / below_thresh.size();

	//laplace add-1 smoothing - tends towards uniform
	//double sum = 0;
	//for (int i = 0; i < pmf.size(); i++) sum += pmf[i];
	//
	//for (int i = 0; i < pmf.size(); i++) {
	//	if (pmf[i] == 0) {
	//		pmf[i] = (pmf[i] + 1) / (sum + 65536 * 1);
	//	}
	//}

	//nearest neighbor
	//int h = 10;
	//double neighborsum = 0;
	//int have_before = 0;
	//int have_after = 0;
	//for (int x = 0; x < pmf.size(); x++) {
	//	int val_before = h / 2;
	//	int val_after = h / 2;
	//
	//	//sum values before x
	//	if (have_before < val_before) {
	//		//get as many as possible from before x
	//		for (int b = have_before + 1; have_before < val_before; b++) {
	//			if (x - b < 0) break;
	//			neighborsum += pmf[x - b];
	//			have_before++;
	//		}
	//	}
	//	else {
	//		//drop the last from the window and add previous val
	//		neighborsum -= pmf[x - val_before];
	//		neighborsum += pmf[x - 1];
	//	}
	//
	//	//sum values after x
	//	if (have_after < val_after) {
	//		//get as many as possible from after x
	//		for (int b = have_after + 1; have_after < val_after; b++) {
	//			if (x + b >= pmf.size()) break;
	//			neighborsum += pmf[x + b];
	//			have_after++;
	//		}
	//	}
	//	else {
	//		//drop the current from the window and advance window if possible
	//		neighborsum -= pmf[x];
	//		if (x + val_after < pmf.size())	neighborsum += pmf[x + val_after];
	//	}
	//
	//	pmf[x] = neighborsum / (have_before + have_after);
	//}

}

void writeIterationState(char* filename, int iteration, vector<double>& alphas,
	vector<double>& owd_pmf, vector<double>& t_pmf, vector<vector<double>>& loss_pmfs,
	vector<vector<double>>& win_pmf, vector<vector<double>>& ttl_pmf, vector<vector<double>>& df_pmf,
	vector<vector<double>>& mss_pmf, vector<vector<double>>& opt_pmf, vector<vector<double>>& rst_pmf) {

	//write state to file to restart later from this iteration

	FILE* fout;
	int ret = fopen_s(&fout, filename, "wb");
	if (ret != 0) {
		printf("Error opening file %s!\n", filename);
		return;
	}

	//write iteration
	fwrite(&iteration, sizeof(int), 1, fout);

	//write alphas
	int sz = alphas.size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < sz; i++) fwrite(&(alphas[i]), sizeof(double), 1, fout);
	
	//write OWD
	sz = owd_pmf.size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < sz; i++) fwrite(&(owd_pmf[i]), sizeof(double), 1, fout);
	
	//write T
	sz = t_pmf.size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < sz; i++) fwrite(&(t_pmf[i]), sizeof(double), 1, fout);

	//write loss
	sz = loss_pmfs.size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < loss_pmfs.size(); i++) {
		for (int j = 0; j < i + 1; j++) {
			fwrite(&(loss_pmfs[i][j]), sizeof(double), 1, fout);
		}
	}

	//write constants
	int OSes = win_pmf.size();
	fwrite(&OSes, sizeof(int), 1, fout);

	sz = win_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(win_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	sz = ttl_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(ttl_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	sz = df_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(df_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	sz = mss_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(mss_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	sz = opt_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(opt_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	sz = rst_pmf[0].size();
	fwrite(&sz, sizeof(int), 1, fout);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fwrite(&(rst_pmf[i][j]), sizeof(double), 1, fout);
		}
	}

	fclose(fout);
}

void readIterationState(char* filename, int* iteration, vector<double>& alphas,
	vector<double>& owd_pmf, vector<double>& t_pmf, vector<vector<double>>& loss_pmfs,
	vector<vector<double>>& win_pmf, vector<vector<double>>& ttl_pmf, vector<vector<double>>& df_pmf,
	vector<vector<double>>& mss_pmf, vector<vector<double>>& opt_pmf, vector<vector<double>>& rst_pmf) {

	printf("Reading start state from from %s...\n", filename);

	FILE* fin;
	int ret = fopen_s(&fin, filename, "rb");
	if (ret != 0) {
		printf("Error opening file %s!\n", filename);
		exit(-1);
	}

	fread(iteration, sizeof(int), 1, fin);

	//read alphas
	int sz;
	fread(&sz, sizeof(int), 1, fin);
	if (sz != alphas.size()) {
		printf("Size to read from %s does not match alpha! (%d != %d!)\n", filename, sz, alphas.size());
		exit(-1);
	}
	for (int i = 0; i < sz; i++) fread(&(alphas[i]), sizeof(double), 1, fin);

	//read OWD
	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < sz; i++) fread(&(owd_pmf[i]), sizeof(double), 1, fin);

	//read T
	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < sz; i++) fread(&(t_pmf[i]), sizeof(double), 1, fin);

	//read loss
	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < loss_pmfs.size(); i++) {
		for (int j = 0; j < i + 1; j++) {
			fread(&(loss_pmfs[i][j]), sizeof(double), 1, fin);
		}
	}

	//read constants
	int OSes;
	fread(&OSes, sizeof(int), 1, fin);
	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(win_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(ttl_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(df_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(mss_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(opt_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fread(&sz, sizeof(int), 1, fin);
	for (int i = 0; i < OSes; i++) {
		for (int j = 0; j < sz; j++) {
			fread(&(rst_pmf[i][j]), sizeof(double), 1, fin);
		}
	}

	fclose(fin);
}

void Faulds_Thread(LPVOID* tv) {
	ThreadVars* vars = (ThreadVars*)tv;
	map<u_int, vector<Signature>>::iterator osj;

	//give each thread an id
	WaitForSingleObject(vars->mutex, INFINITE);
	int myid = vars->thread_id;
	vars->thread_id++;
	ReleaseMutex(vars->mutex); // RELEASE Mutex

	//set up temp vector and pmfs
	vector<double> row(vars->db_signatures->size(), 0);
	vector<double> new_owd_pmf(OWD_NUM_BINS, 0);
	vector<double> new_t_pmf(T_NUM_BINS, 0);	
	vector<vector<double>> new_loss_pmfs = (*vars->global_new_loss_pmfs);
	vector<bool> opt_flags(vars->db_signatures->size(), true); //flag saying whether to use full option string (true) or only order (false)

	while (true) {

		if (WaitForMultipleObjects(2, vars->events, FALSE, INFINITE) == WAIT_OBJECT_0) return;

		//go through all samples
		for (int i = myid; i < vars->sample_signatures->size(); i += NUM_THREADS) {
			
			int osnum = 0; //the os we are on
			double rowsum = 0; //summation of all the p(x|j)*alpha_j for all j
			double highest_prob = 0; //to keep track of highest
			int classification = 0;
			int winning_osnum = -1;
			
			fill(row.begin(), row.end(), 0);
			fill(new_t_pmf.begin(), new_t_pmf.end(), 0);
			fill(new_owd_pmf.begin(), new_owd_pmf.end(), 0);
			for (int k = 0; k < new_loss_pmfs.size(); k++) fill(new_loss_pmfs[k].begin(), new_loss_pmfs[k].end(), 0);			

			//for each OS j
			for (osj = vars->db_signatures->begin(); osj != vars->db_signatures->end(); osj++){

				int lost_packets = osj->second[0].packet_arrival_timestamp.size() - (*vars->sample_signatures)[i].packet_arrival_timestamp.size();
				//if (lost_packets >= 0){
				if (lost_packets >= 0 && lost_packets < 10) {
					//get probability of sample i from OS j using H+
					HershelPlusEM hp((*vars->sample_signatures)[i], osj->second, osnum, (*vars->alphas)[osnum], &new_t_pmf, &new_owd_pmf, &opt_flags);					
					double prob_i_j = hp.totalprob;

					//HershelPlus hpold((*vars->sample_signatures)[i], osj->second);
					//double prob_i_j = hpold.totalprob;

					//gather highest classification
					if (prob_i_j > highest_prob) {
						highest_prob = prob_i_j;
						classification = osj->first;
						winning_osnum = osnum;
					}

					//ignore samples with 1 packet for the EM update
					//if ((*vars->sample_signatures)[i].packet_arrival_timestamp.size() <= 1) prob_i_j = 0;

					//update loss probability					
					new_loss_pmfs[osj->second[0].packet_arrival_timestamp.size()][lost_packets] += prob_i_j;

					//update EM matrix
					//(*vars->prob_array)[i][osnum] = prob_i_j; 
					//rowsum += (*vars->prob_array)[i][osnum];					

					//update row
					row[osnum] = prob_i_j;
					rowsum += prob_i_j;
				}

				osnum++;
			}

			//store winning prob
			(*vars->sample_signatures)[i].prob = highest_prob;
			(*vars->sample_signatures)[i].predicted_ip_int = classification;
			(*vars->sample_signatures)[i].predicted_winner_index = winning_osnum;

			WaitForSingleObject(vars->mutex, INFINITE);
			ResetEvent(vars->events[1]); //reset group_ready event (could cause race if one thread gets here before others enter iteration, but generally should not happen)

			//check if there is a valid probability in the row
			if (rowsum > 0){
			//if (rowsum > 1e-6){ //try having a threshold on rowsum
				//now that we have denominator, calculate each value in row i for OS j
				//for (int j = 0; j < osnum; j++) (*vars->prob_array)[i][j] /= rowsum;
				for (int j = 0; j < osnum; j++) {
					row[j] /= rowsum;
					(*vars->prob_array)[j] += row[j];
				}

				//update the global T/OWD pmfs after normalizing the temp pmf
				for (int index = 0; index < new_t_pmf.size(); index++) (*vars->global_new_t_pmf)[index] += new_t_pmf[index] / rowsum;
				//for (int index = 0; index < new_owd_pmf.size(); index++) (*vars->global_new_owd_pmf)[index] += new_owd_pmf[index] / (rowsum * (*vars->sample_signatures)[i].packet_arrival_timestamp.size());
				for (int index = 0; index < new_owd_pmf.size(); index++) (*vars->global_new_owd_pmf)[index] += new_owd_pmf[index] / rowsum; //based on Daren's recommendation					
					
				//normalize(new_t_pmf);
				//normalize(new_owd_pmf);
				//for (int index = 0; index < new_t_pmf.size(); index++) (*vars->global_new_t_pmf)[index] += new_t_pmf[index];
				//for (int index = 0; index < new_owd_pmf.size(); index++) (*vars->global_new_owd_pmf)[index] += new_owd_pmf[index];

				//update the global loss pmfs after normalizing loss pmfs for this sample
				for (int n = 0; n < new_loss_pmfs.size(); n++) {
					for (int k = 0; k < new_loss_pmfs[n].size(); k++) {
						(*vars->global_new_loss_pmfs)[n][k] += new_loss_pmfs[n][k] / rowsum;
					}
				}

				//update the global constants probs
				for (int j = 0; j < osnum; j++) {

					//only update if beyond a threshold
					//if ((*vars->prob_array)[i][j] > UPDATE_THRESHOLD) {

					(*vars->global_new_windowsize_pmf)[j][(*vars->sample_signatures)[i].win] += row[j]; // (*vars->prob_array)[i][j];
					(*vars->global_new_ttl_pmf)[j][ttl_map[(*vars->sample_signatures)[i].ttl]] += row[j]; // (*vars->prob_array)[i][j];
					(*vars->global_new_df_pmf)[j][(*vars->sample_signatures)[i].df] += row[j]; // (*vars->prob_array)[i][j];
					(*vars->global_new_mss_pmf)[j][(*vars->sample_signatures)[i].mss] += row[j]; // (*vars->prob_array)[i][j];

					if (opt_flags[j]) (*vars->global_new_opt_pmf)[j][opt_map[(*vars->sample_signatures)[i].options_int]] += row[j]; // (*vars->prob_array)[i][j];
					else (*vars->global_new_opt_pmf)[j][opt_map[(*vars->sample_signatures)[i].option_order]] += row[j]; //(*vars->prob_array)[i][j];

					(*vars->global_new_rst_pmf)[j][(*vars->sample_signatures)[i].rst_int] += row[j]; //(*vars->prob_array)[i][j];
					//}
				}
				
			}

			vars->hosts_done++;
			//keep track of accuracy
			if (classification == (*vars->sample_signatures)[i].ip_int){
				vars->correct++;
				(*vars->sample_signatures)[i].was_correct = true;
			}

			vars->classifications[classification]++;

			if (vars->hosts_done == vars->sample_signatures->size()){
				SetEvent(vars->iteration_done_event); //set group done event
			}
			ReleaseMutex(vars->mutex); // RELEASE Mutex
		}
	}
}

void Faulds(map<u_int, vector<Signature>> original_sigs, vector<Signature> disturbed, int max_rto_length, bool read_iteration_data){

	//reset files
	FILE *falpha, *fowd, *ft, *fl, *fconstants, *fclasses;
	//fopen_s(&falpha, "alphas.txt", "w");
	fopen_s(&fowd, "owd_dist.txt", "w");
	fopen_s(&ft, "t_dist.txt", "w");
	fopen_s(&fl, "loss_dist.txt", "w");
	fopen_s(&fclasses, "class_dist.txt", "w");

	fclose(fclasses);
	fclose(fl);
	fclose(fowd);
	fclose(ft);

	fopen_s(&fconstants, "winningprob_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "pi.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "windowsize_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "ttl_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "df_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "opt_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "mss_dist.txt", "w");
	fclose(fconstants);
	fopen_s(&fconstants, "rst_dist.txt", "w");
	fclose(fconstants);
	//fclose(falpha);
	
	////////////////////START SETTING INITIAL VALUES/////////////////
	int em_iteration = 0;

	vector<double> alphas(original_sigs.size(), 1.0 / original_sigs.size());
	vector<double> alpha_updates(original_sigs.size());
	double alpha_change = DBL_MAX;
	int total_time = 0;

	vector<double> owd_pmf;
	vector<double> t_pmf;

	//set up OWD probabilities for Hershel+
	double sum = 0;
	for (int s = 0; s < OWD_NUM_BINS; s++) {
		double x = (s + 1) * OWD_GUESS_LIMIT_IN_SECS / OWD_NUM_BINS;
		double prob = exp(-HP_OWD_LAMBDA * x);
		owd_pmf.push_back(prob);
		sum += prob;
	}
	//normalize array
	for (int i = 0; i < owd_pmf.size(); i++){
		owd_pmf[i] = owd_pmf[i] / sum;
	}

	//set up PMFs
	/////////////////////////DL's method//////////////////
	// now the T model: erlang(k); mean = k / lambda 
	//double sum = 0;
	vector<double> parray;
	double common = pow(HP_ERLANG_NU, HP_ERLANG_K) / factorial(HP_ERLANG_K - 1);
	for (int s = 0; s < T_NUM_BINS; s++) {
		double x = (s + 1) * T_GUESS_LIMIT_IN_SECS / T_NUM_BINS;
		double prob = common * pow(x, HP_ERLANG_K - 1) * exp(-HP_ERLANG_NU * x);
		parray.push_back(prob);
		sum += prob;
	}
	//normalize array
	for (int i = 0; i < parray.size(); i++){
		parray[i] = parray[i] / sum;
	}

	//owd_pmf = parray;
	t_pmf = parray;
	
	//vector<double> parray(T_NUM_BINS, 1.0 / T_NUM_BINS);
	//owd_pmf = parray;
	//t_pmf = parray;
	
	HershelPlus::t_pmf = t_pmf;
	HershelPlus::owd_pmf = owd_pmf;
	vector<double> global_new_t_pmf(t_pmf.size(), 0);	
	vector<double> global_new_owd_pmf(owd_pmf.size(), 0);	

	/////////////////////////////////////////////////////////////

	//generate matrix for (n choose k) values
	vector<vector<double>> n_choose_k(max_rto_length + 1, vector<double>());
	for (int n = 0; n < n_choose_k.size(); n++) {
		for (int k = 0; k <= n; k++) {
			double p = (factorial(n) / (factorial(k) * factorial(n - k)));
			n_choose_k[n].push_back(p);
		}
	}

	/////////////////////////////////////////////////////////////
	//Start loss pmfs as binomial
	vector<vector<double>> loss_pmfs(max_rto_length + 1, vector<double>());

	for (int n = 0; n < loss_pmfs.size(); n++){
		//resize each pmf accordingly
		loss_pmfs[n].resize(n + 1);

		//then set up for losing k packets out of n
		for (int k = 0; k <= n; k++){
			//binomial pmf
			double prob_lose_k_n = n_choose_k[n][k] * pow(HP_LOSS_PROB, k) * pow(1 - HP_LOSS_PROB, n - k);			
			loss_pmfs[n][k] = prob_lose_k_n;
		}
	}

	//copy structure to the update pmfs, but reset the update pmf at the beginning of the iteration
	vector<vector<double>> global_new_loss_pmfs = loss_pmfs;

	//output starting loss
	fopen_s(&fl, "loss_dist.txt", "a");
	for (int n = 0; n < global_new_loss_pmfs.size(); n++) {
		for (double d : global_new_loss_pmfs[n]) fprintf(fl, "%.9g\t", d);
		fprintf(fl, "\n");
	}
	fprintf(fl, "\n_________________________________________\n");
	fclose(fl);

	/////////////////////////////////////////////////////////////

	//go through OS database and assign starting probabilities (.9 to known value, .1 / MAX-1 to everything else)
	double prob_win_nochange = 0.9; 
	double prob_ttl_nochange = 0.9;
	double prob_df_nochange = 0.9; 
	double prob_mss_nochange = 0.9;
	double prob_opt_nochange = 0.9;
	double prob_opt_subset = 0.1;
	//double prob_opt_subset = 0.09;
	double prob_rst_nochange = 0.9;
	//vector<vector<double>> windowsize_pmf(original_sigs.size(), vector<double>(65536, (1-prob_win_nochange) / 65535));
	vector<vector<double>> windowsize_pmf(original_sigs.size(), vector<double>(65536, (1 - prob_win_nochange) / (65536 - 1)));
	vector<vector<double>> ttl_pmf(original_sigs.size(), vector<double>(4, (1-prob_ttl_nochange) / 3));
	vector<vector<double>> df_pmf(original_sigs.size(), vector<double>(2, 1-prob_df_nochange));
	vector<vector<double>> mss_pmf(original_sigs.size(), vector<double>(65536, (1 - prob_mss_nochange) / (65536 - 1)));
	vector<vector<double>> opt_pmf(original_sigs.size(), vector<double>(opt_map.size(), 0));
	vector<vector<double>> rst_pmf(original_sigs.size(), vector<double>(32, (1 - prob_rst_nochange) / (32- 1)));
	map<u_int, vector<Signature>>::iterator os_it;
	int osnum = 0;
	for (os_it = original_sigs.begin(); os_it != original_sigs.end(); os_it++, osnum++){		
		int win = os_it->second[0].win;
		int ttl = os_it->second[0].ttl;
		int df = os_it->second[0].df;
		int mss = os_it->second[0].mss;
		long long int opt = os_it->second[0].options_int;
		int rst = os_it->second[0].rst_int;

		windowsize_pmf[osnum][win] = prob_win_nochange;
		ttl_pmf[osnum][ttl_map[ttl]] = prob_ttl_nochange;
		df_pmf[osnum][df] = prob_df_nochange;
		mss_pmf[osnum][mss] = prob_mss_nochange;

		opt_pmf[osnum][opt_map[opt]] = prob_opt_nochange;		
		//set values for subsets in pmf		
		sum = prob_opt_nochange;
		for (int subset : option_subsets[os_it->second[0].option_order]) {			
			opt_pmf[osnum][opt_map[subset]] = prob_opt_subset / option_subsets[os_it->second[0].option_order].size();
			sum += opt_pmf[osnum][opt_map[subset]];
		}
		//go through again and set the remaining values
		//int remaining_bins = opt_pmf[osnum].size() - 1 - option_subsets[os_it->second[0].option_order].size();
		//double prob_remaining = 1 - sum;
		//for (int i = 0; i < opt_pmf[osnum].size(); i++){
		//	if (opt_pmf[osnum][i] == 0) opt_pmf[osnum][i] = prob_remaining / remaining_bins;		
		//}		
		
		rst_pmf[osnum][rst] = prob_rst_nochange;

	}	

	vector<vector<double>> global_new_windowsize_pmf(original_sigs.size(), vector<double>(65536, 0));
	vector<vector<double>> global_new_ttl_pmf(original_sigs.size(), vector<double>(4, 0));
	vector<vector<double>> global_new_df_pmf(original_sigs.size(), vector<double>(2, 0));
	vector<vector<double>> global_new_mss_pmf(original_sigs.size(), vector<double>(65536, 0));
	vector<vector<double>> global_new_opt_pmf(original_sigs.size(), vector<double>(opt_map.size(), 0));
	vector<vector<double>> global_new_rst_pmf(original_sigs.size(), vector<double>(32, 0));
	
	/////////////////////////////////////////////////////////////

	if (read_iteration_data) {
		//read iteration from disk
		readIterationState("iteration_state.bin", &em_iteration, alphas,
			owd_pmf, t_pmf, loss_pmfs,
			windowsize_pmf, ttl_pmf, df_pmf, mss_pmf, opt_pmf, rst_pmf);
	}

	/////////////////////////////////////////////////////////////

	//matrix to gather up all p(xi|j)*alpha_j for all j
	//to optimize space, dont really need to store the entire matrix, just column sums
	//vector<vector<double>> prob_array(disturbed.size(), vector<double>(original_sigs.size(), 0)); 
	vector<double> prob_array(original_sigs.size(), 0);

	//create common variables
	ThreadVars tv;
	tv.db_signatures = &original_sigs;
	tv.sample_signatures = &disturbed;
	tv.prob_array = &prob_array;
	tv.global_new_owd_pmf = &global_new_owd_pmf;
	tv.global_new_t_pmf = &global_new_t_pmf;
	tv.global_new_loss_pmfs = &global_new_loss_pmfs;
	tv.global_new_windowsize_pmf = &global_new_windowsize_pmf;
	tv.global_new_ttl_pmf = &global_new_ttl_pmf;
	tv.global_new_df_pmf = &global_new_df_pmf;
	tv.global_new_mss_pmf = &global_new_mss_pmf;
	tv.global_new_opt_pmf = &global_new_opt_pmf;
	tv.global_new_rst_pmf = &global_new_rst_pmf;
	tv.alphas = &alphas;
	tv.thread_id = 0;
	tv.mutex = CreateMutex(NULL, 0, NULL);
	tv.iteration_done_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	tv.events[0] = CreateEvent(NULL, TRUE, FALSE, NULL);
	tv.events[1] = CreateEvent(NULL, TRUE, FALSE, NULL);

	HANDLE *handles = new HANDLE[NUM_THREADS];
	//Split Threads
	for (int i = 0; i < NUM_THREADS; i++){
		handles[i] = CreateThread(NULL, 4096, (LPTHREAD_START_ROUTINE)Faulds_Thread, &tv, 0, NULL);
		SetThreadPriority(handles[i], THREAD_PRIORITY_LOWEST);
	}
	printf("Started %d threads...\n\n", NUM_THREADS);

	//set up static lookup tables
	HershelPlusEM::owd_pmf = &owd_pmf;
	HershelPlusEM::t_pmf = &t_pmf;
	HershelPlusEM::loss_pmfs = &loss_pmfs;
	HershelPlusEM::n_choose_k = &n_choose_k;
	HershelPlusEM::windowsize_pmf = &windowsize_pmf;
	HershelPlusEM::ttl_pmf = &ttl_pmf;
	HershelPlusEM::df_pmf = &df_pmf;
	HershelPlusEM::mss_pmf = &mss_pmf;
	HershelPlusEM::opt_pmf = &opt_pmf;
	HershelPlusEM::rst_pmf = &rst_pmf;

	double acc_start, acc_max = 0, acc_last;
	em_iteration++;

	while (em_iteration <= EM_MAX_ITERATIONS){	
		//set up the iteration
		//reset vars, zero out the vectors
		tv.correct = 0;
		tv.hosts_done = 0;
		tv.classifications.clear();
		std::fill(prob_array.begin(), prob_array.end(), 0);
		std::fill(alpha_updates.begin(), alpha_updates.end(), 0);
		std::fill(global_new_t_pmf.begin(), global_new_t_pmf.end(), 0);
		std::fill(global_new_owd_pmf.begin(), global_new_owd_pmf.end(), 0);		
		for (int k = 0; k < global_new_loss_pmfs.size(); k++) std::fill(global_new_loss_pmfs[k].begin(), global_new_loss_pmfs[k].end(), 0);
		for (int j = 0; j < global_new_windowsize_pmf.size(); j++) std::fill(global_new_windowsize_pmf[j].begin(), global_new_windowsize_pmf[j].end(), 0);
		for (int j = 0; j < global_new_ttl_pmf.size(); j++) std::fill(global_new_ttl_pmf[j].begin(), global_new_ttl_pmf[j].end(), 0);
		for (int j = 0; j < global_new_df_pmf.size(); j++) std::fill(global_new_df_pmf[j].begin(), global_new_df_pmf[j].end(), 0);
		for (int j = 0; j < global_new_mss_pmf.size(); j++) std::fill(global_new_mss_pmf[j].begin(), global_new_mss_pmf[j].end(), 0);
		for (int j = 0; j < global_new_opt_pmf.size(); j++) std::fill(global_new_opt_pmf[j].begin(), global_new_opt_pmf[j].end(), 0);
		for (int j = 0; j < global_new_rst_pmf.size(); j++) std::fill(global_new_rst_pmf[j].begin(), global_new_rst_pmf[j].end(), 0);

		printf("------------------------------------------\n");

		//set iteration ready to start
		SetEvent(tv.events[1]); //set group ready

		int step = 0;
		int old_done = 0;
		int sleep_interval = 1000;
		while (WaitForSingleObject(tv.iteration_done_event, sleep_interval) == WAIT_TIMEOUT){
			step++;
			int remaining_hosts = tv.sample_signatures->size() - tv.hosts_done;
			double rate = tv.hosts_done / (double)(step * (double)(sleep_interval / 1000));
			double remaining_time = (remaining_hosts / rate) / 60;

			if (tv.hosts_done > old_done) //only print if we had some successes
				printf("IPs left: %d, Success: %d at %.3f per sec. %.2f min to go\r", remaining_hosts, tv.hosts_done, rate, remaining_time);

			old_done = tv.hosts_done;
		}	

		int seconds_for_iteration = step * (sleep_interval / 1000);
		total_time += seconds_for_iteration;
		printf("\nDone with iteration in about %d seconds\n", seconds_for_iteration);

		///////////////////////----------------ITERATION PROCESSING-----------//////////////////////////

		//once we've built the table, we can calculate the updates
		//add up all the values across the columns in alpha_updates
		//for (int i = 0; i < prob_array.size(); i++){
			//for (int j = 0; j < prob_array[i].size(); j++){
				//alpha_updates[j] += prob_array[i][j];
			//}
		//}
		sum = 0;
		for (int j = 0; j < prob_array.size(); j++) {
			alpha_updates[j] += prob_array[j];
			sum += prob_array[j];
		}
		

		//once done adding up all values average the new alphas and calculate change
		alpha_change = 0;
		for (int j = 0; j < alpha_updates.size(); j++){
			//alpha_updates[j] /= disturbed.size();
			alpha_updates[j] /= sum;
			alpha_change += fabs(alpha_updates[j] - alphas[j]);
			alphas[j] = alpha_updates[j];
		}
		
		//output to files

		int injected_count = 0;
		//keep track of injected samples
		//unordered_map<u_int, int> injected_class;
		//fopen_s(&fconstants, "injected.txt", "a");
		//fopen_s(&fclasses, "winningprob_dist.txt", "a");
		//for (int d = 0; d < disturbed.size(); d++) {
			//if (original_sigs.count(disturbed[d].ip_int) == 0) {
				//injected_count++;
				//injected_class[disturbed[d].predicted_ip_int]++;		
				//fprintf(fconstants, "%e\n", disturbed[d].prob);
			//}
		//	else fprintf(fclasses, "%e\n", disturbed[d].prob);
		//}
		//printf("Injected Count: %d\n", injected_count);
		//fclose(fconstants);
		//fclose(fclasses);
		

		fopen_s(&falpha, "alphas.txt", "a");
		printf("Iteration: %d\n", em_iteration);
		double acc = (double)tv.correct / (disturbed.size() - injected_count);		
		fprintf(falpha, "\n%.4lf | ", acc);
		printf("\nalpha change: %lf, accuracy: %.4lf (ignore acc if real data)\n", alpha_change, acc);

		map<int, double> osnum_to_alpha; //store a mapping from os index to alpha
		for (int j = 0; j < alphas.size(); j++){
			fprintf(falpha, "%.9lf ", alphas[j]);
			osnum_to_alpha[j] = alphas[j];
		}
		printf("\n");


		//normalize new distributions
		normalize(global_new_t_pmf);
		normalize(global_new_owd_pmf);		
		for (int n = 0; n < global_new_loss_pmfs.size(); n++) normalize(global_new_loss_pmfs[n]);			
			
		for (int n = 0; n < global_new_windowsize_pmf.size(); n++) normalize(global_new_windowsize_pmf[n]);
		for (int n = 0; n < global_new_ttl_pmf.size(); n++) normalize(global_new_ttl_pmf[n]);
		for (int n = 0; n < global_new_df_pmf.size(); n++) normalize(global_new_df_pmf[n]);
		for (int n = 0; n < global_new_mss_pmf.size(); n++) normalize(global_new_mss_pmf[n]);
		for (int n = 0; n < global_new_opt_pmf.size(); n++) normalize(global_new_opt_pmf[n]);
		for (int n = 0; n < global_new_rst_pmf.size(); n++) normalize(global_new_rst_pmf[n]);		

		//output constant feature change
		//printf("constant feature change: ");
		fopen_s(&fconstants, "pi.txt", "a");
		double totalnochange = 0;
		double totalwin_nochange = 0, totalttl_nochange = 0, totaldf_nochange = 0, totalmss_nochange = 0, totalopt_nochange = 0, totalrst_nochange = 0;		
		int osnum = 0;
		for (os_it = original_sigs.begin(); os_it != original_sigs.end(); os_it++, osnum++) {
			int win = os_it->second[0].win;
			int ttl = os_it->second[0].ttl;
			int df = os_it->second[0].df;
			int mss = os_it->second[0].mss;
			long long int opt = os_it->second[0].options_int;
			int rst = os_it->second[0].rst_int;
			
			double win_no_change = global_new_windowsize_pmf[osnum][win];
			double ttl_no_change = global_new_ttl_pmf[osnum][ttl_map[ttl]];
			double df_no_change = global_new_df_pmf[osnum][df];
			double mss_no_change = global_new_mss_pmf[osnum][mss];
			double opt_no_change = global_new_opt_pmf[osnum][opt_map[opt]];
			double rst_no_change = global_new_rst_pmf[osnum][rst];

			//printf("(%.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf) ", (1 - win_no_change), (1 - ttl_no_change), (1-df_no_change), (1-mss_no_change), (1-opt_no_change), (1-rst_no_change));			
			if (acc < 0.01 && IP_to_OS_map.count(os_it->first) > 0) fprintf(fconstants, "%s, %d, %g, %d, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", IP_to_OS_map[os_it->first].c_str(), osnum, osnum_to_alpha[osnum], tv.classifications[os_it->first], (win_no_change), (ttl_no_change), (df_no_change), (mss_no_change), (opt_no_change), (rst_no_change));
			else fprintf(fconstants, "%d, %g, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf\n", osnum, osnum_to_alpha[osnum], (win_no_change), (ttl_no_change), (df_no_change), (mss_no_change), (opt_no_change), (rst_no_change));

			totalnochange += osnum_to_alpha[osnum] * (win_no_change + ttl_no_change + df_no_change + mss_no_change + opt_no_change + rst_no_change) / 6;
			totalwin_nochange += osnum_to_alpha[osnum] * win_no_change;
			totalttl_nochange += osnum_to_alpha[osnum] * ttl_no_change;
			totaldf_nochange += osnum_to_alpha[osnum] * df_no_change;
			totalmss_nochange += osnum_to_alpha[osnum] * mss_no_change;
			totalopt_nochange += osnum_to_alpha[osnum] * opt_no_change;
			totalrst_nochange += osnum_to_alpha[osnum] * rst_no_change;
		}
		fclose(falpha);
		fprintf(fconstants, "\n\n");
		fclose(fconstants);
		printf("Average feature change prob: %.4lf (%.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf) \n", (totalnochange),
			(totalwin_nochange), (totalttl_nochange), (totaldf_nochange), (totalmss_nochange), (totalopt_nochange), (totalrst_nochange));

		//average ploss for 3 - 10 packets, take k=0 case for easy calculation
		//if y = (n k) p^k (1-p)^n-k, then for k = 0: y = (1-p)^n, and  p = 1 - y^(1/n)	
		if (original_sigs.size() > 10) {
			/*
			double totalploss = 0;
			for (int n = 3; n <= 10; n++) {
				double ploss = 1 - pow(global_new_loss_pmfs[n][0], 1.0 / n);
				totalploss += ploss;
			}
			totalploss /= 8;
			printf("Average loss prob using old avg: %lf\n", totalploss);
			*/

			double numerator = 0;
			double denominator = 0;
			for (int k = 3; k < 10; k++) {
				double plossk = 0;
				for (int l = 0; l < global_new_loss_pmfs[k].size(); l++) {
					plossk += l * global_new_loss_pmfs[k][l];
				}

				for (os_it = original_sigs.begin(), osnum = 0; os_it != original_sigs.end(); os_it++, osnum++) {
					if (os_it->second[0].packet_arrival_timestamp.size() == k) {
						numerator += osnum_to_alpha[osnum] * plossk;
						denominator += osnum_to_alpha[osnum] * k;
					}
				}				
			}
			double ploss = numerator / denominator;

			printf("new loss prob: %lf\n", ploss);
		}
		
		/*
		if (em_iteration == 1) {
			int zp_count = 0;
			fopen_s(&fclasses, "first_classification.txt", "w");
			for (int i = 0; i < disturbed.size(); i++) {
				if (disturbed[i].predicted_ip_int == 0) zp_count++;
				else {
					fprintf(fclasses, "%u,%d,%s,%d,%d,%d,%s,%d,%d,%d,%d,%d",
						disturbed[i].ip_int, disturbed[i].predicted_winner_index, IP_to_OS_map[disturbed[i].predicted_ip_int].c_str(), disturbed[i].win, disturbed[i].ttl, disturbed[i].df, disturbed[i].options, disturbed[i].mss, disturbed[i].rst, disturbed[i].rst_ack, disturbed[i].rst_win, disturbed[i].rst_nonzero);
					for (double d : disturbed[i].packet_arrival_timestamp) fprintf(fclasses, ",%lf", d);
					fprintf(fclasses, "\n");
				}
			}
			fclose(fclasses);
			printf("Zero probability: %d\n", zp_count);
		}
		*/

		//apply smoothing to constant feature pmfs
		//for (int n = 0; n < global_new_windowsize_pmf.size(); n++) {
		//	printf("Smoothing windowsize pmf for OS %d\r", n);
		//	smooth(global_new_windowsize_pmf[n], 20); //pmf to smooth, bandwidth of kernel
		//}
		//for (int n = 0; n < global_new_mss_pmf.size(); n++) {
		//	printf("Smoothing MSS pmf for OS %d\r", n);
		//	smooth(global_new_mss_pmf[n], 20);
		//}

		//output distributions
		printf("\nOutputting files...\n");
		fopen_s(&fowd, "owd_dist.txt", "a");		
		for (double d : global_new_owd_pmf) fprintf(fowd, "%.9g\t", d);
		fprintf(fowd, "\n");
		fclose(fowd);

		fopen_s(&ft, "t_dist.txt", "a");
		for (double d : global_new_t_pmf) fprintf(ft, "%.9g\t", d);
		fprintf(ft, "\n");
		fclose(ft);

		fopen_s(&fl, "loss_dist.txt", "a");
		for (int n = 0; n < global_new_loss_pmfs.size(); n++){
			for (double d : global_new_loss_pmfs[n]) fprintf(fl, "%.9g\t", d);
			fprintf(fl, "\n");
		}
		fprintf(fl, "\n_________________________________________\n");
		fclose(fl);

		//write iteration to disk
		writeIterationState("iteration_state.bin", em_iteration, alphas,
			global_new_owd_pmf, global_new_t_pmf, global_new_loss_pmfs,
			global_new_windowsize_pmf, global_new_ttl_pmf, global_new_df_pmf,
			global_new_mss_pmf, global_new_opt_pmf, global_new_rst_pmf);

		owd_pmf = global_new_owd_pmf;
		t_pmf = global_new_t_pmf;
		loss_pmfs = global_new_loss_pmfs;
		windowsize_pmf = global_new_windowsize_pmf;
		ttl_pmf = global_new_ttl_pmf;
		df_pmf = global_new_df_pmf;
		mss_pmf = global_new_mss_pmf;
		opt_pmf = global_new_opt_pmf;
		rst_pmf = global_new_rst_pmf;

		if (em_iteration == 0) acc_start = acc;
		if (acc_max < acc) acc_max = acc;
		acc_last = acc;

		em_iteration++;
				
	}

	//Set Quit and wait for threads to return
	printf("\n---Quit Condition Reached!---\nWaiting for all threads to end..\n");
	SetEvent(tv.events[0]); //quit event
	WaitForMultipleObjects(NUM_THREADS, handles, TRUE, INFINITE);

	// Close handles
	for (int i = 0; i < NUM_THREADS; i++) CloseHandle(handles[i]);
	CloseHandle(tv.mutex);


	//output files
	fopen_s(&fconstants, "winningprob_dist.txt", "a");
	fopen_s(&fclasses, "final_classification.txt", "w");
	int zp_count = 0;
	for (int i = 0; i < disturbed.size(); i++) {
		fprintf(fconstants, "%e\n", disturbed[i].prob);
		
		if (disturbed[i].predicted_ip_int == 0) zp_count++;
		else {
			fprintf(fclasses, "%u,%d,%s,%d,%d,%d,%s,%d,%d,%d,%d,%d",
				disturbed[i].ip_int, disturbed[i].predicted_winner_index, IP_to_OS_map[disturbed[i].predicted_ip_int].c_str(), disturbed[i].win, disturbed[i].ttl, disturbed[i].df, disturbed[i].options, disturbed[i].mss, disturbed[i].rst, disturbed[i].rst_ack, disturbed[i].rst_win, disturbed[i].rst_nonzero);
			for (double d : disturbed[i].packet_arrival_timestamp) fprintf(fclasses, ",%lf", d);
			fprintf(fclasses, "\n");
		}
			
	}
	//fclose(fl);
	fclose(fclasses);
	fclose(fconstants);
	printf("Zero probability: %d\n", zp_count);
	
	int n = 0;
	fopen_s(&fconstants, "windowsize_dist.txt", "a");
	for (int n = 0; n < global_new_windowsize_pmf.size(); n++) {
		for (int w = 0; w < global_new_windowsize_pmf[n].size(); w++) {
			//if (global_new_windowsize_pmf[n][w] > 0)
			fprintf(fconstants, "%.4g\t", global_new_windowsize_pmf[n][w]);
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);


	fopen_s(&fconstants, "mss_dist.txt", "a");
	for (int n = 0; n < global_new_mss_pmf.size(); n++) {
		for (int w = 0; w < global_new_mss_pmf[n].size(); w++) {
			fprintf(fconstants, "%.4g\t", global_new_mss_pmf[n][w]);
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);


	fopen_s(&fconstants, "ttl_dist.txt", "a");
	for (int n = 0; n < global_new_ttl_pmf.size(); n++) {
		for (int w = 0; w < global_new_ttl_pmf[n].size(); w++) {
			fprintf(fconstants, "%.4g\t", global_new_ttl_pmf[n][w]);
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);

	fopen_s(&fconstants, "df_dist.txt", "a");
	for (int n = 0; n < global_new_df_pmf.size(); n++) {
		for (int w = 0; w < global_new_df_pmf[n].size(); w++) {
			fprintf(fconstants, "%.4g\t", global_new_df_pmf[n][w]);
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);
	

	fopen_s(&fconstants, "opt_dist.txt", "a");
	unordered_map<int, long long int> reversemap; //build reverse optmap to iterate through pmf
	for (auto p : opt_map) {
		reversemap[p.second] = p.first;
	}
	for (int n = 0; n < global_new_opt_pmf.size(); n++) {
		for (int w = 0; w < global_new_opt_pmf[n].size(); w++) {
			//only print if there is a value > 0
			if (global_new_opt_pmf[n][w] > 0) {
				char optstring[20];
				getOptionString(reversemap[w], optstring, 20);
				fprintf(fconstants, "%s\t%.4g\t", optstring, global_new_opt_pmf[n][w]);
			}
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);
	
	fopen_s(&fconstants, "rst_dist.txt", "a");
	for (int n = 0; n < global_new_rst_pmf.size(); n++) {
		for (int w = 0; w < global_new_rst_pmf[n].size(); w++) {
			//only print the non-zero probs for opt
			if (global_new_rst_pmf[n][w] > 0) {
				int rst, ra, rw, rn;
				getRSTbits(w, rst, ra, rw, rn);
				fprintf(fconstants, "%d%d%d%d\t%.4g\t", rst, ra, rw, rn, global_new_rst_pmf[n][w]);
			}
		}
		fprintf(fconstants, "\n");
	}
	fclose(fconstants);

	
	fopen_s(&fclasses, "class_dist.txt", "a");
	for (auto p : tv.classifications) {
		if (p.first == -1) fprintf(fclasses, "Unknown, %d\n", p.second);
		else fprintf(fclasses, "%s, %d\n", IP_to_OS_map[p.first].c_str(), p.second);
	}
	fprintf(fclasses, "\n---------------------\n");
	fclose(fclasses);


	printf("Accuracy (Start - Max - End): %lf - %lf - %lf\n", acc_start, acc_max, acc_last);
	printf("Time taken for EM iterations: %d seconds\n", total_time);

}
