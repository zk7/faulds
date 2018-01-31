/*
Implementation of Hard EM using Hershel+ as the classifier for a performance comparison 
Author: Zain Shamsi
*/

#include "classes.h"

#define ITERATIONS 50

struct SimpleThreadVars {
	int thread_id;
	HANDLE mutex;
	HANDLE events[2];
	HANDLE iteration_done_event;

	map<u_int, vector<Signature>>* db_signatures;
	vector<Signature>* sample_signatures;
	vector<double>* prob_array;
	vector<double>* alphas;
	map<unsigned int, int> classifications;
	vector<vector<double>>* calculated_probs;
	int iteration;
	int hosts_done;
	int correct;
};

inline size_t key(int i, int j) { return (size_t)i << 32 | (unsigned int)j; }

void HardEM_Thread(LPVOID* tv) {
	SimpleThreadVars* vars = (SimpleThreadVars*)tv;
	map<u_int, vector<Signature>>::iterator osj;

	//give each thread an id
	WaitForSingleObject(vars->mutex, INFINITE);
	int myid = vars->thread_id;
	vars->thread_id++;
	ReleaseMutex(vars->mutex); // RELEASE Mutex

	while (true) {

		if (WaitForMultipleObjects(2, vars->events, FALSE, INFINITE) == WAIT_OBJECT_0) return;

		//go through all samples
		for (int i = myid; i < vars->sample_signatures->size(); i += NUM_THREADS) {

			int osnum = 0; //the os we are on
			double rowsum = 0; //summation of all the p(x|j)*alpha_j for all j
			double highest_prob = 0; //to keep track of highest
			int classification = -1;
			int winning_osindex = -1;

			//for each OS j
			for (osj = vars->db_signatures->begin(); osj != vars->db_signatures->end(); osj++) {

				double prob_i_j = 0;
				if (vars->iteration == 1){		
					//calculate the prob for sample i, OS j
					int lost_packets = osj->second[0].packet_arrival_timestamp.size() - (*vars->sample_signatures)[i].packet_arrival_timestamp.size();

					if (lost_packets >= 0) {
						//get probability of sample i from OS j using H+D
						HershelPlus hpd((*vars->sample_signatures)[i], osj->second);
						prob_i_j = hpd.totalprob;						
					}

					//store it
					(*vars->calculated_probs)[i][osnum] = prob_i_j;
				}
				else {
					//get it from matrix
					prob_i_j = (*vars->calculated_probs)[i][osnum];
				}

				double finalprob = (*vars->alphas)[osnum] * prob_i_j;

				//gather highest classification
				if (finalprob > highest_prob) {
					highest_prob = finalprob;
					classification = osj->first;
					winning_osindex = osnum;
				}			

				osnum++;
			}

			WaitForSingleObject(vars->mutex, INFINITE);
			ResetEvent(vars->events[1]); //reset group_ready event (could cause race if one thread gets here before others enter iteration, but generally should not happen)

			vars->hosts_done++;
			//keep track of accuracy
			if (classification == (*vars->sample_signatures)[i].ip_int) {
				vars->correct++;
			}

			if (winning_osindex >= 0) (*vars->prob_array)[winning_osindex]++;
			vars->classifications[classification]++;

			if (vars->hosts_done == vars->sample_signatures->size()) {
				SetEvent(vars->iteration_done_event); //set group done event
			}
			ReleaseMutex(vars->mutex); // RELEASE Mutex
		}
	}
}

void HardEM_MT(map<u_int, vector<Signature>> original_sigs, vector<Signature> disturbed) {
	FILE* falpha;

	//set up alpha arraySimpleIter_Thread
	vector<double> alphas(original_sigs.size(), 1.0 / original_sigs.size());
	int total_time = 0;

	//set up OWD probabilities
	double sum = 0;
	vector<double> parray;
	for (DWORD s = 0; s < OWD_NUM_BINS; s++) {
		double x = (s + 1) * OWD_GUESS_LIMIT_IN_SECS / OWD_NUM_BINS;
		double prob = exp(-HP_OWD_LAMBDA * x);
		parray.push_back(prob);
		sum += prob;
	}
	//normalize array
	for (int i = 0; i < parray.size(); i++) {
		parray[i] = parray[i] / sum;
	}
	HershelPlus::owd_pmf = parray;

	parray.clear();
	double common = pow(HP_ERLANG_NU, HP_ERLANG_K) / factorial(HP_ERLANG_K - 1);
	for (int s = 0; s < T_NUM_BINS; s++) {
		double x = (s + 1) * T_GUESS_LIMIT_IN_SECS / T_NUM_BINS;
		double prob = common * pow(x, HP_ERLANG_K - 1) * exp(-HP_ERLANG_NU * x);
		parray.push_back(prob);
		sum += prob;
	}
	//normalize array
	for (int i = 0; i < parray.size(); i++) {
		parray[i] = parray[i] / sum;
	}

	HershelPlus::t_pmf = parray;

	//matrix to gather up the new probabilities for p(OSj)
	vector<double> prob_array(original_sigs.size(), 0);

	//matrix to store probabilities for H+
	vector<vector<double>> calculated_probs(disturbed.size(), vector<double>(original_sigs.size(), 0));

	//create common variables
	SimpleThreadVars tv;
	tv.prob_array = &prob_array;
	tv.db_signatures = &original_sigs;
	tv.sample_signatures = &disturbed;
	tv.calculated_probs = &calculated_probs;
	tv.alphas = &alphas;
	tv.iteration = 1;
	tv.thread_id = 0;
	tv.mutex = CreateMutex(NULL, 0, NULL);
	tv.iteration_done_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	tv.events[0] = CreateEvent(NULL, TRUE, FALSE, NULL);
	tv.events[1] = CreateEvent(NULL, TRUE, FALSE, NULL);

	HANDLE *handles = new HANDLE[NUM_THREADS];
	//Split Threads
	for (int i = 0; i < NUM_THREADS; i++) {
		handles[i] = CreateThread(NULL, 4096, (LPTHREAD_START_ROUTINE)HardEM_Thread, &tv, 0, NULL);
		SetThreadPriority(handles[i], THREAD_PRIORITY_LOWEST);
	}
	printf("Started %d threads...\n\n", NUM_THREADS);

	double acc_start, acc_max = 0, acc_last;
	while (tv.iteration <= ITERATIONS) {
		//set up the iteration
		//reset vars, zero out the vectors
		tv.correct = 0;
		tv.hosts_done = 0; 
		std::fill(prob_array.begin(), prob_array.end(), 0);
		tv.classifications.clear();
		
		printf("------------------------------------------\n");

		//set iteration ready to start
		SetEvent(tv.events[1]); //set group ready

		int step = 0;
		int old_done = 0;
		int sleep_interval = 1000;
		while (WaitForSingleObject(tv.iteration_done_event, sleep_interval) == WAIT_TIMEOUT) {
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

		//calculate and update probability for each OS
		for (int j = 0; j < prob_array.size(); j++) {
			prob_array[j] /= disturbed.size();			
			alphas[j] = prob_array[j];
		}

		printf("\naccuracy: %.4lf\n", (double)tv.correct / disturbed.size());

		//output to files
		fopen_s(&falpha, "alphas.txt", "a");
		printf("Iteration: %d\n", tv.iteration);
		double acc = (double)tv.correct / disturbed.size();
		fprintf(falpha, "\n%.4lf | ", acc);
		for (int j = 0; j < alphas.size(); j++) {
			fprintf(falpha, "%.4lf ", alphas[j]);
		}
		printf("\n");
		fclose(falpha);


		if (tv.iteration == 1) acc_start = acc;
		if (acc_max < acc) acc_max = acc;
		acc_last = acc;

		tv.iteration++;

	}

	//Set Quit and wait for threads to return
	printf("\n---Quit Condition Reached!---\nWaiting for all threads to end..\n");
	SetEvent(tv.events[0]); //quit event
	WaitForMultipleObjects(NUM_THREADS, handles, TRUE, INFINITE);

	// Close handles
	for (int i = 0; i < NUM_THREADS; i++) CloseHandle(handles[i]);
	CloseHandle(tv.mutex);

	printf("Accuracy (Start - Max - End): %lf - %lf - %lf\n", acc_start, acc_max, acc_last);
	printf("Time taken for EM iterations: %d seconds\n", total_time);

	_fcloseall();
}
