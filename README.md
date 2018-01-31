# Faulds
Faulds is an iterative OS fingerprinting algorithm for measurement over large scans. It uses the Expectation-Maximization (EM) framework along with an existing classifier (Hershel+ in this case) to iteratively keep re-classifying the same dataset while learning the underlying parameters of the data. In the problem of single-probe OS fingerprinting, Faulds recovers empirical distributions of network delays, packet loss, and modification of TCP header fields by end-users seen in the scan.

**NOTE:** Faulds relies on having a large enough dataset for EM to work and can only improve accuracy over previous methods if it has enough data to learn empirical distributions. If you are working with a small set of observations, or are just interested in a live tool for OS fingerprinting remote hosts, please see [Hershel+](https://github.com/zk7/hershelplus).

# Project files

Faulds is written in C++ compiled on Windows using VS2015.

The Visual Studio project files use the Windows library for multi-threading and hence compile in Win32/64. It also includes all the databases used to test Faulds (3 signature toy DB from the paper, 116 OS Hershel DB, 420 OS Hershel+ DB). Other files included are winvals, mssvals, optvals, rstvals -- which contain probabilities of values observed for these header fields in our Internet scan. These files are used to generate simulated samples.
 
Running: `Faulds.exe <restart_previous_run> <database_file> <observations_file>` 

All the commandline args are optional, it will default to start a new run with the files present in the folder. Simulation parameters and number of EM iterations have to be changed in the code. 

### Generated files

Running Faulds will generate several files. If running with simulated samples, the generation code will output several files prefixed with `real_`. These contain the "real" distributions of the simulated dataset for ground truth. It will also write the ground truth alphas to `alphas.txt`.
The Faulds algorithm itself will output files for each parameter (e.g., `owd_dist.txt` contains recovered one-way delay distributions). Each row of the file corresponds to one iteration. Distributions for user features and final classification results are output after the last iteration. 

After each iteration, Faulds writes the currently recovered distribution to `iteration_state.bin`. This can be used as a checkpoint to restart Faulds from this point (the same dataset and database must be provided). 

### OS Signatures

The data files containing the OS and Internet signatures have the same text format:

    int id
    int tcp_window
    int ip_ttl
    int ip_df
    string tcp_options
    longlong tcp_options_encoded
    int mss
    int rst_present
    int rst_ack flag
    int rst_window
    int rst_sequence
    int rst_nonzero
    double RTT (0 value for database)
    double RTO1_timestamp
    double RTO2_timestamp
    double RTO3_timestamp
    ...

`3-signatures.txt`, `116OS_db.txt`, and `420OS_db.txt` contain the database signatures, `420OS_mapping.txt` maps plain text labels to signature data in the database for the 420 database, and `test_observations.txt` contains sample observed signatures for testing. Note that Faulds does not work optimally with this small number of samples. 


# Publication
### Conference
Z. Shamsi, D.B.H. Cline, and D. Loguinov, "Faulds: A Non-Parametric Iterative Classifier for Internet-Wide OS Fingerprinting", ACM CCS, Nov 2017.

    @inproceedings{shamsi2017,
      author = {Zain Shamsi and Daren B.H. Cline and Dmitri Loguinov},
      title = {Faulds: A Non-Parametric Iterative Classifier for Internet-Wide OS Fingerprinting},
      booktitle = {ACM CCS},
      year = {2017},	
      organization = {ACM},
      location = {Dallas, Texas, USA},
      pages = {971--982},
      numpages = {12},
      doi = {10.1145/3133956.3133963},
    }
  
[ACM Portal](https://dl.acm.org/citation.cfm?id=3133956.3133963) 

[Direct Paper Link](http://irl.cs.tamu.edu/people/zain/papers/ccs2017.pdf)



