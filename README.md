# Faulds
Faulds is an iterative OS fingeprinting algorithm for measurement over large scans. It uses the Expectation-Maximization (EM) framework along with an existing classifier (Hershel+ in this case) to iteratively keep re-classifying the same dataset. Using this procedure it can rebuild the probability distributions of the unknown parameters. In the problem of single-probe OS fingerprinting, Faulds recovers empirical distribution of network delays, packet loss, and modification of TCP header fields by end-users seen in the captured data.

**NOTE:** Faulds relies on having a large enough dataset for EM to work and can only improve accuracy over previous methods if it has enough data to learn empirical distributions. If you are working with a small set of observations, or are just interested in a live tool for OS fingerprinting remote hosts, please see [Hershel+](https://github.com/zk7/hershelplus).

# Project files

Faulds is written in C++ compiled on Windows using VS2015.

The Visual Studio project files use the Windows library for multi-threading and hence compile in Win32/64. It also includes all the databases used to test Faulds (3 signature toy database, 116 OS Hershel DB, 420 OS Hershel+ DB). Other files included are winvals, mssvals, optvals, rstvals -- which contain probabilities of values observed for these header fields in our Internet scan. These files are used to generate simulated samples.
 
Running: `Faulds.exe <restart_previous_run> <database_file> <observations_file>` 

All the command line args are optional, it will default to start a new run with the files present in the folder. Simulation parameters and number of EM iterations have to be changed in the code. 

### Generated files

Running Faulds will generate several files. If running with simulated samples, the generation code will output several files prefixed with `real_`. These contain the "real" distributions of the simulated dataset for ground truth. It will also write the ground truth alphas to `alphas.txt`.
The Faulds algorithm itself will output files for each parameter (e.g., `owd_dist.txt` contains recovered one-way delay distributions). Each row of the file corresponds to one iteration. Distributions for user features and final classification results are output after the last iteration. 

### OS Signatures

The data files containing the OS and Internet signatures have mostly the same text format. For the files in the multi-platform folder, this is format:

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

420OS_db.txt contains the database signatures, 420OS_mapping.txt maps plain text labels to signature data in the database, and observations.txt contains sample observed signatures that are to be classified using the Hershel+ algorithm. 


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



