Simera v2.1
Ben Nichols 2016
bnichols1979@hotmail.com

*** INTRODUCTION ***

The Simera software was written as part of my Ph.D. project at the University of Glasgow. This version executes the Simera 2 algorithm which simulates PCR and has been shown to produce more realistic chimeras than other PCR simulation software (http://theses.gla.ac.uk/6801/). The algorithm is described in the accompanying file simera_2_algorithm.pdf.


*** DISCLAIMER ***

The software has only been tested on Mac OS X (version 10.9.4 and older). I believe that it has been tested sufficiently and that it is free of errors but this is not guaranteed.


*** INSTALLATION ***

First it is required that the gsl library (http://www.gnu.org/software/gsl/) is installed. Once this has been done, install Simera from this directory using the commands:

	make clean
	make


*** INSTRUCTIONS ***

Usage: Simera [OPTIONS]

Required options
  -i FILENAME	Input file in FASTA format.

Additional options
  -o DIRNAME	Directory for output files. If directory doesn't exist then it will be created. Default 'simera_output'.
  -n INT	Number of rounds of PCR to be simulated. Default 25.
  -s INT	Overall sequence abundance to be sampled. Default 10000.
  -l REAL	Value of the parameter lambda. Default 0.00005.
  -c INT	Number of potential chimeras to generate. Default 10000.
  -f STRING	Choice of forward primer. Default GTGNCAGCMGCCGCGGTAA (515f 16S forward primer).
  -r STRING	Choice of reverse primer. Default GGACTACHVGGGTWTCTAAT (806r 16S reverse primer).
  -x		Extra output files (see OUTPUT FILES) generated if this option is selected.
  -h		Displays help information.


*** NOTES ***

The sequence abundance information must be contained in the input FASTA file as part of the sequence name. The abundance should be a positive integer at the end of each name and separated from the rest of the name by the '_' character. The two scripts, convert_abundance.pl and revert_abundance.pl may be used to convert from and to UCHIME (http://drive5.com/usearch/manual/uchime_algo.html) format respectively.

To produce sensible results the input FASTA file must contain, from each sequence, only the region that is to be amplified. Primer Prospector (http://pprospector.sourceforge.net/) is one application that can be used to select the amplicons from the full sequences. The sequence parts matching the forward and reverse primers must then be re-added to the resulting FASTA file of amplicons which can be done using the substr.pl script provided.

The sample size should be changed depending on how many reads are expected from the sequencing platform used. The default of 10000 reads was typical for the 454 pyrosequencing data used for testing.

The default value for lambda (recommended) was found to best reproduce the number of chimeras found in experimental datasets.

During testing the default value for the number of potential chimeras (option -c) was found to work well, with no noticable difference in results for higher values. However, in terms of results, a higher value will never be worse than a lower one (but may cause running time to be much longer).

New primers must be specified if the default primers are inappropriate for your data.


*** OUTPUT FILES ***

- info.txt		Input information.
- all_seqs.fa		All output sequences in FASTA format: full set.
- samp_all_seqs.fa	All output sequences in FASTA format: sampled set.
- good.fa		All good sequences in FASTA format: full set.
- samp_good.fa		All good sequences in FASTA format: sampled set.
- chimeras.fa		All chimeras in FASTA format: full set.
- samp_chimeras.fa	All chimeras in FASTA format: sampled set.
- abund.txt		Output sequence abundances: full set.
- samp_abund.txt	Output sequence abundances: sampled set.
- breaks.txt		Chimera break points: full set.
- samp_breaks.txt	Chimera break points: sampled set.
- parents.fa		Chimera parents in FASTA format: full set.
- samp_parents.txt	Chimera parents in FASTA format: sampled set.
- summary.txt		Summary of output: full set.
- samp_summary.txt	Summary of output: sampled set.

All of these files will be written to the chosen output directory if option -x is included. If option -x is not included then only info.txt and and samp_all_seqs.fa will be written.

All output FASTA files are formatted such that a sequence's name indicates if the sequence is a chimera or not. A sequence's name also contains the sequence's abundance. If a sequence is a chimera then its name will begin with either 'chimeraFWD' or 'chimeraREV' depending on whether the chimera was formed using a forward or reverse primer. The abundance will be given at the end of the sequence name as a positive integer separated from the rest of the name by the '_' character. The output FASTA files may be converted to UCHIME format using the supplied revert_abundance.pl script.


*** SCRIPTS ***

- convert_abundance.pl
Converts a FASTA file with the abundances represented in UCHIME format to the format required for use with Simera.
Usage: convert_abundance.pl input.fa > output.fa

- revert_abundance.pl
Converts a FASTA file with the abundances represented in the format required for use with Simera to UCHIME format.
Usage: revert_abundance.pl input.fa > output.fa

- substr.pl
Used to add the sequence parts matching the forward and reverse primers to the FASTA file of amplicons. The script takes in two FASTA files - one with the full sequences and one with the amplicons output from e.g. Primer Prospector - and two integers X and Y, where X is the length (number of bases) of the forward primer and Y is the length of the reverse primer. The output is the amplicon data in FASTA format with X and Y base pairs of the full sequence (the parts matching the primers) re-attached to the beginning and end of each sequence as required for Simera input.
Usage: substr.pl amplicons.fa sequences.fa X Y > output.fa


*** TUTORIAL ***

This tutorial uses the files present in the tutorial directory and assumes that Simera and all scripts can be executed from within that directory.

Primer Prospector must be installed to complete all steps of this tutorial.

The output files for intermediate steps are included if you want to skip a step.

The FASTA file, example_16S.fa, is a subset of 1000 sequences taken from the Silva 16S dataset (http://www.arb-silva.de/). LogNormal(0,1) random variables, rounded up to the nearest integer, have been used to add abundances to the end of the sequence names.

- Generate 515f_806r_amplicons.fasta from example_16S.fa and primers.txt using Primer Prospector.

	analyze_primers.py -f example_16S.fa -P primers.txt
	get_amplicons_and_reads.py -f example_16S.fa -i 515f_example_16S_hits.txt:806r_example_16S_hits.txt

- Generate example_input.fa using substr.pl. The lengths of the forward and reverse primers are 19 and 20 bps respectively.

	substr.pl 515f_806r_amplicons.fasta example_16S.fa 19 20 > example_input.fa

- Run Simera on example_input.fa.

	Simera -i example_input.fa -o example_output -x

- This will run Simera using example_input.fa as input and will write all files detailed in the OUTPUT FILES section to the example_output directory. Note that this is the same as running the more verbose command below in which all default options have been explicitly included.

	Simera -i example_input.fa -o example_output -n 25 -s 10000 -l 0.00005 -c 10000 -f GTGNCAGCMGCCGCGGTAA -r GGACTACHVGGGTWTCTAAT -x

Following the above steps will result in a simulation of 25 rounds of PCR on a sample containing 1000 different species of bacteria with varying abundances for which the 16S gene has been targeted for amplification. The main output file is example_output/samp_all_seqs.fa which will contain 10000 reads, i.e. the total abundance of all good sequences and chimeras will be 10000.
