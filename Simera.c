/* This software is copyright © by Ben Nichols 2016. */
/* Permission is granted for anyone to copy, use or modify this software for the purposes of research or education provided that this copyright notice is retained and note is given of any modifications. */
/* Contact: bnichols1979@hotmail.com */

/****System includes****/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>
#include "Simera.h"
#include <iostream>

#include "stocc/stocc.h" 

#include "randomc/mersenne.cpp"           // code for random number generator
#include "stocc/stoc1.cpp"                // random library source code
#include "stocc/stoc3.cpp"                // random library source code
#include "stocc/fnchyppr.cpp"             // calculate probabilities of Fisher's distribution
#include "stocc/wnchyppr.cpp"             // calculate probabilities of Wallenius distribution
#include "randomc/userintf.cpp"			  // define system specific user interface

/*global constants*/

static int  nLines = 3;

static char szSequence[] = "ACGTNacgtn-";

using namespace std;

int main(int argc, char *argv[]) {
	
	char* szVersion = "Simera_v2.1";
	
	t_Data tSeqData, tPrimerData, tOutData, tChimData; //Data objects to store sequence info
			                   
	int64_t nPrimAbund;					//Abundance of primers
	char *szInput = NULL;				//Input file must be specified
	char *dirName = "simera_output";	//Default output directory
	int nRounds = 25;					//Default number of rounds
	int nSamp = 10000;					//Default sample size 
	double dLambda = 0.00005;			//Default lambda value
	
	char *PrimerFWD = "GTGNCAGCMGCCGCGGTAA";				//forward primer - default 515f
	char *PrimerREV = Rev_Comp("GGACTACHVGGGTWTCTAAT");		//reverse primer (Simera requires reverse comp of actual reverse primer) - default 806r
	
	int nTotal	= 10000;				//Default number of potential chimeras to generate
	int c;
	int xflag = 0;
	int nDir;
	
	printf("\n========================================\n");
	printf("  Simera - simulates chimera formation\n");
	printf("  Use option -h for help\n");
	printf("========================================\n\n");
	
	while ((c = getopt (argc, argv, "hxi:o:n:s:l:c:f:r:")) != -1)	//Read in command line options
		switch (c)
	{
		case 'h':
			printf("Usage: Simera [OPTIONS]\n\n");
			printf("Required options\n");
			printf("  -i FILENAME\tInput file in FASTA format.\n\n");
			printf("Additional options\n");
			printf("  -o DIRNAME\tDirectory for output files. If directory doesnt exist then it will be created. Default 'simera_output'.\n");
			printf("  -n INT\tNumber of rounds of PCR to be simulated. Default 25.\n");
			printf("  -s INT\tOverall sequence abundance to be sampled. Default 10000.\n");
			printf("  -l REAL\tValue of the parameter lambda. Default 0.00005.\n");
			printf("  -c INT\tNumber of potential chimeras to generate. Default 10000.\n");
			printf("  -f STRING\tChoice of forward primer. Default 'GTGNCAGCMGCCGCGGTAA' (515f 16S forward primer).\n");
			printf("  -r STRING\tChoice of reverse primer. Default 'GGACTACHVGGGTWTCTAAT' (806r 16S reverse primer).\n");
			printf("  -x\t\tExtra output files (see readme.txt) generated if this option is selected.\n");
			printf("  -h\t\tDisplays help information.\n\n");
			return 1;
		case 'x':
			xflag = 1;
			break;
		case 'i':
			szInput = optarg;
			break;
		case 'o':
			dirName = optarg;
			break;
		case 'n':
			nRounds = atoi(optarg);
			break;
		case 's':
			nSamp = atoi(optarg);
			break;
		case 'l':
			dLambda = atof(optarg);
			break;
		case 'c':
			nTotal = atoi(optarg);
			break;
		case 'f':
			PrimerFWD = optarg;
			break;
		case 'r':
			PrimerREV = Rev_Comp(optarg);
			break;
		case '?':
			if (optopt == 'i' || optopt == 'o' || optopt == 'n' || optopt == 's' || optopt == 'l' || optopt == 'c' || optopt == 'f' || optopt == 'r'){
				fprintf (stderr, "Error! Option -%c requires an argument.\n", optopt);
			}
			else if (isprint (optopt)){
				fprintf (stderr, "Error! Unknown option '-%c'.\n", optopt);
			}
			else{
				fprintf (stderr, "Error! Unknown option character '\\x%x'.\n", optopt);
			}
			return 1;
		default:
			abort ();
	}
	
	if(szInput == NULL){
		fprintf (stderr, "Error! Option -i missing.\n");
		return 1;
	}
	
	if(PrimerFWD == NULL){
		fprintf (stderr, "Error! Option -f missing.\n");
		return 1;
	}
	
	if(PrimerREV == NULL){
		fprintf (stderr, "Error! Option -r missing.\n");
		return 1;
	}
	
	if(nRounds < 1){
		fprintf (stderr, "Error! Option -n must be greater than zero.\n");
		return 1;
	}
	
	if(nSamp < 1){
		fprintf (stderr, "Error! Option -s must be greater than zero.\n");
		return 1;
	}
	
	if(nTotal < 1){
		fprintf (stderr, "Error! Option -c must be greater than zero.\n");
		return 1;
	}
	
	if(dLambda > 1 || dLambda < 0){
		fprintf (stderr, "Error! Option -l must be a value between zero and one.\n");
		return 1;
	}
	
	int nSeed = time(0);	//Arbitrary variable RNG seed
	//int nSeed = 10;		//Constant RNG seed for testing
	
	const gsl_rng_type *T;	//RNG set up
	gsl_rng *rng;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	rng = gsl_rng_alloc(T); 
	gsl_rng_set (rng, nSeed);
	
	readData(szInput, &tSeqData);	//Sequences from fasta file
	
	initData(&tPrimerData);			//Primers
	initData(&tOutData);			//Output data (after sampling)
	initData(&tChimData);			//Chimeras
	
	nPrimAbund = tSeqData.nTotAbund*((int64_t)pow(2,nRounds));	//Set primer abundance high enough to not run out
	
	addSeq(PrimerFWD, (char *)"", (char *)"", nPrimAbund/2, 0, 0, 0, TRUE, 0, 0, 0, &tPrimerData);	//Add FWD primers
	addSeq(PrimerREV, (char *)"", (char *)"", nPrimAbund/2, 0, 0, 0, TRUE, 0, 0, 0, &tPrimerData);	//Add REV primers
	
	PCR(rng, &tSeqData, &tPrimerData, &tChimData, dLambda, nPrimAbund, nRounds, nSeed, nTotal);		//Simulate PCR amplification

	nDir = mkdir(dirName,S_IRWXU);		//Create directory for output
	nDir = chdir(dirName);
	writeInputInfo(nSeed, szInput, szVersion, nRounds, dLambda, nSamp, PrimerFWD, PrimerREV, "info.txt");		//Write output data files
	if(xflag == 1){
		
		writeData(&tSeqData,"all_seqs.fa");
		
		writeGood(&tSeqData,"good.fa");
		
		writeChim(&tSeqData,"chimeras.fa");
		
		writeAbund(&tSeqData,"abund.txt");
		
		writeBP(&tSeqData,"breaks.txt");
		
		writeParents(&tSeqData,"parents.fa");
		
		writeSummary(&tSeqData,"summary.txt");
		
	}
	
	Sample(rng, &tSeqData, &tOutData, nSamp);	//Sample from output
	writeData(&tOutData,"samp_all_seqs.fa");	//Write sampled output data files
	
	
	if(xflag == 1){
		
		writeGood(&tOutData,"samp_good.fa");
		
		writeChim(&tOutData,"samp_chimeras.fa");
		
		writeAbund(&tOutData,"samp_abund.txt");
		
		writeBP(&tOutData,"samp_breaks.txt");
		
		writeParents(&tOutData,"samp_parents.fa");
		
		writeSummary(&tOutData,"samp_summary.txt");
		
	}
	
	freeData(&tOutData);
	freeData(&tSeqData);
	freeData(&tPrimerData);
	freeData(&tChimData);
	
	gsl_rng_free(rng);

}


void PCR(gsl_rng *ptRNG, t_Data *ptData, t_Data *ptPrimer, t_Data *ptChimeras, double dLambda, int64_t nPrimAbund, int nRounds, int nSeed, int nTotal){
	
	int i = 0, j = 0, k = 0, m = 0, f = 0;
	
	char *sFrag;
	char szChimera[MAX_LINE_LENGTH];
	char *SeqA;
	char *SeqB;
	
	int nLenA, nLenB, Posit, maxi, nDiffer, stop, n, nDenom, nProp, nExists, nFragWarn, nTotFragLen;
	int64_t nSamp, nSeqAbund, nFragAbund, nFragCount, nNewFragAbund, nSamp_temp, nAbund_temp, nExTot;
	double dWeight, dProb, dAbundA, dAbundB, dMeanWeight, dDiffer, dWeightSum, dPrimerWeightSum, dPrimerMeanWeight, dEff, dFail, dBeta, dMeanFragLen;
	
	int nLenPrFWD = ptPrimer->anLen[0];
	int nLenPrREV = ptPrimer->anLen[1];
	
	nFragAbund = 0;
	dWeightSum = 0;
	dPrimerWeightSum = 0;
	nFragCount = 0;
	nFragWarn = 0;
	nExTot = 0;
	nTotFragLen = 0;
	
	for(i = 0; i < ptData->nSeq; i++){											//Loop through all sequences 
		
		SeqA = strdup(ptData->aszSeq[i]);										//Set sequence to work with
		nLenA = ptData->anLen[i];												//Length of sequence
		dAbundA = ((double)ptData->anAbund[i])/((double)ptData->nTotAbund);		//Abundance of sequence as a proportion of total abundance
		
		dDiffer = SimsPrimerFWD(ptPrimer->aszSeq[0], SeqA);						//Similarities between primer and sequence at optimal primer position
		dDiffer = (double)nLenPrFWD - dDiffer;									//Convert similarities to differences
		dWeight = exp(-dDiffer);												//Weight for forward primer on current sequence
		//dWeight = 0.1832834;//test
		dPrimerWeightSum = dPrimerWeightSum + (dAbundA*dWeight);				//Increase proportional weight sum based on abundances and individual weights
																				//Note that, upon completion, proportional weight sum is equivalent to mean primer weight
		
		dDiffer = SimsPrimerREV(ptPrimer->aszSeq[1], SeqA);						//Repeat for reverse primer
		dDiffer = (double)nLenPrREV - dDiffer;
		dWeight = exp(-dDiffer);
		//dWeight = 0.1832834;//test
		dPrimerWeightSum = dPrimerWeightSum + (dAbundA*dWeight);
		
	
	}
	dPrimerMeanWeight = dPrimerWeightSum*0.5;		//Mean primer weight is halved because there are two primers
	//dPrimerMeanWeight = 1; //test
	
	while (ptChimeras->nSeq < nTotal){				//Algorithm part 1: Generating nTotal potential chimeras
		nSamp = 0;
		i = -1;
		nDenom = ptData->nTotAbund;
	
		while(nSamp == 0){							//Select a sequence A at random with probability proportional to its abundance
		
			i++;
			dProb = (double)ptData->anAbund[i]/(double)nDenom;
			nDenom = nDenom - ptData->anAbund[i];
			nSamp = gsl_ran_binomial(ptRNG, dProb, 1);
		
		}
		
		SeqA = strdup(ptData->aszSeq[i]);										//Set sequence A to work with
		nLenA = ptData->anLen[i];												//Length of sequence A
		dAbundA = ((double)ptData->anAbund[i])/((double)ptData->nTotAbund);		//Proportional abundance of sequence A
		
		j=0;
		while(j >= nLenA || j < nLenPrFWD){				//Find break point at random position j on sequence A using lambda parameter
														//If no break point occurs (amplification completes) then repeat from start until break point occurs
			j = gsl_ran_geometric(ptRNG, dLambda);
		}
		
		nSamp = 0;
		k = -1;
		nDenom = ptData->nTotAbund;
		
		while(nSamp == 0){								//Select a second sequence B at random with probability proportional to its abundance
		
			k++;
			dProb = (double)ptData->anAbund[k]/(double)nDenom;
			nDenom = nDenom - ptData->anAbund[k];
			nSamp = gsl_ran_binomial(ptRNG, dProb, 1);
		
		}
		sFrag = substring(SeqA, 0, j);											//Generate fragment from sequence A at break point j
		
		SeqB = strdup(ptData->aszSeq[k]);										//Set sequence B to work with
		nLenB = ptData->anLen[k];												//Length of sequence B
		dAbundB = ((double)ptData->anAbund[k])/((double)ptData->nTotAbund);		//Proportional abundance of sequence B
					
		SimsFWD(sFrag, SeqB, &Posit, &maxi, nLenPrFWD);							//Find optimal position and number of similarities between fragment and sequence B
		nDiffer = nLenPrFWD - maxi;												//Convert similarities to differences at optimal position
		dWeight = exp(-(double)nDiffer);										//Find weight for fragment when used with sequence B
		dProb = dAbundA*dAbundB*pow(1 - dLambda, j-nLenPrFWD)*dLambda*dWeight;	//Find probability of chimera occurring based on the above weight, lambda and sequence A and B abundances 
		dWeightSum = dWeightSum + ((double)ptData->anAbund[i]*dWeight);
		nFragCount = nFragCount + ptData->anAbund[i];
		strcpy(szChimera,sFrag);
		strcat(szChimera, SeqB + Posit);										//Join fragment and part of sequence to form chimera
		
		addSeq2(szChimera,sFrag,(char *)"chimeraFWD",1, 1, ptData->anID[i], ptData->anID[k], FALSE, j, dWeight, dProb, ptChimeras); //Add chimera and info to pool of chimeras
		nTotFragLen = nTotFragLen + j;
		f++;
		
		free(sFrag);
		
		if(ptChimeras->nSeq < nTotal){					//Repeat to generate reverse fragments and chimeras in the same way
			nSamp = 0;
			i = -1;
			nDenom = ptData->nTotAbund;
		
			while(nSamp == 0){
			
				i++;
				dProb = (double)ptData->anAbund[i]/(double)nDenom;
				nDenom = nDenom - ptData->anAbund[i];
				nSamp = gsl_ran_binomial(ptRNG, dProb, 1);
			
			}
			SeqA = strdup(ptData->aszSeq[i]);
			nLenA = ptData->anLen[i];
			dAbundA = ((double)ptData->anAbund[i])/((double)ptData->nTotAbund);

			j=0;
			while(j >= nLenA || j < nLenPrREV){				//Find break point at random position j on sequence A using lambda parameter
															//If no break point occurs (amplification completes) then repeat from start until break point occurs
				j = gsl_ran_geometric(ptRNG, dLambda);
			}
			
			nSamp = 0;
			k = -1;
			nDenom = ptData->nTotAbund;
		
			while(nSamp == 0){
			
				k++;
				dProb = (double)ptData->anAbund[k]/(double)nDenom;
				nDenom = nDenom - ptData->anAbund[k];
				nSamp = gsl_ran_binomial(ptRNG, dProb, 1);
			
			}
			sFrag = substring(SeqA, nLenA-j-1, nLenA-1);
		
			SeqB = strdup(ptData->aszSeq[k]);
			nLenB = ptData->anLen[k];
			dAbundB = ((double)ptData->anAbund[k])/((double)ptData->nTotAbund);
				
			SimsREV(sFrag, SeqB, &Posit, &maxi, nLenPrREV);
			nDiffer = nLenPrREV - maxi;
			dWeight = exp(-(double)nDiffer);
			dProb = dAbundA*dAbundB*pow(1 - dLambda, j-nLenPrREV)*dLambda*dWeight;
			dWeightSum = dWeightSum + ((double)ptData->anAbund[i]*dWeight);
			nFragCount = nFragCount + ptData->anAbund[i];
			strncpy (szChimera, SeqB, Posit);
			szChimera[Posit] = '\0';
			strcat(szChimera, sFrag);
		
			addSeq2(szChimera,sFrag,(char *)"chimeraREV",1, 1, ptData->anID[i], ptData->anID[k], FALSE, Posit, dWeight, dProb, ptChimeras);
			nTotFragLen = nTotFragLen + j;
			f++;			
			
			free(sFrag);
		}
		
		printf("\rGenerated %i potential chimeras...", ptChimeras->nSeq);
		fflush(stdout);
	}
	
	printf("\n");
	
	dMeanFragLen = (double)nTotFragLen/(double)f;	//Mean fragment length calculated from the above loop (algorithm part 1)
	dMeanWeight = dWeightSum/((double)nFragCount);	//Mean fragment weighting calculated from the above
	nNewFragAbund = 0;
	
	for(i = 0; i < nRounds; i++){					//Algorithm part 2: Simulating sequence amplification and chimera selection - loop for nRounds
		printf("\rPCR Round %i...", i+1);
		fflush(stdout);
		int64_t anIncSeq[ptData->nSeq];				//Array of sequence abundances
		int64_t anIncSeq_temp[ptData->nSeq];
		int64_t anIncChim[ptChimeras->nSeq];		//Array of chimera abundances
		
		nSeqAbund = ptData->nTotAbund;
		nPrimAbund = nSeqAbund;
		
		dEff = ((double)nPrimAbund + (double)nFragAbund)/((double)nFragAbund + (double)nPrimAbund + (double)nSeqAbund);		//PCR efficiency
		dFail = 1 - dEff;																									//PCR failure rate
		nSamp = gsl_ran_binomial(ptRNG, dFail, nSeqAbund);																	//Number of failures
		
		nSeqAbund = (nSeqAbund*2) - nSamp;					
		nPrimAbund = nSeqAbund;
		
		for(j = 0; j < ptData->nSeq; j++){
			dProb = 1 - pow(1 - dLambda, ptData->anLen[j] - nLenPrFWD);		//Abundance of fragments generated from forward primers
			nSamp = gsl_ran_binomial(ptRNG, dProb, ptData->anAbund[j]);
			nNewFragAbund = nNewFragAbund+nSamp;
			
			dProb = 1 - pow(1 - dLambda, ptData->anLen[j] - nLenPrREV);		//Abundance of fragments generated from reverse primers
			nSamp = gsl_ran_binomial(ptRNG, dProb, ptData->anAbund[j]);
			nNewFragAbund = nNewFragAbund+nSamp;
			
		}
		
		nSeqAbund = nSeqAbund - nNewFragAbund;								//Reduce amplified sequences (fragmented sequences don't amplify)
		
		if(nSeqAbund < 1){
			fprintf(stderr, "\nError! PCR failure. Chosen value for lambda parameter may be too high.\n"); //If no sequences amplify then PCR fails
			exit(EXIT_FAILURE);
		}
		
		dBeta = (dPrimerMeanWeight*(double)nPrimAbund)/((dPrimerMeanWeight*(double)nPrimAbund) + (dMeanWeight*(double)nFragAbund));	//Beta determines probability of picking a primer instead of a fragment	
		
		StochasticLib1 sto((int64_t)nSeed);
		
		stop = 0;
		nAbund_temp = nSeqAbund;
		nSamp = 0;
		
		while(stop == 0){								//Use Beta to choose how many primers are used for amplification, the remainder are fragments 
			
			if(nAbund_temp > 300000000000){
				nSamp_temp = sto.Binomial(300000000000, dBeta);
				nSamp = nSamp + nSamp_temp;
				nAbund_temp = nAbund_temp - 300000000000;
			}
			else {
				nSamp_temp = sto.Binomial(nAbund_temp, dBeta);
				nSamp = nSamp + nSamp_temp;
				stop = 1;
			}
			
		}
		
		nPrimAbund = nPrimAbund - nSamp;
		
		for(j = 0; j < ptData->nSeq; j++){
			incSeq(j, ptData->anAbund[j], ptData);	//Double abundance
		}
		
		
		stop = 0;
		nSamp_temp = nSamp;
		for(j = 0; j < ptData->nSeq; j++){
			anIncSeq[j] = 0;
		}
		
		while(stop == 0){							//Use hypergeometric RV to determine which sequences are amplified using new abundance
			
			if(nSamp_temp > 1000000000){
				sto.MultiHypergeometric(anIncSeq_temp, ptData->anAbund, 1000000000, ptData->nSeq);
				for(j = 0; j < ptData->nSeq; j++){
					anIncSeq[j] = anIncSeq[j] + anIncSeq_temp[j];
				}
				nSamp_temp = nSamp_temp - 1000000000;
			}
			else {
				sto.MultiHypergeometric(anIncSeq_temp, ptData->anAbund, nSamp_temp, ptData->nSeq);
				for(j = 0; j < ptData->nSeq; j++){
					anIncSeq[j] = anIncSeq[j] + anIncSeq_temp[j];
				}
				stop = 1;
			}
			
		}
		
		emptyData(ptData);						//Reduce all sequence abundances to 0
		
		for(j = 0; j < ptData->nSeq; j++){
			incSeq(j, anIncSeq[j], ptData);		//Add new abundances generated from hypergeom
		}
		
		nSeqAbund = nSeqAbund - nSamp;			//Remaining sequences combine with fragments (SeqAbund == Chimera Abund)
		nFragAbund = nFragAbund - nSeqAbund;
		nAbund_temp = nSeqAbund;
		
		
		if(nSeqAbund > 0){
		
			for(j = 0; j < ptData->nSeq; j++){																	//Abundance of fragments generated from other fragments
				nProp = (int)(((double)ptData->anAbund[j]/(double)ptData->nTotAbund)*(double)nAbund_temp);
				dProb = 1 - pow(1 - dLambda, max(((double)ptData->anLen[j] - dMeanFragLen),0.0));				//Probability of fragmentation when fragment is used instead of primer
				nSamp = gsl_ran_binomial(ptRNG, dProb, nProp);
				nNewFragAbund = nNewFragAbund+nSamp;
				nSeqAbund = nSeqAbund - nSamp;
			}
			
			sto.Multinomial(anIncChim, ptChimeras->adProb, nSeqAbund, ptChimeras->nSeq);		//Select chimeras based on their probabilities of forming
			
		
			for(j = 0; j < ptChimeras->nSeq; j++){												//Add chimeras and info to sequence pool
				if(anIncChim[j] > 0){
					addSeq(ptChimeras->aszSeq[j], ptChimeras->aszFrag[j], ptChimeras->aszID[j], anIncChim[j], ptChimeras->anNmera[j], ptChimeras->anParentA[j], ptChimeras->anParentB[j], 0, ptChimeras->anBreak[j], ptChimeras->adWeight[j], ptChimeras->adProb[j], ptData);
					nExists = existSeq(ptChimeras->aszSeq[j], ptData);
					if(nExists == TRUE) nExTot = nExTot + anIncChim[j];
				}
			}
		}
		
		nFragAbund = nNewFragAbund + nFragAbund;
		if(nFragAbund < 0){
			nFragAbund = 0;
			nFragWarn = 1;
		}
		nNewFragAbund = 0;
	}
	printf("\nDone.\n");
	if(nFragWarn == 1) printf("Warning! Fragment:Primer usage was unexpectedly high. Please check that you are using appropriate primers for your data.\n");
}

char* substring(char* str, size_t begin, size_t len)	//Create substring of length len starting from position begin 
{
	char *sub = (char*) malloc(len+1);
	strncpy(sub, str+begin, len);
	sub[len] = '\0';
	
	return sub;
}

void initData(t_Data *ptData)	//New (empty) t_Data object
{
	
	ptData->nSeq		= 0;
	ptData->nSize		= INIT_SIZE;
	ptData->nMaxLen		= 0;
	ptData->nTotAbund	= 0;
	ptData->nMaxAbund	= 0;
	ptData->nNextID		= 1;
	
	ptData->aszID		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->aszSeq		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->aszFrag		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->anLen       = (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anAbund		= (int64_t *) malloc(ptData->nSize*sizeof(int64_t));
	ptData->anNmera		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anID		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anParentA	= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anParentB	= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anBreak		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->adWeight	= (double *) malloc(ptData->nSize*sizeof(double));
	ptData->adProb		= (double *) malloc(ptData->nSize*sizeof(double));
	
}

void freeData(t_Data *ptData)	//Free t_Data object
{	
	int i;
	for(i = 0; i < ptData->nSeq; i++){
		free(ptData->aszID[i]);
		free(ptData->aszSeq[i]);
	}

	free(ptData->aszID);
	free(ptData->aszSeq);
	free(ptData->anLen);
	free(ptData->anAbund);
	free(ptData->anNmera);
	free(ptData->anID);
	free(ptData->anParentA);
	free(ptData->anParentB);
	free(ptData->anBreak);
}
		
void incSeq(int id, int64_t inc, t_Data *ptData)		//Increase abundance of sequence id by inc amount
{
	ptData->nTotAbund = ptData->nTotAbund + inc;
	ptData->anAbund[id] = ptData->anAbund[id] + inc;
	if(ptData->anAbund[id] > ptData->nMaxAbund){
		
		ptData->nMaxAbund = ptData->anAbund[id];
	}
}

int existSeq(char* seq, t_Data *ptData)		//Check if sequence already exists
{
	int i;
	int nSeqLength = strlen(seq);
	
	for(i = 0; i < ptData->nSeq; i++){
    	
		if(nSeqLength == ptData->anLen[i] && strncmp(seq, ptData->aszSeq[i], nSeqLength) == 0){
			
			return TRUE;	//return if exists
		}
		
	}
	
	return FALSE;
}

void addSeq(char* seq, char* fragment, char* name, int64_t abund, int chimera, int parentA, int parentB, int frag, int bp, double weight, double prob, t_Data *ptData)	//Add abund copies of sequence seq to set ptData
{
	int i;
	int nSeqLength = strlen(seq);
	ptData->nTotAbund = ptData->nTotAbund + abund;	//increase abundance of all sequences combined
	
	for(i = 0; i < ptData->nSeq; i++){
    	
		if(nSeqLength == ptData->anLen[i] && strncmp(seq, ptData->aszSeq[i], nSeqLength) == 0){	//check if sequence already exists
			ptData->anAbund[i] = ptData->anAbund[i] + abund;	//increase abundance of sequence
			
			if(ptData->anAbund[i] > ptData->nMaxAbund){
				
				ptData->nMaxAbund = ptData->anAbund[i];
			}
			return;//return if exists
		}
		
	}
	
	//for new sequences
	
	ptData->nSeq++;	//increase number of sequences by one
	
	if(ptData->nSeq == ptData->nSize){	//reallocate memory if there is too much data
		
		ptData->nSize *= 2;
		ptData->aszID		= (char **)realloc(ptData->aszID, ptData->nSize*sizeof(char *));
		ptData->aszSeq		= (char **)realloc(ptData->aszSeq, ptData->nSize*sizeof(char *));
		ptData->aszFrag		= (char **)realloc(ptData->aszFrag, ptData->nSize*sizeof(char *));
		ptData->anLen		= (int *)realloc(ptData->anLen, ptData->nSize*sizeof(int));
		ptData->anAbund		= (int64_t *)realloc(ptData->anAbund, ptData->nSize*sizeof(int64_t));
		ptData->anNmera		= (int *)realloc(ptData->anNmera, ptData->nSize*sizeof(int));
		ptData->anID		= (int *)realloc(ptData->anID, ptData->nSize*sizeof(int));
		ptData->anParentA	= (int *)realloc(ptData->anParentA, ptData->nSize*sizeof(int));
		ptData->anParentB	= (int *)realloc(ptData->anParentB, ptData->nSize*sizeof(int));
		ptData->anBreak		= (int *)realloc(ptData->anBreak, ptData->nSize*sizeof(int));
		ptData->adWeight	= (double *)realloc(ptData->adWeight, ptData->nSize*sizeof(double));
		ptData->adProb		= (double *)realloc(ptData->adProb, ptData->nSize*sizeof(double));
		
	}
	
	if(frag == FALSE){
		
		ptData->anID[ptData->nSeq - 1] = ptData->nNextID;
		ptData->nNextID++;
		
	}
	else{
		
		ptData->anID[ptData->nSeq - 1] = 0;
		
	}
	
	ptData->aszSeq[ptData->nSeq - 1] = strdup(seq);	//add info for new sequence
	ptData->aszFrag[ptData->nSeq - 1] = strdup(fragment);
	ptData->aszID[ptData->nSeq - 1] = strdup(name);
	ptData->anLen[ptData->nSeq - 1] = strlen(seq);
	ptData->anAbund[ptData->nSeq - 1] = abund;
	ptData->anNmera[ptData->nSeq - 1] = chimera;
	ptData->anParentA[ptData->nSeq - 1] = parentA;
	ptData->anParentB[ptData->nSeq - 1] = parentB;
	ptData->anBreak[ptData->nSeq - 1] = bp;
	ptData->adWeight[ptData->nSeq - 1] = weight;
	ptData->adProb[ptData->nSeq - 1] = prob;
	
	if(strlen(seq) > ptData->nMaxLen){
		
		ptData->nMaxLen = strlen(seq);	//if new sequence is longest, change the maximum sequence length record
		
	}
	if(abund > ptData->nMaxAbund){
		
		ptData->nMaxAbund = abund;
		
	}
}

void addSeq2(char* seq, char* fragment, char* name, int64_t abund, int chimera, int parentA, int parentB, int frag, int bp, double weight, double prob, t_Data *ptData)	//Add abund copies of sequence seq to set ptData
{
	int i;
	int nSeqLength = strlen(seq);
	ptData->nTotAbund = ptData->nTotAbund + abund;	//increase abundance of all sequences combined
	
	for(i = 0; i < ptData->nSeq; i++){
    	
		if(nSeqLength == ptData->anLen[i] && strncmp(seq, ptData->aszSeq[i], nSeqLength) == 0 && parentA == ptData->anParentA[i] && parentB == ptData->anParentB[i] && bp == ptData->anBreak[i]){	//check if sequence already exists
			ptData->anAbund[i] = ptData->anAbund[i] + abund;	//increase abundance of sequence
			
			if(ptData->anAbund[i] > ptData->nMaxAbund){
				
				ptData->nMaxAbund = ptData->anAbund[i];
			}
			return;//return if exists
		}
		
	}
	
	//for new sequences
	
	ptData->nSeq++;	//increase number of sequences by one
		
	if(ptData->nSeq == ptData->nSize){	//reallocate memory if there is too much data
		
		ptData->nSize *= 2;
		ptData->aszID		= (char **)realloc(ptData->aszID, ptData->nSize*sizeof(char *));
		ptData->aszSeq		= (char **)realloc(ptData->aszSeq, ptData->nSize*sizeof(char *));
		ptData->aszFrag		= (char **)realloc(ptData->aszFrag, ptData->nSize*sizeof(char *));
		ptData->anLen		= (int *)realloc(ptData->anLen, ptData->nSize*sizeof(int));
		ptData->anAbund		= (int64_t *)realloc(ptData->anAbund, ptData->nSize*sizeof(int64_t));
		ptData->anNmera		= (int *)realloc(ptData->anNmera, ptData->nSize*sizeof(int));
		ptData->anID		= (int *)realloc(ptData->anID, ptData->nSize*sizeof(int));
		ptData->anParentA	= (int *)realloc(ptData->anParentA, ptData->nSize*sizeof(int));
		ptData->anParentB	= (int *)realloc(ptData->anParentB, ptData->nSize*sizeof(int));
		ptData->anBreak		= (int *)realloc(ptData->anBreak, ptData->nSize*sizeof(int));
		ptData->adWeight	= (double *)realloc(ptData->adWeight, ptData->nSize*sizeof(double));
		ptData->adProb		= (double *)realloc(ptData->adProb, ptData->nSize*sizeof(double));
		
	}
	
	if(frag == FALSE){
		
		ptData->anID[ptData->nSeq - 1] = ptData->nNextID;
		ptData->nNextID++;
		
	}
	else{
		
		ptData->anID[ptData->nSeq - 1] = 0;
		
	}
	
	ptData->aszSeq[ptData->nSeq - 1] = strdup(seq);	//add info for new sequence
	ptData->aszFrag[ptData->nSeq - 1] = strdup(fragment);
	ptData->aszID[ptData->nSeq - 1] = strdup(name);
	ptData->anLen[ptData->nSeq - 1] = strlen(seq);
	ptData->anAbund[ptData->nSeq - 1] = abund;
	ptData->anNmera[ptData->nSeq - 1] = chimera;
	ptData->anParentA[ptData->nSeq - 1] = parentA;
	ptData->anParentB[ptData->nSeq - 1] = parentB;
	ptData->anBreak[ptData->nSeq - 1] = bp;
	ptData->adWeight[ptData->nSeq - 1] = weight;
	ptData->adProb[ptData->nSeq - 1] = prob;
	
	if(strlen(seq) > ptData->nMaxLen){
		
		ptData->nMaxLen = strlen(seq);	//if new sequence is longest, change the maximum sequence length record
		
	}
	if(abund > ptData->nMaxAbund){
		
		ptData->nMaxAbund = abund;
		
	}
}

void emptyData(t_Data *ptData){		//Set all abundances to zero
	int i;
	for(i = 0; i < ptData->nSeq; i++){
		ptData->anAbund[i] = 0;
	}
	ptData->nTotAbund = 0;
	ptData->nMaxAbund = 0;
}

void removeSeq(size_t nID, int64_t abund, t_Data *ptData)	//Remove abund copies of sequence nID from set ptData
{
	
	ptData->anAbund[nID-1] = ptData->anAbund[nID-1] - abund;	//decrease sequence abundance by abund
	ptData->nTotAbund = ptData->nTotAbund - abund;			//decrease total abundance by abund
	
}

void refreshData(t_Data *ptData)	//Removes sequences with zero abundance
{
	int nS, i, k;
	nS = ptData->nSeq;	//number of sequences
	ptData->nMaxLen = 0;
	ptData->nMaxAbund = 0;
	
	for(k = 0; k < nS; k++){
		
		if(ptData->anAbund[k] == 0){
			ptData->nSeq--;		//decrease number of sequences by one
			
			for(i = k; i < ptData->nSeq; i++){	//starting at sequence to be removed - shift info for all later sequences up by one

				ptData->anLen[i] = ptData->anLen[i+1];
				ptData->anNmera[i] = ptData->anNmera[i+1];
				ptData->anAbund[i] = ptData->anAbund[i+1];
				ptData->anID[i] = ptData->anID[i+1];
				ptData->anParentA[i] = ptData->anParentA[i+1];
				ptData->anParentB[i] = ptData->anParentB[i+1];
				ptData->adWeight[i] = ptData->adWeight[i+1];
				ptData->adProb[i] = ptData->adProb[i+1];
				ptData->anBreak[i] = ptData->anBreak[i+1];
				
				free(ptData->aszSeq[i]);
				ptData->aszSeq[i] = strdup(ptData->aszSeq[i+1]);
				free(ptData->aszID[i]);
				ptData->aszID[i] = strdup(ptData->aszID[i+1]);
				free(ptData->aszFrag[i]);
				ptData->aszFrag[i] = strdup(ptData->aszFrag[i+1]);
				
			}
			free(ptData->aszSeq[ptData->nSeq]);
			free(ptData->aszID[ptData->nSeq]);
			free(ptData->aszFrag[ptData->nSeq]);
		}
		else{
			if(ptData->anLen[k] > ptData->nMaxLen){
				
				ptData->nMaxLen = ptData->anLen[k];
				
			}
			if(ptData->anAbund[k] > ptData->nMaxAbund){
				
				ptData->nMaxAbund = ptData->anAbund[k];
			}
		}
	}
}

void printData(t_Data *ptData)		//Outputs data in fasta format
{
	int i;
	int chicount = 0;
	for(i=0; i < ptData->nSeq; i++){
		
		printf(">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
		if(ptData->anNmera[i] > 0){
			chicount++;
		}
	
	}
	printf("# Sequences = %i\n# Chimeras = %i\nMax Length = %i\nTotal Abundance = %lli\n\n", ptData->nSeq, chicount, ptData->nMaxLen, ptData->nTotAbund);
}

void printChim(t_Data *ptData)		//Outputs chimeras only
{
	int i;
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] > 0){
			
			printf(">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
		
		}
	}
	printf("# Sequences = %i\nMax Length = %i\nTotal Abundance = %lli\n\n", ptData->nSeq, ptData->nMaxLen, ptData->nTotAbund);
}

void printGood(t_Data *ptData)		//Outputs good sequences only
{
	int i;
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] == 0){

			printf(">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
			
		}
	}
	printf("# Sequences = %i\nMax Length = %i\nTotal Abundance = %lli\n\n", ptData->nSeq, ptData->nMaxLen, ptData->nTotAbund);
}

void printParents(t_Data *ptData)		//Outputs chimera and its parents
{
	int i;
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] > 0){
			printf(">chimera_%s_%i\n\%s\n", ptData->aszID[i], ptData->anID[i], ptData->aszSeq[i]);
			printf(">parentA_%s_%i\n\%s\n", ptData->aszID[ptData->anParentA[i]], ptData->anID[ptData->anParentA[i]], ptData->aszSeq[ptData->anParentA[i]]);
			printf(">parentB_%s_%i\n\%s\n", ptData->aszID[ptData->anParentB[i]], ptData->anID[ptData->anParentB[i]], ptData->aszSeq[ptData->anParentB[i]]);
		}
		
	}
	
}

void printBP(t_Data *ptData)		//Outputs a list of chimera break points
{
	int i;
	if(ptData->anBreak[0] > 0){
		printf("%i\n",ptData->anBreak[0]);
	}
	for(i = 1; i < ptData->nSeq; i++){
		if(ptData->anBreak[i] > 0){
			printf("%i\n",ptData->anBreak[i]);
		}
	}
}

void printAbund(t_Data *ptData)		//Outputs a list of abundances
{
	int i;
	double dAbund;
	for(i = 0; i < ptData->nSeq; i++){
		
		dAbund = ((double)ptData->anAbund[i])/((double)ptData->nTotAbund);
		printf("%f\n",dAbund);
	
	}
}

void readData(char* szInputFile, t_Data *ptData)	//Reads in sequences from a fasta file
{													//Adapted from Perseus function
	
	FILE *ifp = NULL;
	char szLine[MAX_LINE_LENGTH];
	int  nPos = 0, i = 0, j = 0, k = 0, nM = 0, nSequences = 0, nNext =0;
	char *szBrk;  
	char *szRet;
	char *sTemp;
	char *sTemp2;
	
	/*first count sequences and get length*/
	ptData->nSeq    = 0;
	ptData->nMaxLen = 0;
	ptData->nMaxAbund = 0;
	
	ifp = fopen(szInputFile, "r");
	
	if(ifp){
		
		while(fgets(szLine, MAX_LINE_LENGTH, ifp)){
			
			if(szLine[0] == '>'){
				
				if(nPos > ptData->nMaxLen){
					
					ptData->nMaxLen = nPos;
					
				}
				
				ptData->nSeq++;
				
				nPos = 0;
				
			}
			
			else{
				
				i = 0;
				
				while(strrchr(szSequence,szLine[i]) != NULL){
					
					i++;
					nPos++;
				}
			}
		}
		
		fclose(ifp);
	}
	else{
		
		fprintf(stderr, "Can't open input file %s\n", szInputFile);
		exit(EXIT_FAILURE);
	}
	
	if(nPos > ptData->nMaxLen){
		
		ptData->nMaxLen = nPos;
		
	}
	
	ptData->nNextID = ptData->nSeq + 1;
	
	nM = ptData->nMaxLen;
	ptData->nSize = 10000;
	ptData->nTotAbund = 0;
	
	ptData->aszID		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->aszSeq		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->aszFrag		= (char **) malloc(ptData->nSize*sizeof(char *));
	ptData->anLen       = (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anAbund		= (int64_t *) malloc(ptData->nSize*sizeof(int64_t));
	ptData->anNmera		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anID		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anParentA	= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anParentB	= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->anBreak		= (int *)  malloc(ptData->nSize*sizeof(int));
	ptData->adWeight	= (double *) malloc(ptData->nSize*sizeof(double));
	ptData->adProb		= (double *) malloc(ptData->nSize*sizeof(double));
	
	ifp = fopen(szInputFile, "r");
	
	if(ifp){
		
		while(szRet = fgets(szLine, MAX_LINE_LENGTH, ifp)){
			
			if(szLine[0] == '>'){
				
				if(nSequences > 0){
					
					ptData->anLen[nSequences - 1] = nPos;
					
				}
				
				j = 0;
				
				while(szLine[j] != '\0'){
					
					if(szLine[j] == '_'){
						
						k = j+1;
						
					}
					
					j++;
					
				}
				
				sTemp = strdup(szLine + k);
				ptData->anAbund[nSequences] = atoi(sTemp);
				
				if(ptData->anAbund[nSequences] > ptData->nMaxAbund){
					
					ptData->nMaxAbund = ptData->anAbund[nSequences];
					
				}
				
				ptData->nTotAbund = ptData->nTotAbund + atoi(sTemp);
				ptData->anNmera[nSequences] = 0;
				ptData->anParentA[nSequences] = 0;
				ptData->anParentB[nSequences] = 0;
				ptData->anBreak[nSequences] = 0;
				ptData->anID[nSequences] = nSequences + 1;
				ptData->aszFrag[nSequences] = "";
				ptData->adWeight[nSequences] = 0;
				ptData->adProb[nSequences] = 0;
				
				szBrk = strpbrk(szLine, " \n");
				(*szBrk) = '\0';
				ptData->aszID[nSequences] = strdup(szLine + 1);
				nPos = 0;
				nSequences++;
				
				free(sTemp);
				
			}
			
			i = 0;
			
			sTemp2 = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
			
			while(szLine[i] != '\0' && strrchr(szSequence,szLine[i]) != NULL){
				
				sTemp2[i] = toupper(szLine[i]);
				nPos++; i++;
				
			}
			
			if(szLine[0] != '>'){
				
				sTemp2[i] = '\0';
				ptData->aszSeq[nSequences-1] = strdup(sTemp2);
				
			}
			
			free(sTemp2);
		}
		
		ptData->anLen[nSequences - 1] = nPos;
		fclose(ifp);
		
	}
	
	else{
		
		fprintf(stderr, "Can't open input file %s\n", szInputFile);
		exit(EXIT_FAILURE);
		
	}
}


void SimsFWD(char *SeqA, char *SeqB, int *nPos, int *nMax, int nBase)//Number of similarities between SeqA and SeqB if we start at the left
{
	int nA = strlen(SeqA);			//length of seq A
	int nB = strlen(SeqB);			//length of seq B
	int i, j;
	int nScore = 0;
	*nMax = 0;
	*nPos = 0;
	
	for(i = 0; i <= nB-nBase; i++){
		for(j = 0; j < nBase; j++){
			if(SeqA[nA - nBase + j] == SeqB[i+j]){
				nScore++;
			}
			
		}
		if(nScore > *nMax){
			*nMax = nScore;
			*nPos = i + nBase;
		}
		nScore = 0;
	}
}

void SimsREV(char *SeqA, char *SeqB, int *nPos, int *nMax, int nBase)//Number of similarities between SeqA and SeqB if we start at the right
{
	int nA = strlen(SeqA);			//length of seq A
	int nB = strlen(SeqB);			//length of seq B
	int i, j;
	int nScore = 0;
	*nMax = 0;
	*nPos = 0;
	
	for(i = 0; i <= nB-nBase; i++){
		for(j = 0; j < nBase; j++){
			if(SeqA[j] == SeqB[i+j]){
				nScore++;
			}
		}
		if(nScore > *nMax){
			*nMax = nScore;
			*nPos = i;
		}
		nScore = 0;
	}
}

double SimsPrimerFWD(char *SeqA, char *SeqB)//Number of similarities between Primer A and Seq B if we start at the left 
{
	int nA = strlen(SeqA);			//Length of seq A
	int i;
	double dScore = 0;				//Records similarities
	
	for(i = 0; i < nA; i++){
		
		if(SeqB[i] == 'A'){
			if(SeqA[i] == 'A'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'R' || SeqA[i] == 'W' || SeqA[i] == 'M'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'D' || SeqA[i] == 'H' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		
		
		if(SeqB[i] == 'C'){
			if(SeqA[i] == 'C'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'Y' || SeqA[i] == 'S' || SeqA[i] == 'M'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'H' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		
		if(SeqB[i] == 'G'){
			if(SeqA[i] == 'G'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'R' || SeqA[i] == 'S' || SeqA[i] == 'K'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'D' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		
		if(SeqB[i] == 'T'){
			if(SeqA[i] == 'T'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'Y' || SeqA[i] == 'W' || SeqA[i] == 'K'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'D' || SeqA[i] == 'H'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
	}
	return dScore;
}

double SimsPrimerREV(char *SeqA, char *SeqB)//Number of similarities between Primer A and Seq B if we start at the left
{
	int nA = strlen(SeqA);	//Length of seq A
	int nB = strlen(SeqB);	//Length of primer
	int i;
	double dScore = 0;		//Records similarities
	
	for(i = 0; i < nA; i++){
		if(SeqB[nB - nA + i] == 'A'){
			if(SeqA[i] == 'A'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'R' || SeqA[i] == 'W' || SeqA[i] == 'M'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'D' || SeqA[i] == 'H' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		if(SeqB[nB - nA + i] == 'C'){
			if(SeqA[i] == 'C'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'Y' || SeqA[i] == 'S' || SeqA[i] == 'M'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'H' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		if(SeqB[nB - nA + i] == 'G'){
			if(SeqA[i] == 'G'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'R' || SeqA[i] == 'S' || SeqA[i] == 'K'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'D' || SeqA[i] == 'V'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
		
		if(SeqB[nB - nA + i] == 'T'){
			if(SeqA[i] == 'T'){
				dScore = dScore+1;
			}
			if(SeqA[i] == 'Y' || SeqA[i] == 'W' || SeqA[i] == 'K'){
				dScore = dScore + 0.5;
			}
			if(SeqA[i] == 'B' || SeqA[i] == 'D' || SeqA[i] == 'H'){
				dScore = dScore + (1/3);
			}
			if(SeqA[i] == 'N'){
				dScore = dScore + 0.25;
			}
		}
	}
	return dScore;
}

void Sample(gsl_rng *ptRNG, t_Data *ptIn, t_Data *ptOut, int64_t nChoose){		//Generate sample with total abundance of nChoose
	
	int i;
	unsigned int nChoose_B;
	unsigned int *anSamp;
	double *adProb;
	
	nChoose_B = (unsigned int)nChoose;
	anSamp = (unsigned int *) malloc(ptIn->nSeq*sizeof(unsigned int));
	adProb = (double *) malloc(ptIn->nSeq*sizeof(double));
	
	for(i = 0; i < ptIn->nSeq; i++){
		adProb[i] = (double)ptIn->anAbund[i];
	}
	
	gsl_ran_multinomial(ptRNG, ptIn->nSeq, nChoose_B, adProb, anSamp);
	
	for(i = 0; i < ptIn->nSeq; i++){
		if(anSamp[i] > 0){
			
			addSeq(ptIn->aszSeq[i], ptIn->aszFrag[i], ptIn->aszID[i], (int64_t)anSamp[i], ptIn->anNmera[i], ptIn->anParentA[i], ptIn->anParentB[i], FALSE, ptIn->anBreak[i], ptIn->adWeight[i], ptIn->adProb[i], ptOut);
		
		}
	}
	
	free(anSamp);
	
}

void writeData(t_Data *ptData, char* fileName)		//Outputs data in fasta format
{
	int i;
	int chicount = 0;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i=0; i < ptData->nSeq; i++){
		fprintf(out_file,">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
		if(ptData->anNmera[i] > 0){
			chicount++;
		}
	}
	fclose(out_file);
}

void writeChim(t_Data *ptData, char* fileName)		//Outputs chimeras only
{
	int i;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] > 0){
			
			fprintf(out_file,">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
			
		}
	}
	fclose(out_file);
}

void writeGood(t_Data *ptData, char* fileName)		//Outputs good sequences only
{
	int i;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] == 0){
			
			fprintf(out_file,">%s_%i_%lli\n\%s\n", ptData->aszID[i],ptData->anID[i],ptData->anAbund[i],ptData->aszSeq[i]);
			
		}
	}
	fclose(out_file);
}

void writeParents(t_Data *ptData, char* fileName)		//Outputs chimera and its parents
{
	int i, j;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] > 0){
			fprintf(out_file,">chimera_%s_%i\n\%s\n", ptData->aszID[i], ptData->anID[i], ptData->aszSeq[i]);
			for(j=0; j < ptData->nSeq; j++){
				if(ptData->anParentA[i] == ptData->anID[j]){
					fprintf(out_file,">parentA_%s_%i\n\%s\n", ptData->aszID[j], ptData->anID[j], ptData->aszSeq[j]);
				}
			}
			for(j=0; j < ptData->nSeq; j++){
				if(ptData->anParentB[i] == ptData->anID[j]){
					fprintf(out_file,">parentB_%s_%i\n\%s\n", ptData->aszID[j], ptData->anID[j], ptData->aszSeq[j]);
				}
			}
		
		}
	}
	fclose(out_file);
}

void writeBP(t_Data *ptData, char* fileName)		//Outputs a list of chimera break points
{
	int i;
	
	FILE *out_file = fopen(fileName, "w");
	
	if(ptData->anBreak[0] > 0){
		fprintf(out_file,"%i",ptData->anBreak[0]);
	}
	for(i = 1; i < ptData->nSeq; i++){
		if(ptData->anBreak[i] > 0){
			fprintf(out_file,"\n%i",ptData->anBreak[i]);
		}
	}
	fclose(out_file);
}

void writeAbund(t_Data *ptData, char* fileName)		//Outputs a list of abundances
{
	int i;
	double dAbund;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i = 0; i < ptData->nSeq; i++){
		
		dAbund = ((double)ptData->anAbund[i])/((double)ptData->nTotAbund);
		fprintf(out_file,"%f\n",dAbund);
		
	}
	
	fclose(out_file);
}

void writeInputInfo(int nSeed, char* szInput, char* szVersion, int nRounds, double dLambda, int nSamp, char* PrimerFWD, char* PrimerREV, char* fileName)		//Input information
{
	FILE *out_file = fopen(fileName, "w");
	fprintf(out_file,"Version: %s\n", szVersion);
	fprintf(out_file,"RNG seed: %i\n", nSeed);
	fprintf(out_file,"Input file: %s\n", szInput);
	fprintf(out_file,"Simulated rounds: %i\n", nRounds);
	fprintf(out_file,"Lambda: %e\n", dLambda);
	fprintf(out_file,"Sampled reads: %i\n", nSamp);
	fprintf(out_file,"Forward primer: %s\n", PrimerFWD);
	fprintf(out_file,"Reverse primer: %s\n", Rev_Comp(PrimerREV));
	
	fclose(out_file);
}

void writeSummary(t_Data *ptData, char* fileName)		//Output statistics
{
	int i;
	int chicount = 0;
	int64_t chiabund = 0;
	
	FILE *out_file = fopen(fileName, "w");
	
	for(i=0; i < ptData->nSeq; i++){
		
		if(ptData->anNmera[i] > 0){
			chicount++;
			chiabund = chiabund + ptData->anAbund[i];
		}
	}
	fprintf(out_file,"# Sequences = %i\n# Chimeras = %i\nMax length = %i\nTotal abundance = %lli\nChimera abundance = %lli", ptData->nSeq, chicount, ptData->nMaxLen, ptData->nTotAbund, chiabund);
	
	fclose(out_file);
}

char* Rev_Comp(char* SeqA)		//Returns reverse complement of SeqA
{
	int i;
	int nA = strlen(SeqA);			//length of seq A
	char *SeqB = strdup(SeqA);
	
	for(i = 0; i < nA; i++){
		
		if(SeqA[i] == 'A'){
			SeqB[nA-i-1] = 'T';
		}
		if(SeqA[i] == 'C'){
			SeqB[nA-i-1] = 'G';
		}
		if(SeqA[i] == 'G'){
			SeqB[nA-i-1] = 'C';
		}
		if(SeqA[i] == 'T'){
			SeqB[nA-i-1] = 'A';
		}
		if(SeqA[i] == 'Y'){
			SeqB[nA-i-1] = 'R';
		}
		if(SeqA[i] == 'R'){
			SeqB[nA-i-1] = 'Y';
		}
		if(SeqA[i] == 'S'){
			SeqB[nA-i-1] = 'S';
		}
		if(SeqA[i] == 'W'){
			SeqB[nA-i-1] = 'W';
		}
		if(SeqA[i] == 'K'){
			SeqB[nA-i-1] = 'M';
		}
		if(SeqA[i] == 'M'){
			SeqB[nA-i-1] = 'K';
		}
		if(SeqA[i] == 'B'){
			SeqB[nA-i-1] = 'V';
		}
		if(SeqA[i] == 'D'){
			SeqB[nA-i-1] = 'H';
		}
		if(SeqA[i] == 'H'){
			SeqB[nA-i-1] = 'D';
		}
		if(SeqA[i] == 'V'){
			SeqB[nA-i-1] = 'B';
		}
		if(SeqA[i] == 'N'){
			SeqB[nA-i-1] = 'N';
		}
	}
	return SeqB;
}