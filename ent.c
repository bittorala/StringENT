/*
	ENT  --  Entropy calculation and analysis of putative
		 random sequences.

        Designed and implemented by John "Random" Walker in May 1985.

	Multiple analyses of random sequences added in December 1985.

	Bit stream analysis added in September 1997.

	Terse mode output, getopt() command line processing,
	optional stdin input, and HTML documentation added in
	October 1998.
	
	Documentation for the -t (terse output) option added
	in July 2006.
	
	Replaced table look-up for chi square to probability
	conversion with algorithmic computation in January 2008.

	For additional information and the latest version,
	see http://www.fourmilab.ch/random/

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

#include "iso8859.h"
#include "randtest.h"

#define UPDATE  "June 24th, 2020"

#define FALSE 0
#define TRUE  1

#ifdef M_PI
#define PI	 M_PI
#else
#define PI	 3.14159265358979323846
#endif

extern double pochisq(const double ax, const int df);
extern double poz(const double z);

/*  HELP  --  Print information on how to call	*/

static void help(void)
{
        printf("StringENT --  Calculate entropy and basic randomness p-values of file.  Call");
        printf("\n        with ent [options] [input-file]");
        printf("\n");
        printf("\n        Options:   -b   Treat input as a stream of bits");
        printf("\n                   -c   Print occurrence counts");
        printf("\n                   -f   Fold upper to lower case letters");
        printf("\n                   -t   Terse output in CSV format");
        printf("\n                   -u   Print this message\n");
        printf("\nThis is an extension made by Bittor Alana on the battery ENT");
		printf("\nBy John Walker");
        printf("\n   http://www.fourmilab.ch/");
        printf("\n   %s\n", UPDATE);
}

/*  GETOPT  --	Dumb version of getopt for brain-dead Windows.  */

#ifdef _WIN32	
static int optind = 1;

static int getopt(int argc, char *argv[], char *opts)
{
    static char *opp = NULL;
    int o;
    
    while (opp == NULL) {
        if ((optind >= argc) || (*argv[optind] != '-')) {
	   return -1;
	}
	opp = argv[optind] + 1;
	optind++;
	if (*opp == 0) {
	    opp = NULL;
	}	
    }
    o = *opp++;
    if (*opp == 0) { 
	opp = NULL;
    }
    return strchr(opts, o) == NULL ? '?' : o;
}
#endif

/*  Main program  */

int main(int argc, char *argv[])
{
	int i, oc, opt;
	long ccount[256];	      /* Bins to count occurrences of values */
	long totalc = 0;	      /* Total character count */
	long nruns = 0; 		  /* Count of runs */
	long n1 = 0, n2 = 0;	  /* Count of - and + in runs test */
	char *samp;
	double montepi, pipv, chip, meanp, sccp, runp,
	       scc, ent, mean, chisq, lmx2, lmpv;
	int lmnblocks;
	FILE *fp = stdin;
	int counts = FALSE,	      /* Print character counts */
	    fold = FALSE,	      /* Fold upper to lower */
	    binary = FALSE,	      /* Treat input as a bitstream */
	    terse = FALSE;	      /* Terse (CSV format) output */

        while ((opt = getopt(argc, argv, "bcftu?BCFTU")) != -1) {
	    switch (toISOlower(opt)) {
                 case 'b':
		    binary = TRUE;
		    break;

                 case 'c':
		    counts = TRUE;
		    break;

                 case 'f':
		    fold = TRUE;
		    break;

                 case 't':
		    terse = TRUE;
		    break;

                 case '?':
                 case 'u':
		    help();
		    return 0;
	    }
	}
	if (optind < argc) {
	   if (optind != (argc - 1)) {
              printf("Duplicate file name.\n");
	      help();
	      return 2;
	   }
           if ((fp = fopen(argv[optind], "rb")) == NULL) {
              printf("Cannot open file %s\n", argv[optind]);
	      return 2;
	   }
	}
	
#ifdef _WIN32

	    /** Warning!  On systems which distinguish text mode and
		binary I/O (MS-DOS, Macintosh, etc.) the modes in the open
        	statement for "fp" should have forced the input file into
        	binary mode.  But what if we're reading from standard
		input?  Well, then we need to do a system-specific tweak
        	to make sure it's in binary mode.  While we're at it,
        	let's set the mode to binary regardless of however fopen
		set it.

		The following code, conditional on _WIN32, sets binary
		mode using the method prescribed by Microsoft Visual C 7.0
        	("Monkey C"); this may require modification if you're
		using a different compiler or release of Monkey C.	If
        	you're porting this code to a different system which
        	distinguishes text and binary files, you'll need to add
		the equivalent call for that system. */

	    _setmode(_fileno(fp), _O_BINARY);
#endif

        samp = binary ? "bit" : "byte";
	memset(ccount, 0, sizeof ccount);

	/* Initialise for calculations */

	rt_init(binary);

	/* Scan input file and count character occurrences */

	while ((oc = fgetc(fp)) != EOF) {
	   unsigned char ocb;

	   if (fold && isISOalpha(oc) && isISOupper(oc)) {
	      oc = toISOlower(oc);
	   }
	   ocb = (unsigned char) oc;
	   totalc += binary ? 8 : 1;
	   if (binary) {
	    int b;
	    unsigned char ob = ocb;

	    for (b = 0; b < 8; b++) {
		ccount[ob & 1]++;
		ob >>= 1;
	    }
	   } else {
	       ccount[ocb]++;	      /* Update counter for this bin */
	   }
	   rt_add(&ocb, 1);
	}
	fclose(fp);

	/* Complete calculation and return sequence metrics */

	rt_end(&ent, &chisq, &mean, &montepi, &scc, &nruns, &n1, &n2, &lmx2, &lmnblocks);

	/* Calculate probability of observed distribution occurring from
	   the results of the Chi-Square test */

    	chip = pochisq(chisq, (binary ? 1 : 255));

	/* Calculate probability of observed distribution occurring from
	   the results of the Arithmetic Mean test (CLT) */
	   double expmean = binary ? 0.5 : 127.5;
	   double expvariance = binary ? 0.25 : (256.0*256.0 -1) / 12;
	   double aux = sqrt(totalc) * (mean - expmean) / sqrt(expvariance);
	   meanp = 2 * (1 - poz(fabs(aux)));

	/* Calculate probability of serial correlation being this far from zero */
	double t;
	t = scc * sqrt(totalc-2) / sqrt(1-scc*scc); // Get t-statistic from serial correlation
	t = fabs(t);
	
	/*
	Yerukala, Ramu & Boiroju, Naveen & Malkareddy, Krishna. (2013). Approximations to the t-distribution. 
	*/
	aux = t * (1 - 1/(4*totalc)) / sqrt(1 + t*t / (2*totalc)); // For t-distribution CDF approximation
	if (scc >= -99999) {
		sccp = 2 * (1 - poz(aux));
	}
	else {
		sccp = 0;
	}

	/* Runs test p-value */
	double expected_runs = 2*n1*n2 / (n1+n2) + 1;
	double stdev_runs = 2*n1*n2 / ((double)(n1+n2)) *(2*n1*n2-n1-n2) / (n1 + n2) / (n1+n2-1);
	double z = (nruns - expected_runs) / sqrt(stdev_runs);
	z = fabs(z);
	runp = 2*(1-poz(z));

	/* Local means test */
	lmpv = pochisq(lmx2, lmnblocks);

	/* P-value for pi estimation 
		Note that there are bytes / MONTEN coordinates, each of which has probability pi/4
		of being inside a circle. Hence we have a binomial of p=pi/4 and n=bytes/MONTEN
	*/
	long inside, ncoords;
	ncoords = totalc / (binary ? 8 : 1) / MONTEN; /* Binary mode also uses 6 bytes for Monte Carlo */
	inside = montepi / 4 * ncoords;
	z = (inside - PI / 4 * ncoords) / sqrt(ncoords * PI/4.0 * (1-PI/4.0));
	z = fabs(z);
	pipv = 2*(1-poz(z));


	if (terse) {
           printf("0,File-%ss,Entropy,Chi-square,CS-p-value,Mean,M-p-value,Monte-Carlo-Pi,Pi-p-value,Serial-Correlation,SC-p-value,Runs,R-p-value,LM-p-value\n",
              binary ? "bit" : "byte");
           printf("1,%ld,%1.11f,%f,%f,%f,%f,%f,%f,%f,%f,%ld,%f,%f\n",
	      totalc, ent, chisq, chip, mean, meanp, montepi, pipv, scc, sccp, nruns, runp, lmpv);
	}

	/* Print bin counts if requested */

	if (counts) {
	   if (terse) {
              printf("2,Value,Occurrences,Fraction\n");
	   } else {
              printf("Value Char Occurrences Fraction\n");
	   }
	   for (i = 0; i < (binary ? 2 : 256); i++) {
	      if (terse) {
                 printf("3,%d,%ld,%f\n", i,
		    ccount[i], ((double) ccount[i] / totalc));
	      } else {
		 if (ccount[i] > 0) {
                    printf("%3d   %c   %10ld   %f\n", i,
		       /* The following expression shows ISO 8859-1
			  Latin1 characters and blanks out other codes.
			  The test for ISO space replaces the ISO
			  non-blanking space (0xA0) with a regular
                          ASCII space, guaranteeing it's rendered
                          properly even when the font doesn't contain
			  that character, which is the case with many
			  X fonts. */
                       (!isISOprint(i) || isISOspace(i)) ? ' ' : i,
		       ccount[i], ((double) ccount[i] / totalc));
		 }
	      }
	   }
	   if (!terse) {
              printf("\nTotal:    %10ld   %f\n\n", totalc, 1.0);
	   }
	}

	/* Print calculated results */

	if (!terse) {
		printf("StringENT | Results report\n");
		printf("--------------------------------------------------------------------\n");
		printf("Entropy = %1.11f bits per %s.\n", ent, samp);
		printf("\nOptimum compression would reduce the size\n");
		printf("of this %ld %s file by %d percent.\n", totalc, samp,
		(short) ((100 * ((binary ? 1 : 8) - ent) /
				(binary ? 1.0 : 8.0))));
	   printf("--------------------------------------------------------------------\n");
	   printf("Chi square statistic for %ld samples is %1.2f.\n",totalc,chisq);
	   printf("p-value\t\t\t%f  ",chip);
	   if (chip < 0.05) printf("*");
	   if (chip < 0.01) printf("*");
	   if (chip < 0.001) printf("*");
	   printf("\n--------------------------------------------------------------------\n");
	   printf("Arithmetic mean value of data %ss is %1.4f (%.1f = random).\n",samp,mean, binary ? 0.5:127.5);
	   printf("p-value\t\t\t%f  ",meanp);
	   if (meanp < 0.05) printf("*");
	   if (meanp < 0.01) printf("*");
	   if (meanp < 0.001) printf("*");
	   printf("\n--------------------------------------------------------------------\n");
		printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
	      montepi, 100.0 * (fabs(PI - montepi) / PI));
		printf("p-value\t\t\t%f  ",pipv);
		if (pipv < 0.05) printf("*");
	    if (pipv < 0.01) printf("*");
		if (pipv < 0.001) printf("*");
		printf("\n--------------------------------------------------------------------\n");
           printf("Serial correlation coefficient is ");
	   if (scc >= -99999) {
              printf("%1.6f (totally uncorrelated = 0.0).\n", scc);
			  printf("p-value\t\t\t%f  ",sccp);
			  if (sccp < 0.05) printf("*");
		   	  if (sccp < 0.01) printf("*");
			  if (sccp < 0.001) printf("*");
	   } else {
              printf("undefined (all values equal!).\n");
	   }
	   printf("\n--------------------------------------------------------------------\n");
	   printf("The number of runs is %ld.\n", nruns);
	   if (stdev_runs != 0) /* stdev_runs could be zero leading to nan */
	   {
		   printf("p-value\t\t\t%f  ",runp);
		   if (runp < 0.05) printf("*");
		   if (runp < 0.01) printf("*");
		   if (runp < 0.001) printf("*");
	   }
		printf("\n--------------------------------------------------------------------\n");
		printf("The local means test's X^2 statistic is %1.3f for %d blocks.\n", lmx2, lmnblocks);
		printf("p-value\t\t\t%f  ",lmpv);
		if (lmpv < 0.05) printf("*");
		if (lmpv < 0.01) printf("*");
		if (lmpv < 0.001) printf("*");
		printf("\n--------------------------------------------------------------------\n");
		printf("*:\t\tSignificant for α=0.05\n");
		printf("**:\t\tSignificant for α=0.01\n");
		printf("***:\t\tSignificant for α=0.001\n");
	}

	return 0;
}
