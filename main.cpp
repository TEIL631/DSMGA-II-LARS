/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/


#include <math.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"

using namespace std;

double alpha = 0.5;

int
main (int argc, char *argv[], char *envp[]) {
    if (argc != 9 && argc != 11 && argc != 12) {
        printf ("DSMGA2 ell(lvl) nInitial function maxGen maxFe repeat display rand_seed [m] [n] [coeff]\n");
        printf ("function: \n");
        printf ("     ONEMAX     :  0\n");
        printf ("     MK         :  1\n");
        printf ("     FTRAP      :  2\n");
        printf ("     CYC        :  3\n");
        printf ("     NK         :  4\n");
        printf ("     SPIN       :  5\n");
        printf ("     SAT        :  6\n");
	printf ("     MAXCUT     :  7\n");
	printf ("     HTRAP      :  8\n");
	printf ("     HXOR       :  9\n");
	printf ("     HIFF       : 10\n");

        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    int nInitial = atoi (argv[2]); // initial population size
    int fffff = atoi (argv[3]); // function
    int maxGen = atoi (argv[4]); // max generation
    int maxFe = atoi (argv[5]); // max fe
    int repeat = atoi (argv[6]); // how many time to repeat
    int display = atoi (argv[7]); // display each generation or not
    int rand_seed = atoi (argv[8]);  // rand seed
    int m = 1;
    int n = 1;
    int lvl = 0;
    int inst_num = 1;

    if(argc == 11 || argc == 12){
	if(fffff != 0 && fffff != 1 && fffff != 2 && fffff != 3){
	    printf("BBs for leading series only support onemax, mk, cyc, folded currently");
	    return -1;
	}
	lvl = ell;
        m = atoi(argv[9]);
	n = atoi(argv[10]);

	if(argc == 12){
	    double coeff;
	    coeff = atof(argv[11]);
	    ell = myLD.configure(fffff,m,n,lvl,coeff);
	}
	else
	    ell = myLD.configure(fffff,m,n,lvl);

	// printf("effective ell = %d\n", ell);
	//printf("configuration done\n");

	fffff = Chromosome::Function::LEADING;
        //LEADING = true;
    }

    char* inst_num_env = getenv("DSMGA2_INSTANCE_NUMBER");
    if(inst_num_env)
	inst_num = atoi(inst_num_env);

    if (fffff == 4) {

        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 1, inst_num);
        //sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 5 , 1);

        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 5) {
        char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, inst_num);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf", ell, ell, inst_num);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }

    if (fffff == 7){
	char filename[200];
	char opt_filename[200];
	sprintf(filename,"./MAXCUT/w05_%d/w05_%d.%d", ell, ell, inst_num);
	sprintf(opt_filename, "./MAXCUT/g_w05_%d/g_w05_%d.%d", ell, ell, inst_num);
	if (SHOW_BISECTION) printf("Loading: %s\n", filename);
	if (SHOW_BISECTION) printf("Loading groundtruth: %s\n", opt_filename);
	loadMAXCUT(filename, opt_filename, &myMAXCUT);
    }


    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen, stFE, stLSFE, stRE;
    int usedGen;

    int failNum = 0;

    for (i = 0; i < repeat; i++) {

        DSMGA2 ga (ell, nInitial, maxGen, maxFe, fffff);

        if (display == 1)
            usedGen = ga.doIt (true);
        else
            usedGen = ga.doIt (false);


        if (!ga.foundOptima()) {
            failNum++;
            printf ("-");
        } else {
            stFE.record (Chromosome::hitnfe);
            stLSFE.record (Chromosome::lsnfe);
            stGen.record (usedGen);
            stRE.record(ga.revive_count);
            printf ("+");
        }

        fflush (NULL);

    }
    cout<<endl; 
    printf ("%f  %f  %f  %d  %f\n", stGen.getMean (), stFE.getMean(), stLSFE.getMean(), failNum, stRE.getMean());

    if (fffff == 4) freeNKWAProblem(&nkwa);

    return EXIT_SUCCESS;
}
