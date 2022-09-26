/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#define MAX_GEN 2000

int step = 30;

double alpha = 0.9;

using namespace std;

struct Record {
    int n;
    double nfe;
    double gen;
    double revive;
    double reverse;
};



int main (int argc, char *argv[]) {

    if (argc != 4 && argc!=5 && argc !=6 && argc != 7) {
        printf ("sweep ell(lvl) numConvergence function(0~3) [m] [n] [coeff]\n");
        printf ("sweep ell numConvergence 4 [step #] [nk problem #]\n");
        printf ("sweep ell numConvergence 5 [spin problem #]\n");
        printf ("sweep ell numConvergence 6 [sat problem #]\n");
    printf ("sweep ell numConvergence 7 [maxcut problem #]\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");
    printf ("     MAXCUT:  7\n");
    printf ("     HTRAP :  8\n");
    printf ("     HXOR  :  9\n");
    printf ("     HIFF  : 10\n");
        return -1;
    }

    int ell = atoi (argv[1]);
    int numConvergence = atoi (argv[2]); // problem size
    int fffff = atoi(argv[3]);
    int m = 1;
    int n = 1;
    int lvl = 0;

    int problemNum = 0;
    int neighborNum = 0;
    int stepNum = 0;

    if (fffff == 0 || fffff == 1 || fffff == 2 || fffff == 3){
        if (argc >= 6){
            m = atoi(argv[4]);
            n = atoi(argv[5]);
            lvl = ell;

            if (argc == 7){
                double coeff;
                coeff = atof(argv[6]);
                ell = myLD.configure(fffff,m,n,lvl,coeff);
            }
            else{
                ell = myLD.configure(fffff,m,n,lvl);
            }
            fffff = Chromosome::Function::LEADING;
        }
    }

    if (fffff == 4) {
        neighborNum = 4;
        stepNum = atoi (argv[4]);
        problemNum = atoi (argv[5]);
    }

    if (fffff == 5 || fffff == 6 || fffff == 7) {
        problemNum = atoi (argv[4]);
    }


    int nInitial = 10;


    // for debug
    // myRand.seed(123);


    Statistics st;

    Statistics stGen, stLS, stNFE, stRevive, stReverse;


    if (fffff == 5) {
    char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);
    }

    if (fffff == 4) {
        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, neighborNum, stepNum, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf",ell,ell,problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }

    if (fffff == 7){
        char filename[200];
    char opt_filename[200];
    sprintf(filename, "./MAXCUT/w05_%d/w05_%d.%d", ell, ell, problemNum);
    sprintf(opt_filename, "./MAXCUT/g_w05_%d/g_w05_%d.%d", ell, ell, problemNum);
    if(SHOW_BISECTION) printf("Loading: %s\n", filename);
    if(SHOW_BISECTION) printf("Loading groundtruth: %s\n", opt_filename);
    loadMAXCUT(filename, opt_filename, &myMAXCUT);
    }


    bool foundOptima;
    Record rec[3];
    rec[0].n = nInitial;
    rec[1].n = nInitial+step;
    rec[2].n = nInitial+step+step;

    int popu;
    Record q1, q3;

    if (SHOW_BISECTION) printf("Bisection phase 1\n");

    for (int i=0; i<3; ++i) {
        popu = rec[i].n;

        if (SHOW_BISECTION) printf("[%d]: ", popu);

        foundOptima = true;

        stGen.reset();
        stNFE.reset();
        stLS.reset();
        stRevive.reset();
        stReverse.reset();

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt(false);

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);
            stRevive.record(ga.revive_count);
            stReverse.record(ga.reverse_count);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[i].gen = stGen.getMean();

        if (!foundOptima)
            rec[i].nfe = INF;
        else
            rec[i].nfe = stNFE.getMean();

        rec[i].revive = stRevive.getMean();
        rec[i].reverse = stReverse.getMean();

        if (SHOW_BISECTION) printf(" : %f %f %f\n", rec[i].nfe, rec[i].revive, rec[i].reverse);

    }

    while (rec[0].nfe < rec[1].nfe  && ((rec[2].n-rec[0].n)*20 > rec[1].n)) {

        rec[2] = rec[1];
        rec[1].n = (rec[0].n + rec[2].n) / 2;
        step /= 2;
        popu = rec[1].n;

        if (SHOW_BISECTION) printf("[%d]: ", popu);

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt(false);

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);
            stRevive.record(ga.revive_count);
            stReverse.record(ga.reverse_count);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[1].gen = stGen.getMean();

        if (!foundOptima)
            rec[1].nfe = INF;
        else
            rec[1].nfe = stNFE.getMean();

        rec[1].revive = stRevive.getMean();
        rec[1].reverse = stReverse.getMean();

        if (SHOW_BISECTION) printf(" : %f %f %f\n", rec[1].nfe, rec[1].revive, rec[1].reverse);

    }


    while ( (rec[1].nfe >= rec[0].nfe) || (rec[1].nfe >= rec[2].nfe)) {

        popu = rec[2].n + step;

        if (SHOW_BISECTION) printf("[%d]: ", popu);

        foundOptima = true;

        stGen.reset();
        stNFE.reset();
        stLS.reset();
        stRevive.reset();
        stReverse.reset();

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt(false);

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);
            stRevive.record(ga.revive_count);
            stReverse.record(ga.reverse_count);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[0] = rec[1];
        rec[1] = rec[2];
        rec[2].n = popu;
        rec[2].gen = stGen.getMean();

        if (!foundOptima)
            rec[2].nfe = INF;
        else
            rec[2].nfe = stNFE.getMean();

        rec[2].revive = stRevive.getMean();
        rec[2].reverse = stReverse.getMean();

        if (SHOW_BISECTION) printf(" : %f %f %f\n", rec[2].nfe, rec[2].revive, rec[2].reverse);

    }


    if (SHOW_BISECTION) printf("Bisection phase 2\n");

    while ( ((rec[2].n-rec[0].n)*20 > rec[1].n) && (rec[2].n>rec[1].n+1) && (rec[1].n>rec[0].n+1)) {

        q1.n = (rec[0].n + rec[1].n) / 2;

        if (SHOW_BISECTION) printf("[%d]: ", q1.n);

        foundOptima = true;

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, q1.n, MAX_GEN, -1, fffff);
            ga.doIt(false);

            if (!ga.foundOptima()) {
                foundOptima = false;
                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
            if (j==0) {
                stGen.reset();
                stLS.reset();
                stNFE.reset();
                stRevive.reset();
                stReverse.reset();
            }
            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);
            stRevive.record(ga.revive_count);
            stReverse.record(ga.reverse_count);
        }

        q1.gen = stGen.getMean();
        if (foundOptima)
            q1.nfe = stNFE.getMean();
        else
            q1.nfe = INF;
        q1.revive = stRevive.getMean();
        q1.reverse = stReverse.getMean();


        if (SHOW_BISECTION) printf(" : %f %f %f\n", q1.nfe, q1.revive, q1.reverse);


        q3.n = (rec[1].n + rec[2].n) / 2;

        if (SHOW_BISECTION) printf("[%d]: ", q3.n);

        foundOptima = true;

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, q3.n, MAX_GEN, -1, fffff);
            ga.doIt(false);

            if (!ga.foundOptima()) {
                foundOptima = false;
                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
            if (j==0) {
                stGen.reset();
                stLS.reset();
                stNFE.reset();
                stRevive.reset();
                stReverse.reset();
            }
            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);
            stRevive.record(ga.revive_count);
            stReverse.record(ga.reverse_count);
        }

        q3.gen = stGen.getMean();
        if (foundOptima) {
            q3.nfe = stNFE.getMean();
        } else
            q3.nfe = INF;
        q3.revive = stRevive.getMean();
        q3.reverse = stReverse.getMean();

        if (SHOW_BISECTION) printf(" : %f %f %f\n", q3.nfe, q3.revive, q3.reverse);

        if (rec[1].nfe < q1.nfe && rec[1].nfe < q3.nfe) {
            rec[0] = q1;
            rec[2] = q3;
        } else if (q1.nfe < rec[1].nfe && q1.nfe < q3.nfe) {
            rec[2] = rec[1];
            rec[1] = q1;
        } else { // q3nfe smallest
            rec[0] = rec[1];
            rec[1] = q3;
        }
    };



    if (fffff == 4)
        freeNKWAProblem(&nkwa);

    /*printf("population: %d\n", rec[1].n);
    printf("generation: %f\n", rec[1].gen);
    printf("NFE: %f\n", rec[1].nfe);*/

    printf("%d\t", problemNum);
    printf("%d\t", rec[1].n);
    printf("%f\t", rec[1].gen);
    printf("%f\t", rec[1].nfe);
    printf("%f\t", rec[1].revive);
    printf("%f\n", rec[1].reverse);


    return EXIT_SUCCESS;

}

