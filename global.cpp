/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cfloat>
#include "myrand.h"
#include "statistics.h"
#include "doublelinkedlistarray.h"
#include "zkey.h"
#include "chromosome.h"
#include "sat.h"

int maxMemory = 0;

bool GHC = true;
bool SELECTION = true; // FIXME
bool CACHE = false;
bool SHOW_BISECTION = true;
double Partial_GHC_p = 0.3;
int STALL_THRES = 50; // initial value for stall threshold
bool HARD_STOP = false; // forced termination if recent HARD_STOP_COEFF*maxnfe nfe do not yield better fitness
double HARD_STOP_COEFF = 0.4; // % of maxnfe


char outputFilename[100];
Chromosome::Function Chromosome::function;
int Chromosome::nfe;
int Chromosome::lsnfe;
int Chromosome::hitnfe;
bool Chromosome::hit;
unordered_map<unsigned long, double> Chromosome::cache;
double Chromosome::peak_fitness;
int Chromosome::last_improved_nfe;

ZKey zKey;
MyRand myRand;
BitwiseDistance myBD;
SPINinstance mySpinGlassParams;
NKWAProblem nkwa;
SATinstance mySAT;
MAXCUTinstance myMAXCUT;
LEADINGinstance myLD;


void outputErrMsg(const char *errMsg) {
    printf("%s\n", errMsg);
    exit(1);
}

int pow2(int x) {
    return (1 << x);
}

