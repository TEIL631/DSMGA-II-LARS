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
// bool BACKTRACK = true;
// int GHC_type = 2; // 0:GHC, 1:partial GHC, 2: ordered partial GHC 
double Partial_GHC_p = 0.3;
// int Trimming = 0; // 0: disabled, 1: (i-min)/max, 2: sqrt(max-min/pop)*(i-min)/(max-min), 
// bool FIX_CRITERIA = 0; // 0: cover rate, 1: success rate
// int STALL_CRITERIA = 2; // 0: same nfe, 1: pop unchanged, 2: same nfe + non-increasing mi
// bool STALL_DYN_THRES = 0; // 0: constant stall threshold, 1: decreasing (halved after every stall)
int STALL_THRES = 50; // initial value for stall threshold
// int STALL_THRES_MIN = 10; // minimum value for decreasing stall threshold
// bool RM_REINIT = false; // only select reinit bit as startNode in RM
// bool SC = false; // similarity check
bool HARD_STOP = false; // forced termination if recent HARD_STOP_COEFF*maxnfe nfe do not yield better fitness
double HARD_STOP_COEFF = 0.4; // % of maxnfe

bool to_txt = false;
FILE* pop_out;
FILE* mi_out;
FILE* reinit_out;
FILE* time_out;
clock_t start_time;


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

