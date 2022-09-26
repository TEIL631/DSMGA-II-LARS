/***************************************************************************
 *   Copyright (C) 2015 by TEIL                                            *
 ***************************************************************************/

#include <cstdio>
#include <cstring>
#include <time.h>
#include "spin.h"
#include "chromosome.h"
#include "nk-wa.h"
#include "sat.h"

#define TRAP_K 5
#define HTRAP_K 3

extern double alpha;

Chromosome::Chromosome () {
    length = 0;
    lengthLong = 0;
    gene = NULL;
    evaluated = false;
}

Chromosome::Chromosome (int n_length) {
    gene = NULL;
    init (n_length);
}


Chromosome::~Chromosome () {
    if (gene != NULL) delete []gene;
}

void Chromosome::init (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    evaluated = false;
}

void Chromosome::init0 (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];

    for (int i=0; i<lengthLong; ++i)
        gene[i] = 0;

    key = 0;
    evaluated = false;
}

void Chromosome::initR (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    key = 0;
    for (int i=0; i<length; ++i) {

        int val = myRand.flip();
        setValF(i, val);
        if (val == 1)
            key ^= zKey[i];
    }

    evaluated = false;
}

double Chromosome::getFitness () {
    if (evaluated)
        return fitness;
    else {
        fitness = evaluate();
        if (!hit && fitness > getMaxFitness()) {
            hit = true;
            hitnfe = nfe+lsnfe;
        }
	if(fitness - EPSILON > peak_fitness){
	    peak_fitness = fitness - EPSILON;
	    last_improved_nfe = nfe + lsnfe;
	}

        return fitness;
    }
}

bool Chromosome::isEvaluated () const {
    return evaluated;
}

bool Chromosome::hasSeen() const {

    unordered_map<unsigned long, double>::iterator it = cache.find(key);
    if (it != cache.end())
        return true;
    return false;
}

double Chromosome::evaluate () {

    if (CACHE)
        if (hasSeen()) {
            evaluated = true;
            return cache[key];
        }

    ++nfe;
    evaluated = true;
    double accum = 0.0;

    /*if(LEADING){
        char* c = new char[length];
        for(int i=0;i<length;i++){
	    c[i] = getVal(i);
            //printf("%d ",c[i]);
        }
        //printf("\n");
        accum = myLD.evaluate(c);
        //printf("%f\n",f);
	delete[] c;
    }*/

    switch (function) {
        case ONEMAX:
            accum = oneMax();
            break;
        case MKTRAP:
            accum = mkTrap(1, 0.8);
            break;
        case CYCTRAP:
            accum = cycTrap(1, 0.8);
            break;
        case FTRAP:
            accum = fTrap();
            break;
        case SPINGLASS:
            accum = spinGlass();
            break;
        case NK:
            accum = nkFitness();
            break;
        case SAT:
            accum = satFitness();
            break;
	case MAXCUT:
	    accum = maxcutFitness();
	    break;
	case HTRAP:
	    accum = hTrap(1,0.9);
	    break;
	case HXOR:
	    accum = hxor(2);
	    break;
	case HIFF:
	    accum = hiff(2);
	    break;
    	case LEADING:
    	    float tmp_f;
    	    {
    	        char* c = new char[length];
    	        for(int i=0;i<length;i++){
    	            c[i] = getVal(i);
                        //printf("%d ",c[i]);
                    }
                    //printf("\n");
    		tmp_f = myLD.evaluate(c);
    		delete[] c;
    	    }
    	    accum = tmp_f;
    	    break;
    	case LOFT:
    	    accum = loft();
    	    break;
    	case LOFT_OVERLAP:
    	    accum = loft_overlap();
    	    break;
    	default:
    	    accum = mkTrap(1, 0.8);
    	    break;
    }

    if (CACHE)
        cache[key]=accum;


    if((nfe+lsnfe)%100 == 0)
        fprintf(time_out, "%d\t%lf\n", nfe+lsnfe, (clock()-start_time)/(double)(CLOCKS_PER_SEC));

    return accum;

}


double Chromosome::loft () const{
    int len_FT = length/2;
    len_FT -= (len_FT%6);

    double total_score = 0;

    for(int i=0;i<len_FT/6;i++){
        int ones_in_BB = 0;
	for(int j=0;j<6;j++){
	    ones_in_BB += (getVal(i*6+j) == 1);
	}

	int diff = (ones_in_BB >= 3)? ones_in_BB-3 : 3-ones_in_BB;
	switch(diff){
	    case 3: total_score += 1.0 * 6; break;
	    case 2: total_score += 0.0 * 6; break;
	    case 1: total_score += 0.4 * 6; break;
	    case 0: total_score += 0.8 * 6; break;
	}
    }

    for(int i=len_FT;i<length;i++){
        if(getVal(i) == 1)
	    total_score += 1;
	else
	    break;
    }

    return total_score;
}

double Chromosome::loft_overlap () const{
    double total_score = 0;
    total_score += alpha * fTrap();
    
    for(int i=0;i<length;i++){
        if(getVal(i) == 1)
	    total_score += (1-alpha);
	else
	    break;
    }
    return total_score;
}

double
Chromosome::spinGlass () const {

    int *x = new int[length];
    double result;

    for (int i=0; i<length; i++)
        if (getVal(i) == 1)
            x[i] = 1;
        else
            x[i] = -1;

    result = evaluateSPIN(x, &mySpinGlassParams);

    delete []x;

    return result;
}

double Chromosome::nkFitness() const {
    char *x = new char[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = (char) getVal(i);
    }

    double result = evaluateNKProblem(x, &nkwa);
    //double result = evaluateNKWAProblem(x, &nkwa);
    delete []x;
    return result;
}

// OneMax
double Chromosome::oneMax () const {

    double result = 0;

    for (int i = 0; i < length; ++i)
        result += getVal(i);

    return result;
}

bool Chromosome::operator== (const Chromosome& c) const {
    if (length != c.length)
        return false;

    for (int i=0; i<lengthLong; i++)
        if (gene[i] != c.gene[i])
            return false;

    return true;
}

Chromosome& Chromosome::operator= (const Chromosome& c) {

    if (length != c.length) {
        length = c.length;
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;

    memcpy(gene, c.gene, sizeof(long) * lengthLong);

    return *this;
}

double Chromosome::trap (int unitary, double fHigh, double fLow, int trapK) const {
    if (unitary > trapK)
        return 0;

    if (unitary == trapK)
        return fHigh;
    else
        return fLow - unitary * fLow / (trapK-1);
}


double Chromosome::fTrap() const {

    double result = 0.0;

    for (int i=0; i<length/6; ++i) {
        int u=0;
        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    return result;
}

double Chromosome::cycTrap(double fHigh, double fLow) const {
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);
    if (length % (TRAP_K-1) != 0)
        outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        int idx = i * TRAP_K - i;
        for (j = 0; j < TRAP_K; j++) {
            int pos = idx + j;
            if (pos == length)
                pos = 0;
            else if (pos > length)
                outputErrMsg ("CYCLIC BUG");
            //
            u += getVal(pos);
        }
        result += trap (u, fHigh, fLow, TRAP_K);
    }
    return result;
}

double Chromosome::mkTrap (double fHigh, double fLow) const {
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap (u, fHigh, fLow, TRAP_K);
    }

    return result;
}

double Chromosome::hTrap (double high, double low) const {
    int u, _u;

    int t = HTRAP_K;
    while (length>t) {
        t=t*HTRAP_K;
    }
    if (length!=t)
        outputErrMsg ("length != HTRAP_K^k");

    double result = 0;
    
    t = HTRAP_K;
    while (length>=t) {
        int HTRAP_M = length / t;
        int HTRAP_T = int(t/HTRAP_K);

        for (int _i=0; _i<HTRAP_M; _i++) {
            u = 0;
            for (int _j=0; _j<HTRAP_K; _j++) {
                _u = 0;
                for (int _k=0; _k<HTRAP_T; _k++) {
                    _u += getVal(_i*t+_j*HTRAP_T+_k);
                }
                if (_u == HTRAP_T) u += 1;
                else if (_u > 0) {
                    u += (HTRAP_K+1);
                    break;
                }
            }
            if (length>t) {
                result += (int(t/HTRAP_K)*trap (u, high, high, HTRAP_K));
            }
            else {
                result += (int(t/HTRAP_K)*trap (u, high, low, HTRAP_K));
            }
        }
        t=t*HTRAP_K;
    }

    return result;
}

int Chromosome::vxor(int low, int high) const {
    if (low==high) return 1;
    int middle = ((high+1-low)/2) + low;
    if (vxor(low,middle-1) == 0 || vxor(middle,high) == 0) return 0;
    for(int i=0; i<middle-low; i++) if (getVal(low+i) == getVal(middle+i)) return 0;
    return 1;
}

int Chromosome::viff(int low, int high) const {
    if (low==high) return 1;
    int middle = ((high+1-low)/2) + low;
    if (viff(low,middle-1) == 0 || viff(middle,high) == 0) return 0;
    for(int i=0; i<middle-low; i++) if (getVal(low+i) != getVal(middle+i)) return 0;
    return 1;
}

double Chromosome::hxor (int htrap_k) const {
    int t = htrap_k;
    int height = 0;
    while (length>t) {
        t=t*htrap_k;
    }
    if (length!=t) {
        printf("length != HTRAP_K^k\n");
        exit( 0 );
    }

    double result = 0;

    t = 1;
    while (length>=t) {
        int HTRAP_M = length / t;

        for (int _i=0; _i<HTRAP_M; _i++) {
            result += (t*vxor(_i*t,(_i+1)*t-1));
        }
        t=t*htrap_k;
        height++;
    }

    return result;
}

double Chromosome::hiff (int htrap_k) const {
    int t = htrap_k;
    int height = 0;
    while (length>t) {
        t=t*htrap_k;
    }
    if (length!=t) {
        printf("length != HTRAP_K^k\n");
        exit( 0 );
    }

    double result = 0;

    t = 1;
    while (length>=t) {
        int HTRAP_M = length / t;

        for (int _i=0; _i<HTRAP_M; _i++) {
            result += (t*viff(_i*t,(_i+1)*t-1));
        }
        t=t*htrap_k;
        height++;
    }

    return result;
}

int Chromosome::getLength () const {
    return length;
}

double Chromosome::getMaxFitness () const {

    double maxF;
    int t, HTRAP_M;

    switch (function) {
        case ONEMAX:
            maxF = length;
            break;
        case MKTRAP:
            maxF = length/TRAP_K;
            break;
        case FTRAP:
            maxF = length/6;
            break;
        case CYCTRAP:
            maxF =  length/(TRAP_K - 1);
            break;
        case SPINGLASS:
            maxF = mySpinGlassParams.opt;
            break;
        case NK:
            maxF = nkwa.maxF;
            break;
        case SAT:
            maxF = 0;
            break;
	case MAXCUT:
	    maxF = myMAXCUT.opt;
	    break;
	case HTRAP:
	    t = HTRAP_K;
	    maxF = 0;
	    while(length >= t){
	        HTRAP_M = int(length/t);
		maxF += HTRAP_M * int(t/HTRAP_K);
		t = t * HTRAP_K;
	    }
	    break;
	case HXOR:
	    t = 2;
	    HTRAP_M = 1;
	    maxF = 0;
	    while(length >= t){
	        t *= 2;
		HTRAP_M++;
	    }
	    maxF = HTRAP_M * length;
	    break;
	case HIFF:
	    t = 2;
	    HTRAP_M = 1;
	    maxF = 0;
	    while(length >= t){
	        t *= 2;
		HTRAP_M++;
	    }
	    maxF = HTRAP_M * length;
	    break;
	case LEADING:
	    maxF = myLD.get_max();
	    break;
	case LOFT:
	    maxF = length;
	    break;
	case LOFT_OVERLAP:
	    maxF = alpha*length/6 + (1-alpha)*length;
	    break;
	default:
            // Never converge
            maxF = INF;
    }

    return maxF - EPSILON;

}

// contribute to lsnfe
bool Chromosome::tryFlipping(int index) {
    int oldNFE = nfe;

    double oldF = getFitness();
   // printf("id: %d, oldF: %f\n", index,oldF);
    flip(index);

    //2016-10-21
    if (getFitness() - EPSILON <= oldF) {
    //if (getFitness() <= oldF) {
        flip(index);
        evaluated = true;
        fitness = oldF;

        lsnfe += nfe - oldNFE;
        nfe = oldNFE;

        return false;
    } else {

        lsnfe += nfe - oldNFE;
        nfe = oldNFE;

        return true;
    }

}

bool Chromosome::GHC() {

    int* order = new int [length];
    myRand.uniformArray(order, length, 0, length-1);

    bool flag = false;
    for (int i=0; i<length; ++i) {

        if (tryFlipping(order[i])){
    	    flag = true;
    	}
    }

    delete []order;
    return flag;

}

bool Chromosome::PartialGHC(bool* fixed) {

    int* order = new int [length];
    myRand.uniformArray(order, length, 0, length-1);

    bool flag = false;
    for (int i=0; i<length; ++i) {
        if(fixed[order[i]]){
            if(myRand.uniform() < Partial_GHC_p){
                if (tryFlipping(order[i])){
                    flag = true;
                }
            }
        }
        else{
            if (tryFlipping(order[i])){
                flag = true;
            }
        }
    }
    delete []order;
    return flag;
}

bool Chromosome::OrderedPartialGHC(bool* fixed) {

    int* order = new int [length];
    myRand.uniformArray(order, length, 0, length-1);

    bool flag = false;
    for (int i=0; i<length; ++i) {
        if(fixed[order[i]]){
            if(myRand.uniform() < Partial_GHC_p){
                if (tryFlipping(order[i])){
                    flag = true;
                }
            }
        }
    }

    for( int i = 0; i < length; i++){
        if(!fixed[order[i]]){
            if (tryFlipping(order[i])){
                flag = true;
            }
        }
    }

    delete []order;
    return flag;
}

double Chromosome::satFitness() const {
    int *x = new int[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = getVal(i);
    }

    double result = evaluateSAT(x, &mySAT);
    delete []x;
    return result;
}

double Chromosome::maxcutFitness() const {
    char *x = new char[length];

    for(int i=0;i<length;i++)
        x[i] = getVal(i);
    
    double result = evaluateMAXCUT(x, &myMAXCUT);
    delete []x;
    return result;
}

bool Chromosome::similarityCheck(const Chromosome& c) const{
    int dist = 0;
    for(int i=0;i<lengthLong;i++){
        dist += myBD.countOne(gene[i] ^ c.gene[i]);
    }
    return dist <= lengthLong/2;
}
