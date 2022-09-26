#ifndef LEADING_H
#define LEADING_H

#include <vector>
#include <assert.h>

class BB{
public:
	int type;
	std::vector<int> indices;
	std::vector<BB*> dependent_on;
	bool is_solved;
	bool enable_zeromax;
	double zeromax_coeff;

	BB() = delete;
	BB(int type, std::vector<int> indices, double zeromax_coeff):type(type),indices(indices),is_solved(false),enable_zeromax(true),zeromax_coeff(zeromax_coeff){
	    assert((type==0)?(indices.size()==1):(type==2)?(indices.size()==6):(indices.size()==5));
	}
	BB(int type, std::vector<int> indices):type(type),indices(indices),is_solved(false),enable_zeromax(false){
	    assert((type==0)?(indices.size()==1):(type==2)?(indices.size()==6):(indices.size()==5));
	}

	double evaluate(char* c);
	double evaluate_zeromax(char* c);
};


class LEADINGinstance{
private:
	int type;
	int m;
	int n;
	int lvl;
	std::vector<BB*> bbs;
	int power(int base, int exp){
	    int val=1;
	    for(int i=0;i<exp;i++)
		val*=base;
	    return val;
	}
	double maxf;
	bool enable_zeromax;
	double zeromax_coeff;

public:
	LEADINGinstance():enable_zeromax(false){};
	~LEADINGinstance(){
	    for(auto bb : bbs)
		delete(bb);
	};
	int configure(int type, int m, int n, int lvl, double zeromax_coeff);
	int configure(int type, int m, int n, int lvl);
	double evaluate(char* c);
	double get_max(){return maxf;};
};

#endif
