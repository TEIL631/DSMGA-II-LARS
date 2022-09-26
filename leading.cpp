#include "leading.h"
#include "chromosome.h"

double BB::evaluate_zeromax(char* c){
    double f = 0;
    for(int i=0;i<(int)indices.size();i++){
        if(!c[indices[i]])
	    f += zeromax_coeff;
    }
    return f;
}

double BB::evaluate(char* c){
    //printf("this : %x\n", this);
    for(auto BB_ptr : dependent_on){
	//printf("checked : %d\n",BB_ptr->is_solved);
        if(!(BB_ptr->is_solved)){
	    is_solved = false;
	    return (enable_zeromax)? evaluate_zeromax(c) : 0; // or weak zero max here
	}
    }
    //printf("ckpt\n");

    int onebb_ell = (type==0)?1:(type==2)?6:5;
    Chromosome onebb(onebb_ell);
    for(int i=0;i<onebb_ell;i++){
        onebb.setVal(i,c[indices[i]]);
    }

    Chromosome::function = (Chromosome::Function)type;
    if((Chromosome::Function)type == Chromosome::Function::CYCTRAP) // each bb in cyc is actually just mk
	Chromosome::function = Chromosome::Function::MKTRAP;

    double fitness = onebb.getFitness();
    Chromosome::hit = false;
    Chromosome::nfe--;

    //printf("solved : %d\n", fitness > onebb.getMaxFitness() - EPSILON);
    if(fitness > onebb.getMaxFitness() - EPSILON)
	is_solved = true;
    else
	is_solved = false;

    Chromosome::function = Chromosome::Function::LEADING;

    return fitness;
}

int LEADINGinstance::configure(int type, int m, int n, int lvl, double coeff){
    zeromax_coeff = coeff;
    enable_zeromax = true;
    return configure(type,m,n,lvl);
}

int LEADINGinstance::configure(int type, int m, int n, int lvl){
    this->type = type;
    this->m = m;
    this->n = n;
    this->lvl = lvl;

    assert(m == 1 || n == 1);

    int prev_lvl_BB_start_index = 0;
    int prev_lvl_BBs = 0;
    int current_lvl_BBs = (m<n)? 1:pow(m,lvl-1);
    int ell = 0;
    int current_lvl_first_id;

    maxf = 0;
    for(int current_lvl=0;current_lvl<lvl;current_lvl++){
	current_lvl_first_id = ell;
        for(int i=0;i<current_lvl_BBs;i++){
	    if(type==3 && i) // overlapping cyclic BBs
	        ell--;

	    std::vector<int> indices;
	    int BB_size = (type==0)?1:(type==2)?6:5;
	    for(int j=0;j<BB_size;j++){
	        if(type == 3 && i == current_lvl_BBs - 1 && j == BB_size - 1){ // for same level, last bit of the last BB shoulf circle back
		    indices.push_back(current_lvl_first_id);
		}
		else
		    indices.push_back(ell++);
	    }

	    BB* bb;
	    if(enable_zeromax)
		bb = new BB(type,indices,zeromax_coeff);
	    else
	        bb = new BB(type,indices);
	    
	    bbs.push_back(bb);
	    //bbs.push_back(BB(type,indices));
	    maxf += 1;

	    if(current_lvl){
	        if(m==1){ // 1-to-many
		    //printf("dep : %d\n",prev_lvl_BB_start_index + i/n);
	            bbs.back()->dependent_on.push_back(bbs[prev_lvl_BB_start_index + i/n]);
	        }
	        else if(n==1){ // many-to-1
		    for(int tmp=0;tmp<m;tmp++)
	                bbs.back()->dependent_on.push_back(bbs[prev_lvl_BB_start_index + i*m + tmp]);
	        }
	    }
        }

	prev_lvl_BB_start_index = bbs.size() - current_lvl_BBs;
	prev_lvl_BBs = current_lvl_BBs;

	if(m<n)
	    current_lvl_BBs *= n;
	else
	    current_lvl_BBs /= m;
    }

    return ell;
}


double LEADINGinstance::evaluate(char* c){
    double fitness = 0;
    for(auto bb: bbs){
	//printf("%f\n",bb->evaluate(c));
        fitness += bb->evaluate(c);
    }
    return fitness;
}

