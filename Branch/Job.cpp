#include "Problem.h"
#include "Job.h"
#include <iostream>
const Job* Job::jobscmp = NULL;
// Functors used for sorting sequences
bool cmpJobP0(const Job &j1, const Job &j2)
{
    if (j1.p != j2.p) return j1.p < j2.p;
    if (j1.d != j2.d) return j1.d < j2.d;
    return j1.id < j2.id;
}
bool cmpJobD0(const Job &j1, const Job& j2)
{
	if (j1.d != j2.d) return j1.d < j2.d;
	if (j1.p != j2.p) return j1.p < j2.p;
	return j1.id < j2.id;
}
bool cmpJobD(const Job* j1, const Job* j2)
{
    if (j1->d != j2->d) return j1->d < j2->d;
    if (j1->p != j2->p) return j1->p < j2->p;
    return j1->id < j2->id;
}
bool cmpJobPD(const Job* j1, const Job* j2){ return j1->p + j1->d < j2->p + j2->d; }

// !!!! REMEBER to set jobscmp before calling !!!!!
bool cmpJobPIndex(short j1, short j2){ return cmpJobP0(Job::jobscmp[j1], Job::jobscmp[j2]); }
bool cmpJobDIndex(short j1, short j2){ return cmpJobD(Job::jobscmp + j1, Job::jobscmp + j2); }
bool cmpJobPDIndex(short j1, short j2){ return cmpJobPD(Job::jobscmp + j1, Job::jobscmp + j2); }

// Support printing job id
ostream & operator << (ostream & o, const Job& j){
    o << "("<<j.id<<","<<j.p<<","<<j.d<<")";
    return o;
}
bool operator== (const Job& j1, const Job& j2){ return j1.id == j2.id; }; 
bool operator!= (const Job& j1, const Job& j2){ return j1.id != j2.id; };

string strIdList(PtrJob p, short n){							   //log job/sol list
	cout << "strIdList called!" << endl;
	stringstream ss;
	ss << "[";
	for_each(p, p + n, [&ss](Job& p){
		ss << p.id << ",";
	});
	ss << "]";
	return ss.str();
}
string strPoList(PtrPos p, short n){							   //log job/sol list
	stringstream ss;
	ss << "[";
	for_each(p, p + n, [&ss](Pos p){
		ss << p << ",";
	});
	ss << "]";
	return ss.str();
}