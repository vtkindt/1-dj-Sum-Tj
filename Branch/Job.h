#ifndef JOB_H
#define JOB_H
#include <iostream>
#include <bitset>
using namespace std;

/////////////////////////
// Class Job
/////////////////////////
class Problem;
class Job{
public:
    short id;
    short p;
    int d;
    static const Job* jobscmp;		// Cmp functions compare job index based on these jobs. 
									// !! To be set before each comparison
									// !! So multi-threading is not supported!
    
    Job(){};						// Used implicitly by std container, should never be called by user
    Job(short ii, short pp, int dd) :id(ii), p(pp), d(dd){}//, idEdd(-1)
	int T(int C)const{
		return (C - d) > 0 ? (C - d) : 0;
	}
};
inline short ID(const Job& j){ return j.id; }
inline short ID(short id){ return id; }

bool cmpJobP0(const Job &j1, const Job &j2);
bool cmpJobD0(const Job & j1, const Job& j2);

// !!! Set Job::jobscmp before each comparison
bool cmpJobPIndex(short j1, short j2);
bool cmpJobDIndex(short j1, short j2);
bool cmpJobPDIndex(short j1, short j2);

ostream & operator << (ostream & o, const Job& j);
bool operator== (const Job& j1, const Job& j2);
bool operator!= (const Job& j1, const Job& j2);

typedef Job* PtrJob;
typedef const Job* PtrJobConst;
typedef short JobIndex;
typedef short* PtrJobIndex;
typedef const short* PtrJobIndexConst;
typedef short Pos;
typedef short* PtrPos;
string strIdList(PtrJob p, short n);
string strPoList(PtrPos p, short n);
#endif
