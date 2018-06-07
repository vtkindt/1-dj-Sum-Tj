// An implementation of memory without list

#ifndef MEM3_H
#define MEM3_H
#include <unordered_map>
#include "h.h"
using namespace std;
extern ofstream logfile;

//////////////////////////////////////////////////////////////////////////
// A sequence that is stored separately in the heap
//////////////////////////////////////////////////////////////////////////
typedef short *Seq;
#define NULL_SEQ (Seq(0)) //or Seq()

//////////////////////////////////////////////////////////////////////////
// Item stored
//////////////////////////////////////////////////////////////////////////
class SolvedPb{
public:
	Seq sol;	// We only store ids
	int tt;
	//short n;
#ifdef MEM_LOG
	short nbHit;
	SolvedPb() :tt(INT_MIN)//, n(0)
		, nbHit(0){}
#else
	SolvedPb() : sol(NULL)//, n(0)
	{}
#endif
	static vector<Job>* ptrCurrIns;
	//SolvedPb(SolvedPb & sp){ this->sol = sp.sol; this->tt = sp.tt; this->n = sp.n; }
	// std::list need a copy constructor with const arg...
	//SolvedPb(const SolvedPb & sp){ this->sol = (const_cast<SolvedPb &> (sp)).sol; this->tt = sp.tt; }//this->n = sp.n;
	//SolvedPb(int ttt, short* soll) :tt(ttt), sol(soll){}
	template<typename JOBIDS>
	SolvedPb(int ttt, JOBIDS soll, short nn) :tt(ttt)//, n(nn)
	{
#ifdef MEM_LOG
		nbHit = 0;
#endif
		// In mem.h we should use native malloc and free instead of Alloc class
		sol = ALLOC0<short>(nn);
		//sol.alloc(nn);
		FOR_E(i, nn) sol[i] = ID(soll[i]);
	}

	// Called after default constructor
	template<typename JOBIDS>
	void setup(int ttt, JOBIDS soll, short nn){
		FREE0<short>(&sol);
		sol = ALLOC0<short>(nn);
		tt = ttt;
		FOR_E(i, nn) sol[i] = ID(soll[i]);
	}
	// Get the job of sol[index]
	const Job & operator [](int index){
		return (*ptrCurrIns)[sol[index] - 1];
	}
	~SolvedPb(){
		//sol.free();
		FREE0<short>(&sol);
	}
};

//////////////////////////////////////////////////////////////////////////
// Memory/Database
//////////////////////////////////////////////////////////////////////////
class Mem{
public:
	// The key will be given by (startT*MAX_P_SUM+P). Sizeof(list<SolvedPb>)=16
	typedef unordered_map<unsigned int, SolvedPb> MapPb;
	struct DeletionPoint{
		short n;
		MapPb::iterator it;
	};
	list< DeletionPoint > listOfPbByTime;
	DeletionPoint * endRmSection;		// Used by stra7. We remove thing on the section from listOfPbByTime.begin to endRmSection
	char rmFlag;						// 0 = remove unused in old; 1 = remove old.
	long long hit, miss, nbentry, bytesMem, DB_LIM;
	int nbSuperSol;
	bool isMemLimitReached;

	// The 1st dimension of the hash table: index by nb jobs. Each item in the array is a map.
	MapPb mapNbJobs[MAX_N_JOBS + 1];

	Mem() :hit(0), miss(0), nbentry(0), nbSuperSol(0), bytesMem(0), isMemLimitReached(false), DB_LIM(0){
	};
	static bool isNInRange(int n){
		if (n > (Config::CURR_SIZE_INS - 1) * Config::MEM_N_UB_RATIO || n > Config::MEM_N_UB || n < Config::MEM_N_LB)
			return false;
		return true;
	}
	// Return sol if found, null otherwise. The caller should make a copy of sol. The returned sol should not be freed by the caller
	//!! Must return by reference otherwise moved
	template<typename JOBIDS>
	SolvedPb* findSol(short nbjobs, int startt, int sump, JOBIDS jobs){
		if (!isNInRange(nbjobs))return NULL;
		auto it = mapNbJobs[nbjobs].find(h(startt, sump));
		if (it == mapNbJobs[nbjobs].end()) {
			return NULL;
		}
		if (jobSetEqual(jobs, it->second.sol, nbjobs)){
			++hit;
#ifdef MEM_LOG
				++ itList->nbHit;
#endif
				if (it->second.tt < 0) it->second.tt = -it->second.tt;
				return &(it->second);
		}
		++miss;
		return NULL;
	}

	template<typename JOBIDS> //JOBIDS must be Job* or short *, both are used to access the id of the list of jobs
	void addSol(short nbjobs, int startt, int sump, int tt, JOBIDS sol){
		if (!isNInRange(nbjobs))return;
		unsigned int key = h(startt, sump);
		long long ram = get_ram_usage();
		// We try to limit the nb of jobs (ids in "short") stored in the memory
		if ((!isMemLimitReached && ram > (Config::RAM_LIM_MO * _1M)) || (isMemLimitReached && bytesMem > DB_LIM)) {
			if (!isMemLimitReached){
				cout << "\n ************* Mem limite reached!****************\n" << nbentry << ", " << bytesMem;
				isMemLimitReached = true;
				DB_LIM = bytesMem;
				if (Config::MEM_LIM_STRATEGY == 7){
					endRmSection = &(*(--listOfPbByTime.end()));
					rmFlag = 0;
				}
				freezeLoadFactor();
				cerr << "No strategy implemented. Pause." << endl;
				getchar();
			}

		}// End if out of memory

		try{
			auto it = mapNbJobs[nbjobs].find(key); 	// Small trick here: use sign to indicate whether it's already queried. Negative: never used
			if (it == mapNbJobs[nbjobs].end()){
				auto & pb = mapNbJobs[nbjobs][key];
				pb.setup(-tt, sol, nbjobs);
			}
			else{//replace since the n is same
				it->second.tt = -tt;
				FOR_E(i, nbjobs)it->second.sol[i] = ID(sol[i]);
				++ Stat::ctrMemReplace;
			}
			++nbentry;
		}
		catch (exception & ex){
			cerr << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			logfile << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			cerr << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			logfile << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			// Augment balance factor
			//mapNbJobs[nbjobs].max_load_factor(mapNbJobs[nbjobs].load_factor() + 1);
			freezeLoadFactor();
			auto & pb = mapNbJobs[nbjobs][key];
			pb.setup(-tt, sol, nbjobs);
			++nbentry;
		}
	}

	void freezeLoadFactor(){
		FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS)
			mapNbJobs[i].max_load_factor(9999999.0f);
	}

	void removeOldDownToLim(bool checkRam){
		cerr << "removeOldDownToLim not implemented in Mem0. Will call resetMem. Pause." << endl;
		getchar();
		resetMem();
	}

	void resetMem(){
		if (nbentry > 0){
			for (MapPb & m : mapNbJobs) m.clear();
		}
		hit = miss = nbentry = nbSuperSol = bytesMem = DB_LIM = isMemLimitReached = 0;
		listOfPbByTime.clear();
	}

	template<typename JOBIDS>
	bool jobSetEqual(const JOBIDS js1, const Seq js2ids, int n){
		bool tester[MAX_N_JOBS + 1] = { 0 };
		assert(tester[MAX_N_JOBS] == 0);
		FOR_E(i, n)tester[ID(js1[i])] = true;
		FOR_E(i, n)if (tester[js2ids[i]] == false) return false;
		return true;
	}

	template<typename JOBIDS>
	bool jobSeqEqual(const JOBIDS js1, const Seq js2ids, int n){
		FOR_E(i, n)if (ID(js1[i]) != js2ids[i]) return false;
		return true;
	}

	// hash function
	unsigned int h(unsigned int startt, unsigned int sump){ return (startt)*MAX_P_SUM + sump; }
	void h_reverse(unsigned int hashvalue, int & t0, int & sump){
		sump = hashvalue % MAX_P_SUM;
		t0 = (hashvalue - sump) / MAX_P_SUM;
	}

	// Not efficient, only used for logging
	string hashJobSeq(short* sol, int n){
		stringstream ss;
		FOR_E(i, n)ss << sol[i] << ",";
		return ss.str();
	}

};

#endif