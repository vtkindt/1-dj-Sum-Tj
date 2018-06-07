// Mem2 without list, where SolvedPb = {sol, map<t0,tt>}

#ifndef MEM4_H
#define MEM4_H
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
// Item stored: a solution sequence and some <t0,tt> associated
//////////////////////////////////////////////////////////////////////////
class SolvedPb{
public:
	Seq sol;	// We only store ids
	map<int, int> t0_tt;
#ifdef MEM_LOG
	short nbHit;
#else
	SolvedPb() : sol(NULL){}
#endif
	static vector<Job>* ptrCurrIns;
	template<typename JOBIDS>
	SolvedPb(int t00, int ttt, JOBIDS soll, short nn) //:n(nn)
	{
#ifdef MEM_LOG
		nbHit = 0;
#endif
		// In mem.h we should use native malloc and free instead of Alloc class
		sol = ALLOC0<short>(nn);
		//sol.alloc(nn);
		FOR_E(i, nn) sol[i] = ID(soll[i]);
		t0_tt[t00] = ttt;
	}

	template<typename JOBIDS>
	void setup(int t00, int ttt, JOBIDS soll, short nn){
		FREE0<short>(&sol);
		sol = ALLOC0<short>(nn);
		//sol.alloc(nn);
		FOR_E(i, nn) sol[i] = ID(soll[i]);
		t0_tt[t00] = ttt;
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
	// The key will be given by P. We may also add sum(ID)
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
		auto it = mapNbJobs[nbjobs].find(h(jobs, nbjobs, sump));
		if (it == mapNbJobs[nbjobs].end()) {
			return NULL;
		}
		auto ittt = it->second.t0_tt.find(startt);
		if (ittt != it->second.t0_tt.end()) {
			if ((jobSetEqual(jobs, it->second.sol, nbjobs)))
			{
				++hit;
#ifdef MEM_LOG
				++ itList->nbHit;
#endif
				if (ittt->second < 0)ittt->second = -ittt->second;
				return &(it->second);
			}
			++miss;
		}
		return NULL;
	}

	int byteNewSol(int n){ return sizeof(SolvedPb) + n*sizeof(short) + 2 * sizeof(int); }
	int byteSol(const SolvedPb & pb){ return int(sizeof(SolvedPb) + 50 * sizeof(short) + pb.t0_tt.size() * 2 * sizeof(int)); cerr << "byteSol() not supported!"; }

	template<typename JOBIDS> //JOBIDS must be Job* or short *, both are used to access the id of the list of jobs
	void addSol(short nbjobs, int startt, int sump, int tt, JOBIDS sol){
		if (!isNInRange(nbjobs))return;
		unsigned int key = h(sol, nbjobs, sump);
		long long ram = get_ram_usage();
		// We try to limit the nb of jobs (ids in "short") stored in the memory
		if ((!isMemLimitReached && ram > (Config::RAM_LIM_MO * _1M)) || (isMemLimitReached && bytesMem > DB_LIM)) {
			if (!isMemLimitReached){
				cout << "\n ************* Mem limite reached!****************\n" << nbentry << ", " << bytesMem;
				isMemLimitReached = true;
				freezeLoadFactor();
				cerr << "No strategy supported. Pause. \n" << endl;
				getchar();
			}
		}// End if out of memory

		try{
			auto it = mapNbJobs[nbjobs].find(key);
			if (it != mapNbJobs[nbjobs].end()){
				int i = 0;
				for (; i < nbjobs; ++i)if (ID(sol[i]) != it->second.sol[i])return;
				if (i==nbjobs){
					it->second.t0_tt[startt] = -tt;
					bytesMem += 2 * sizeof(int);
					return;
				}
				else{
					for (; i < nbjobs; ++i)it->second.sol[i] = ID(sol[i]);
					it->second.t0_tt.clear();
					it->second.t0_tt[startt] = -tt;
					++ Stat::ctrMemReplace;
					return;
				}
			}
			if (nbentry % 10000 == 0)cout << nbentry << endl;
			// New seq
			auto & pb = mapNbJobs[nbjobs][key];
			pb.setup(startt, -tt, sol, nbjobs);		// Small trick here: use sign to indicate whether it's already queried. Negative: never used
			// Now we got space: add the sol
			bytesMem += byteNewSol(nbjobs);
			++nbentry;
			if (Config::MEM_LIM_STRATEGY == 2 || Config::MEM_LIM_STRATEGY == 3 || Config::MEM_LIM_STRATEGY == 7)
				listOfPbByTime.push_back({ nbjobs, mapNbJobs[nbjobs].find(key) });
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
			pb.setup(startt, -tt, sol, nbjobs);		// Small trick here: use sign to indicate whether it's already queried. Negative: never used
			// Now we got space: add the sol
			bytesMem += byteNewSol(nbjobs);
			++nbentry;
		}
	}

	bool removeOldDownToLim(bool checkRam){
		auto itDelList = listOfPbByTime.begin();
		long long oldBytes = bytesMem;
		while (itDelList != listOfPbByTime.end()){
			if (checkRam && get_ram_usage() < (Config::RAM_LIM_MO * _1M))break;
			if (!checkRam && bytesMem < DB_LIM) break;

			bytesMem -= byteSol(itDelList->it->second);
			--nbentry;
			mapNbJobs[itDelList->n].erase(itDelList->it);
			listOfPbByTime.erase(itDelList++);
		}
		++Stat::ctrMemRemove;
		Stat::szMemRemove += (oldBytes - bytesMem);
		if (listOfPbByTime.empty())
		{
			++Stat::ctrMemRemoveFail;
			return false;
		}
		return true;
	}

	void freezeLoadFactor(){
		FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS)
			mapNbJobs[i].max_load_factor(9999999.0f);
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
	
	// hash function
	template<typename JOBIDS>
	unsigned int h(JOBIDS jobs, int n, unsigned int sump){
		//return sump;
		int sumid = 0;
		FOR_E(i, n)sumid += ID(jobs[i]);
		return (sumid)*MAX_P_SUM + sump;
	}

	// Not efficient, only used for logging
	string hashJobSeq(short* sol, int n){
		stringstream ss;
		FOR_E(i, n)ss << sol[i] << ",";
		return ss.str();
	}

	void log(ofstream & o){
		FOR_E(i, Config::CURR_SIZE_INS){
			o << endl << i << endl;
			for (auto & kv : mapNbJobs[i]){
				o << kv.first << "\t" << kv.second.t0_tt.size() << endl;
			}
		}
	}
};

#endif