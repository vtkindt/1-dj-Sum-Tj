// An implementation of memory using single list. It corresponds to the fastest one, without considering memory full strategies.
/*
* This file define main structures to perform memorization. Mainly the Mem class and the SolvedPb class.
* The later represents one node, stored in the "memory".
*/
#ifndef MEM0_H
#define MEM0_H
#include <forward_list>
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
#ifdef MEM_LOG_ML	
	int t0;		// Also log the t0 for training purpose
#endif
	int tt;
#ifdef MEM_FIELD_ID
	unsigned int id;
#endif
#ifdef MEM_FIELD_HIT
	short nbHit;
	SolvedPb() : sol(NULL), nbHit(0){}
#else
	SolvedPb() : sol(NULL){}
#endif
	static vector<Job>* ptrCurrIns;
	//SolvedPb(SolvedPb & sp){ this->sol = sp.sol; this->tt = sp.tt; this->n = sp.n; }
	// std::list need a copy constructor with const arg...
	//SolvedPb(const SolvedPb & sp){ this->sol = (const_cast<SolvedPb &> (sp)).sol; this->tt = sp.tt; }//this->n = sp.n;
	//SolvedPb(int ttt, short* soll) :tt(ttt), sol(soll){}
	template<typename JOBIDS>
	SolvedPb(int ttt, JOBIDS soll, short nn) :tt(ttt)//, n(nn)
	{
#ifdef MEM_FIELD_HIT
		nbHit = 0;
#endif
		// In mem.h we should use native malloc and free instead of Alloc class
		sol = ALLOC0<short>(nn);
		//sol.alloc(nn);
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
	typedef unordered_map<unsigned int, forward_list<SolvedPb>> MapPb;
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
		if (it == mapNbJobs[nbjobs].end() || it->second.empty()) {
			return NULL;
		}
		for (auto itList = it->second.begin(); itList != it->second.end(); ++itList)
			if (jobSetEqual(jobs, itList->sol, nbjobs)){
				++hit;
#ifdef MEM_FIELD_HIT
				++ itList->nbHit;
#endif
				if (itList->tt < 0)itList->tt = -itList->tt;
				return &(*itList);
			}
		++miss;
		return NULL;
	}

	template<typename JOBIDS> //JOBIDS must be Job* or short *, both are used to access the id of the list of jobs
	void addSol(short nbjobs, int startt, int sump, int tt, JOBIDS sol, unsigned int id = -1){
		if (!isNInRange(nbjobs))return;
		unsigned int key = h(startt, sump);
		long long ram = get_ram_usage();
		// We try to limit the nb of jobs (ids in "short") stored in the memory
		if ((!isMemLimitReached && ram > (Config::RAM_LIM_MO * _1M)) || (isMemLimitReached && bytesMem > DB_LIM)
			//|| nbentry==7728436
			) {
			if (!isMemLimitReached){
				cout << "\n ************* Mem limite reached!****************\n" << nbentry << ", " << bytesMem;
				isMemLimitReached = true;
				DB_LIM = bytesMem-1;
				freezeLoadFactor();
			}

			switch (Config::MEM_LIM_STRATEGY){
				case 0: return;
				case 5:{// Replace ou remove all sol with n>N*MEM_RM_RATIO. The ratio is dynamic
					auto itChain = mapNbJobs[nbjobs].find(key);
					//if (itChain != mapNbJobs[nbjobs].end() && !itChain->second.empty())
					//{
					//	replaceLast(itChain->second, nbjobs, sol, tt);
					//	return;
					//}
					//else
					//{
						int nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
						double ratio(Config::MEM_STRA_RATIO), ratioLast(1);
						++Stat::ctrMemRemove;
						long long nbRm(0);
						while ((nbRm += removeByRatio(ratio, ratioLast)) < nbBytesNeeded){
							ratioLast = ratio;
							ratio -= Config::MEM_STRA_RATIO_STEP; // about 12 hashtables.
							if (ratio < MEM_STRA_RATIO_LB){
								Stat::szMemRemove += nbRm;
								++Stat::ctrMemRemoveFail;
								return;
							}
						}
						Stat::szMemRemove += nbRm;
						if (ratio < Stat::minRmRatio) Stat::minRmRatio = ratio;
					//}
					break;
				}
#ifdef MEM_FIELD_HIT
				case 9:{// Remove unused
					int nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
					++Stat::ctrMemRemove;
					long long nbrm = removeUnusedAndUpdateHit();
					Stat::szMemRemove += nbrm;
					if (nbrm < nbBytesNeeded)	return;
					break;
				}
#endif
				default:
					cerr << "Strategy not supported. Pause."<<endl;
					getchar();
			}// End of switch
		}// End if out of memory

		try{
			auto & chain = mapNbJobs[nbjobs][key];
			FAIL_ON(tt < 0);
			chain.emplace_front(-tt, sol, nbjobs);		// Small trick here: use sign to indicate whether it's already queried. Negative: never used
#ifdef MEM_FIELD_ID
			chain.front().id = id;
#endif
#ifdef MEM_LOG_ML
			chain.front().t0 = startt;
#endif
			bytesMem += sizeof(SolvedPb) + nbjobs*(sizeof(short));
			++nbentry;
			//cout << nbentry << endl;
		}
		catch (exception & ex){
			cerr << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			logfile << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			cerr << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			logfile << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			freezeLoadFactor();
			cerr << "Pause";
			getchar();
			//mapNbJobs[nbjobs][key].emplace_front(-tt, sol, nbjobs);
			//bytesMem += sizeof(SolvedPb) + nbjobs*(sizeof(short));
			//++nbentry;
		}
	}

	void freezeLoadFactor(){
		FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS)
			mapNbJobs[i].max_load_factor(9999999.0f);
	}

	bool addSolIfNot(short nbjobs, int startt, int sump, int tt, Job* sol){
		bool isAdded = false;
		auto searchRes = findSol(nbjobs, startt, sump, sol);
		if (searchRes == NULL)
		{
			addSol(nbjobs, startt, sump, tt, sol);
			isAdded = true;
		}
		return isAdded;
	}

	void removeOldDownToLim(bool checkRam ){
		cerr << "removeOldDownToLim not implemented in Mem0. Will call resetMem. Pause." << endl;
		getchar();
		resetMem();
	}

	int removeByRatio(double ratio, double ratioUb = 1){
		if (ratio > ratioUb)	// Circular effect
			return removeByRatio(0, ratioUb) + removeByRatio(ratio, 1);
		int sizeLbToRm = Config::CURR_SIZE_INS * ratio;
		int sizeUbToRm = Config::CURR_SIZE_INS * ratioUb;
		long long oldBytes = bytesMem;
		FOR_BE(i, sizeLbToRm + 1, sizeUbToRm){
			for (auto & kv : mapNbJobs[i]){
				while (!kv.second.empty()){
					kv.second.erase_after(kv.second.before_begin());
					bytesMem -= (sizeof(SolvedPb) + i*sizeof(short));
					-- nbentry;
				}
			}
			mapNbJobs[i].clear();
		}
		return (oldBytes - bytesMem);
	}

#ifdef MEM_FIELD_HIT
	long long removeUnusedAndUpdateHit(){
		long long oldBytes = bytesMem;
		FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS){
			for (auto & kv : mapNbJobs[i]){
				auto itbef = kv.second.before_begin();
				for (auto it = kv.second.begin(); it != kv.second.end();){
					if (it->nbHit-- <= 0){
						++it;
						kv.second.erase_after(itbef);
						--nbentry;
						bytesMem -= (sizeof(SolvedPb) + i*sizeof(short));
					}
					else {
						itbef = it ++ ;
					}
				}
			}
		}
		return oldBytes - bytesMem;
	}
#endif

	void resetMem(){
		if (nbentry > 0){
			for (MapPb & m : mapNbJobs) m.clear();
		}
		hit = miss = nbentry = nbSuperSol = bytesMem = DB_LIM = isMemLimitReached = 0;
		//listOfPbByTime.clear();
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

#ifdef MEM_LOG_ML
	void log(ofstream & o){
		//logChainLen(o);
		logAll(o);
		//logNNbSolHits(o);
	}

	
	void logAll(ofstream & o){
		static auto cmpIdEDD = [](short id1, short id2){
			return cmpJobD0(SolvedPb::ptrCurrIns->at(id1 - 1), SolvedPb::ptrCurrIns->at(id2 - 1));
		};
		for (int n = 6; n < SolvedPb::ptrCurrIns->size(); ++n)
		for (auto & kv : mapNbJobs[n]){
			for (auto & pb : kv.second){
				//if (pb.nbHit >= minhit)
				{
					// Save this sequence
					o << pb.nbHit << " " << n << " " << pb.t0 << endl;
					for (int i = 0; i < n; ++i){
						o<<pb.sol[i]<<" ";
					}

					o<< endl;
				}
			}
		}
	}
#endif

};

#endif