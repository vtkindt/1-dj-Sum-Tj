// A Database for memorizing a list of solved problems 
// (It is different with a database keeping currently best nodes)
// Caution: We do not consider d' version since we do not check d values when matching. But it should be ok if only original d values are used to compute tt
// Use MAX_BYTES_MEM_MO to control the memory size

#ifndef MEM_H
#define MEM_H
#include <list>
#include <unordered_map>
#include "h.h"
using namespace std;
extern ofstream logfile;

//////////////////////////////////////////////////////////////////////////
// A sequence that is stored separately in the heap
//////////////////////////////////////////////////////////////////////////
typedef short *Seq ; 
#define NULL_SEQ (Seq(0)) //or Seq()

//////////////////////////////////////////////////////////////////////////
// Item stored
//////////////////////////////////////////////////////////////////////////
class SolvedPb{
public:
	Seq sol;	// We only store ids
#ifdef MEM_FIELD_ID
	unsigned int id;
#endif
	int tt;
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
		nbHit=0;
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
	typedef unordered_map<unsigned int, list<SolvedPb>> MapPb;
	struct DeletionPoint{
		list<SolvedPb> * l;
		list<SolvedPb>::iterator it;
		short n;
	};
	list< DeletionPoint > listOfPbByTime;
	DeletionPoint * endRmSection;		// Used by stra7. We remove thing on the section from listOfPbByTime.begin to endRmSection
	char rmFlag;						// 0 = remove unused in old; 1 = remove old.
	long long hit, miss, nbentry, bytesMem, DB_LIM;
	int nbSuperSol;
	bool isMemLimitReached;
	
	// The 1st dimension of the hash table: index by nb jobs. Each item in the array is a map.
	MapPb mapNbJobs[MAX_N_JOBS+1];	

	Mem() :hit(0), miss(0), nbentry(0), nbSuperSol(0), bytesMem(0), isMemLimitReached(false), DB_LIM(0){
	};	
	static bool isNInRange(int n){
		if (n > (Config::CURR_SIZE_INS-1) * Config::MEM_N_UB_RATIO || n > Config::MEM_N_UB || n < Config::MEM_N_LB)
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
				++ hit;
#ifdef MEM_FIELD_HIT
				++ itList->nbHit;
#endif
				if (itList->tt < 0)itList->tt = -itList->tt;
				return &(*itList);
			}
		++ miss;
		return NULL;
	}

	template<typename JOBIDS> //JOBIDS must be Job* or short *, both are used to access the id of the list of jobs
	void addSol(short nbjobs, int startt, int sump, int tt, JOBIDS sol, unsigned int idPb){
		if (!isNInRange(nbjobs))return;
		unsigned int key = h(startt, sump);
		long long ram = get_ram_usage();
		// We try to limit the nb of jobs (ids in "short") stored in the memory
		if ((!isMemLimitReached && ram > (Config::RAM_LIM_MO * _1M)) || (isMemLimitReached && bytesMem > DB_LIM) ) {
			static double rationStra8 = 1;
			if (!isMemLimitReached){
				cout << "\n ************* Mem limite reached!****************\n"<< nbentry<<", "<<bytesMem;
				isMemLimitReached = true;
				DB_LIM = bytesMem;
				if (Config::MEM_LIM_STRATEGY == 7){
					endRmSection = &(*(--listOfPbByTime.end()));
					rmFlag = 0;
				}
				rationStra8 = Config::MEM_STRA_RATIO;
				freezeLoadFactor();
			}

			// According to strategies
			switch(Config::MEM_LIM_STRATEGY) {
				case 0: return;
				case 1: {// replace LastUnused Or Last
					auto itChain = mapNbJobs[nbjobs].find(key);
					if (itChain == mapNbJobs[nbjobs].end() || itChain->second.empty())return;//not inserted. Note: if use chain=map[][], a empty chain will be inserted which explains the mem increase
					replaceLast(itChain->second, nbjobs, sol, tt);
					return;
				}
				case 2: {// Replace LastUnused Or Last, otherwise remove global old entries until got enough space.
					auto itChain = mapNbJobs[nbjobs].find(key);
					if (itChain != mapNbJobs[nbjobs].end() && !itChain->second.empty())
					{
						//cerr << "Replace happened!\n"; getchar();
						replaceLast(itChain->second, nbjobs, sol, tt);
						return;
					}
					else{// Cannot remove anything in the current chain (not exist or empty), remove global old pbs
						long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
						if(!removeOldDownToLim(false, nbBytesNeeded))return;
						break;
					}// Else: Remove
				}//Case 2
				case 3: {	// Directly remove old entries
					cerr << "Stra 3 Deprecated." << endl; getchar();
					removeOldDownToLim(false);
					break;
				}//Case 3
				case 4:{// Replace ou remove all sol with n>N*MEM_RM_RATIO
					auto itChain = mapNbJobs[nbjobs].find(key);
					if (itChain != mapNbJobs[nbjobs].end() && !itChain->second.empty())
					{
						//cerr << "Replace happened!\n"; getchar();
						replaceLast(itChain->second, nbjobs, sol, tt);
						return;
					}
					else{
						++Stat::ctrMemRemove;
						long long nbRm;
						if ( (nbRm = removeByRatio(Config::MEM_STRA_RATIO)) == 0)
						{
							++Stat::ctrMemRemoveFail;
							return;
						}
						Stat::szMemRemove += nbRm;
					}
					break;
				}
				case 5:{// Replace ou remove all sol with n>N*MEM_RM_RATIO. The ratio is dynamic
					auto itChain = mapNbJobs[nbjobs].find(key);
					if (itChain != mapNbJobs[nbjobs].end() && !itChain->second.empty())
					{
						replaceLast(itChain->second, nbjobs, sol, tt);
						return;
					}
					else{
						long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
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
					}
					break;
				}
				case 6:{// Replace ou remove big sols
					auto itChain = mapNbJobs[nbjobs].find(key);
					if (itChain != mapNbJobs[nbjobs].end() && !itChain->second.empty())
					{
						replaceLast(itChain->second, nbjobs, sol, tt);
						return;
					}
					else{
						++Stat::ctrMemRemove;
						long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
						long long nbRm(0);
						for (int i = Config::CURR_SIZE_INS - 1; i > Config::CURR_SIZE_INS*MEM_STRA_RATIO_LB; --i){
							nbRm += removeByN(i);
							if (nbRm > nbBytesNeeded){
								double ratio = (double)i / Config::CURR_SIZE_INS;
								//cout << i << "\t" << ratio<<endl;
								if (ratio < Stat::minRmRatio) Stat::minRmRatio = ratio;
								break;
							}
						}
						Stat::szMemRemove += nbRm;
						if (nbRm < nbBytesNeeded){
							++Stat::ctrMemRemoveFail;
							return;
						}
					}
					break;
				}
				case 7:{// Consider the current global list without update, remove unused among old, then old. Then take the new current global list and to same thing.
					++Stat::ctrMemRemove;
					long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
					long long nbRm(0);
					for (auto itDelList = listOfPbByTime.begin(); &(*itDelList) != endRmSection; ){
						if (rmFlag == 1 || itDelList->it->tt < 0){
							cerr << "Enable SolvedPb::n and the following line first.Pause" << endl;
							//nbRm += (sizeof(SolvedPb) + (itDelList->it->n)*sizeof(short));
							--nbentry;
							itDelList->l->erase(itDelList->it);
							listOfPbByTime.erase(itDelList++);
						}
						else ++itDelList;
						if (&(*itDelList) == endRmSection){
							if (rmFlag == 0){
								rmFlag = 1;
								itDelList = listOfPbByTime.begin();
							}
							else endRmSection = &(*(--listOfPbByTime.end()));
						}
						if (nbRm > nbBytesNeeded)break;
					}
					bytesMem -= nbRm;
					Stat::szMemRemove += nbRm;
					if (nbRm < nbBytesNeeded)
					{// Failed
						++Stat::ctrMemRemoveFail;
						cerr << "Stra 7, nothing to remove. Not normal. " << nbentry<<"\t"<< bytesMem << "\t" << ram <<"\t"<< listOfPbByTime.size()<< endl;
						getchar();
						return;
					}
					break;
				}
				case 8:{// Like stra5 but each time remove sequences on a range of sizes (defined by ratio step)
					auto itChain = mapNbJobs[nbjobs].find(key);
					{
						long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
						++Stat::ctrMemRemove;
						long long nbRm(0);
						double newRatio = (rationStra8-Config::MEM_STRA_RATIO_STEP);
						if (newRatio <= 0)newRatio += 1;
						while ((nbRm += removeByRatio(newRatio, rationStra8)) < nbBytesNeeded){
							rationStra8 = newRatio;
							newRatio = (rationStra8 - Config::MEM_STRA_RATIO_STEP);
							if (newRatio <= 0)newRatio += 1;
						}
						Stat::szMemRemove += nbRm;
						//if (ratio < Stat::minRmRatio) Stat::minRmRatio = ratio;
					}
					break;
				}
#ifdef MEM_FIELD_HIT
				case 9: {// remove all Unused
					long long nbBytesNeeded = (sizeof(SolvedPb) + (sizeof(short)*nbjobs))* MEM_STRA_NEED_FACTOR;
					++Stat::ctrMemRemove;
					long long nbrm = removeUnusedAndUpdateHit();
					Stat::szMemRemove += nbrm;
					if (nbrm<nbBytesNeeded)	return;
					break;
				}
#endif
				default:
					cerr << "Strategy not supported. Pause." << endl;
					getchar();
			}//Switch
		}// End if out of memory

		try{
			auto & chain = mapNbJobs[nbjobs][key];
			chain.emplace_front(-tt, sol, nbjobs);		// Small trick here: use sign to indicate whether it's already queried. Negative: never used
#ifdef MEM_FIELD_ID
			chain.front().id = idPb;
#endif
			// Now we got space: add the sol
			bytesMem +=  sizeof(SolvedPb) + (nbjobs)* (sizeof(short)); //count actural allocation. (sizeof(short)*nbjobs) +
			++nbentry;
			if (Config::MEM_LIM_STRATEGY == 2 || Config::MEM_LIM_STRATEGY==3 || Config::MEM_LIM_STRATEGY==7)
				listOfPbByTime.push_back({ &chain, chain.begin(), nbjobs });
		}
		catch (exception & ex){
			cerr << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			logfile << __FILE__ << ":" << __LINE__ << ":" << ex.what() << endl;
			cerr << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			logfile << nbjobs << "\t" << mapNbJobs[nbjobs].size() << "\t" << mapNbJobs[nbjobs].bucket_count() << "\t" << nbentry << "\t" << bytesMem << endl;
			freezeLoadFactor();
			cerr << "Pause";
			getchar();
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
			addSol(nbjobs, startt, sump, tt, sol,1);
			isAdded = true;
		}
		return isAdded;
	}

	void resetMem(){
		if (nbentry > 0){
			for (MapPb & m : mapNbJobs) m.clear();
		}
		hit = miss = nbentry = nbSuperSol = bytesMem = DB_LIM = isMemLimitReached=0;
		listOfPbByTime.clear();
	}

	template<typename JOBIDS>
	bool replaceLast(list<SolvedPb> & chain, int n, JOBIDS sol, int tt){
		if (chain.empty()) return false;
		auto curr = --chain.end();
		if (curr != chain.begin()){
			chain.splice(chain.begin(), chain, curr);
			curr = chain.begin();
		}
		// Replace "curr"
		FOR_E(i, n) curr->sol[i] = ID(sol[i]);
		curr->tt = -tt;
		++Stat::ctrMemReplace;
		return true;
	}

	template<typename JOBIDS>
	bool replaceLastUnusedOrLast(list<SolvedPb> & chain, int n, JOBIDS sol, int tt){
		if (chain.empty()) return false;
		auto curr = --chain.end();
		while (curr != chain.begin()){
			if (curr->tt < 0)break;
			--curr;
		}
		if (curr != chain.begin()) {
			chain.splice(chain.begin(), chain, curr);
			curr = chain.begin();
		}
		// Replace "curr"
		FOR_E(i, n) curr->sol[i] = ID(sol[i]);
		curr->tt = -tt;
		++Stat::ctrMemReplace;
		return true;
	}

	// When malloc fails, this function may be called to free some space
	bool removeOldDownToLim(bool checkRam, long long nbBytesNeeded = 0 ){
		auto itDelList = listOfPbByTime.begin();
		long long oldBytes = bytesMem;
		while (itDelList != listOfPbByTime.end()){
			if (checkRam && get_ram_usage() < (Config::RAM_LIM_MO * _1M-nbBytesNeeded))break;
			if (!checkRam && bytesMem < DB_LIM - nbBytesNeeded) break;

			//cerr << "Enable SolvedPb::n and the following line first.Pause"<<endl;
			bytesMem -= (sizeof(SolvedPb) + (itDelList->n)*sizeof(short));
			--nbentry;
			//if (nbentry < 1000)cout << nbentry << "\t" << bytesMem << endl;
			itDelList->l->erase(itDelList->it);
			listOfPbByTime.erase(itDelList++);
		}
		++Stat::ctrMemRemove;
		Stat::szMemRemove += (oldBytes - bytesMem);
		if (listOfPbByTime.empty())
		{
			++ Stat::ctrMemRemoveFail;
			FAIL_ON(!string("Global list empty, no entry to remove.").empty());
			return false;
		}
		return true;
	}
#ifdef MEM_FIELD_HIT
	int removeUnusedAndUpdateHit(){
		long long oldBytes = bytesMem;
		FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS){
			for (auto & kv : mapNbJobs[i])
				for (auto it = kv.second.begin(); it != kv.second.end();){
					if (it->nbHit-- <= 0){
						kv.second.erase(it++);
						--nbentry;
						bytesMem -= (sizeof(SolvedPb) + i*sizeof(short));
					}
					else ++it;
				}
		}
		return oldBytes - bytesMem;
	}
#endif

	int removeByRatio(double ratio, double ratioUb = 1){
		if (ratio > ratioUb)	// Circular effect
			return removeByRatio(0, ratioUb) + removeByRatio(ratio,1);
		int sizeLbToRm = Config::CURR_SIZE_INS * ratio;
		int sizeUbToRm = Config::CURR_SIZE_INS * ratioUb;
		long long oldBytes = bytesMem;
		FOR_BE(i, sizeLbToRm + 1, sizeUbToRm){
			for (auto & kv : mapNbJobs[i]){
				bytesMem -= kv.second.size()*(sizeof(SolvedPb) + i*sizeof(short));
				nbentry -= kv.second.size();
			}
			mapNbJobs[i].clear();
		}
		return (oldBytes - bytesMem);
	}

	int removeByN(int n){
		long long oldBytes = bytesMem;
		for (auto & kv : mapNbJobs[n]){
			bytesMem -= kv.second.size()*(sizeof(SolvedPb) + n*sizeof(short));
			nbentry -= kv.second.size();
		}
		mapNbJobs[n].clear();
		return (oldBytes - bytesMem);
	}
	// I tried to use static var and reset flags before returning, but the time is the same, surprisingly
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


#ifdef MEM_LOG
	void log(ofstream & o){
		//logChainLen(o);
		logUsedSol(o, 2);
		//logNNbSolHits(o);
	}

	void logNbEntries(ofstream & o){
		// try to see the distribution of entries: for each N, there are how many entries? which of them are queried?
		for (int i = Config::MEM_N_LB; i < Config::CURR_SIZE_INS; ++i){
			long long entries = 0;
			long long usefulEntries = 0;
			for (auto & kv : mapNbJobs[i])
				for (auto & pb: kv.second)
				{
					++ entries;
					if (pb.nbHit > 0) ++usefulEntries;
				}
			o << i << "\t" << entries <<"\t"<<usefulEntries << "\n";
			//for (aut)
		}
	}

	void logChainLen(ofstream & o){
		// How long are the chains?
		map<int, int> len_ctr;
		for (int i = Config::MEM_N_LB; i < Config::CURR_SIZE_INS; ++i){
			for (auto & kv : mapNbJobs[i]){
				++len_ctr[int(kv.second.size())];
			}
		}
		for (auto & ctr : len_ctr)
			o << ctr.first << "\t" << ctr.second << endl;
	}
	void logNNbSolHits(ofstream& o){
		// Log for given N, the number of solution that are hit at minimum minhit times
		for (int i = Config::MEM_N_LB; i < Config::CURR_SIZE_INS; ++i){
			int ctr1(0), ctr2(0), ctr3(0), ctr4(0);
			for (auto & kv : mapNbJobs[i]){
				for (auto & pb : kv.second)
					if (pb.nbHit >= 4)++ctr1, ++ctr2, ++ctr3, ++ctr4;
					else if (pb.nbHit >= 3)++ctr2, ++ctr3, ++ctr1;
					else if (pb.nbHit >= 2)++ctr1, ++ctr2;
					else if (pb.nbHit >= 1)++ctr1;
			}
			o << i << "\t" << ctr1 << "\t" << ctr2 << "\t" << ctr3 << "\t" << ctr4 << endl;
		}
	}
	void logUsedSol(ofstream & o, int minhit = 1){
		static auto cmpIdEDD = [](short id1, short id2){
			return cmpJobD0(SolvedPb::ptrCurrIns->at(id1 - 1), SolvedPb::ptrCurrIns->at(id2-1)); 
		};
		int n = 50;
		// Sort global to edd
		vector<Job> currIns = *SolvedPb::ptrCurrIns;
		sort(currIns.begin(), currIns.end(), cmpJobD0);
		Job* jobs = new Job[n];
		for (auto & kv : mapNbJobs[n]){
			for (auto & pb : kv.second){
				if (pb.nbHit >= minhit){
					// Save this sequence
					// Frist sort to edd
					sort(pb.sol, pb.sol + n, cmpIdEDD);
					// Save it showing its position in global edd sequence
					for (int i = 0, j = 0; i < currIns.size(); ++i){
						if (j == n || currIns[i].id != pb.sol[j])
							o << setw(3) << currIns[i].id << "\t";
						else if (currIns[i].id == pb.sol[j]){
							o << '#' << "\t";
							++j;
						}
					}
					o << endl;
				}
			}
		}
		delete []jobs;
	}
#endif

};

#endif