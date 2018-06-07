#ifndef PROBLEM_H
#define PROBLEM_H

/**
* This file defines Problem class, including main branching algorithms.
*/
#include <memory>
#include "h.h"
#include  "Job.h"

#ifdef USE_MEM0
#include "Mem0.h"
#elif defined USE_MEM2
#include "Mem2.h"
#elif defined USE_MEM3
#include "Mem3.h"
#elif defined USE_MEM4
#include "Mem4.h"
#else
#include "Mem.h"
#endif
#define INLINE inline
#pragma warning(disable:4127)	// For FAIL_ON(true): l'expression conditionnelle est une constante

#define EDD(i) (jobs[edd[i]])

inline void msort(PtrJobConst jobs, PtrJobIndex edd, int n, bool func(JobIndex,JobIndex) ) 
{ Job::jobscmp = jobs; sort(edd, edd + n, func); }

inline int sump(PtrJobConst jobs, int n, int init = 0)
{
    while (n-- > 0)init += (jobs++)->p;
    return init;
}

// Helper functions for logging purpose
class Problem;
string strPoList(PtrPos p, short n);
extern ofstream logfile;
template<typename MSG>										//base version
void log(MSG msg){
	logfile << msg << "\n";
}
template<typename MSG, typename VAL1, typename... VALS>
void log(MSG msg, VAL1 v, VALS... vals){					//variadic version
	logfile << msg << "\t" << v;
	log("", vals...);
}
template<typename MSG>
void log(MSG msg, const Problem & p){						//log problem
	logfile << msg << ":" << "ID=" << p.id << ";N=" << p.N 
		<< ";startT=" << p.startT << ";LJ="<<p.jL->id << "\n[";
	FOR_E(i, p.N)logfile << p.jobs[p.edd[(i)]].id << ",";
	logfile << "]\n";
}

/////////////////////////
// Pack of counters
/////////////////////////
class CtrPack{
public:
	CtrPack() : ctrAllNodes(-2){}
	long long ctrAllNodes;  
	int ctrMovL;			
	int ctrMerL;			
	long long ctrSzMerL;	
	int ctrMovR;			
	int ctrMerR;			
	long long ctrSzMerR;	
	long long ctrJobsMem;
	long long ctrNodeSingleKid;
	long long ctrNodeMultiKid;
	long long ctrKidsNotSingle;
};

/////////////////////////
// Class Problem
/////////////////////////
class Problem{
private:
public:
	unsigned int id;
	int P;					// Sum of p
	int startT;
	PtrJob jL;		        // Longest job
	PtrJob jobs;
	PtrJobIndex edd;
	PtrPos rList;           // Reduced List: list of candidate positions to branch on 
#ifdef ENABLE_DB_SIGMA_CODE
	PtrJobIndex sigma;		// Prefixed job seq ids
	short N_sigma;
#endif
	short indN_EDD;			// Index of the longest job in EDD sequence
	short N;
	short N_rlist;
	bool _isRLComputed;		// Whether the reduced list already computed
	bool _isEddSorted;
	//char mmToldMe;

	static unsigned int idNext;
	static vector<Job> currIns; // Original input jobs at the root node. Sorted by id starting from id=1. 
	// Used for reading original due date etc. Set this once for each global instance
	static Mem pbMemo;
	static Mem sigmaMemo;
	static long long ctrAllNodes;		// total number of nodes generated
	static long long ctrMovL;			// number of nodes moved
	static long long ctrMerL;			// number of nodes merged
	static long long ctrSzMerL;			// total size of merged nodes
	static long long ctrMovR;			// number of nodes moved         (right)
	static long long ctrMerR;			// number of nodes merged        (right)
	static long long ctrSzMerR;			// total size of merged nodes    (right)
	static long long ctrJobsMemNotUsed;		// number of nodes left-merged with size > N/2, i.e. those should also be merged when right-m enabled
	static long long ctrNodeSingleKid;	// nb of nodes having only one kid
	static long long ctrNodeMultiKid;	// nb of nodes having only one kid
	static long long ctrKidsNotSingle;	// nb of kids with brothers

	static void resetStatic(){
		Problem::ctrAllNodes = 1;
		Problem::idNext = Problem::ctrMovL = Problem::ctrMerL = Problem::ctrSzMerL =
			Problem::ctrMovR = Problem::ctrMerR = Problem::ctrSzMerR = Problem::ctrJobsMemNotUsed =
			Problem::ctrKidsNotSingle = Problem::ctrNodeSingleKid = Problem::ctrNodeMultiKid = 0;
		Problem::pbMemo = Mem();	// Also reset hashtable inner attributs like load_factor
		Problem::sigmaMemo = Mem();
	}

public:
	Problem(short size, int startTime = 0);
	Problem(PtrJobConst beg,
		int n,
		int startTime = 0,
		bool reuseEL = false,
		bool eddSorted = false);
	Problem(PtrJobConst jobsSrc,
		PtrJobIndexConst indBegSrc,
		int n,
		int startTime = 0,
		bool reuseEL = false,
		bool eddSorted = false);
	Problem(Problem &&p) throw() :Problem(0) {
		*this = move(p);
	}
	~Problem()
	{
		FREE<Job>(&jobs, N);
		FREE<JobIndex>(&edd, N);
		FREE<Pos>(&rList, N);
#ifdef ENABLE_DB_SIGMA_CODE
		if (N_sigma > 0)FREE<JobIndex>(&sigma, N_sigma);
#endif
	}

private:
	// Copy constructor, manage table of pointers: edd spt...
	explicit Problem(Problem &){
		cerr << "copy constructor is disabled.\n";	FAIL_ON(true);
	}
	Problem & operator =(Problem &){
		cerr << "The = operator is disabled.\n";    FAIL_ON(true);
	}
	Problem & operator =(Problem &&){
		cerr << "Move assignment is disabled\n";    FAIL_ON(true);
	}

public:
	// This should be called when all jobs are added
	void init(int startTime = 0);
	
	// Compute the list of positions where the longest job can be branched on
	// Rules implemented: d version: 3,4
	PtrPos reducedList();

	// Compute T(n,i) in the paper: tt of the edd sequence with the longest job switched to pos i
	int* Problem::computeTni(PtrJobIndexConst edd2, int n, short indN, int startT);

	// This is the method to call
	PtrJob solve(int & ttOut, short paraK = 0);

	// Solve by brute force, used for very small ins. Default: 5
	static PtrJob solveBrute(PtrJob jobs, short n, int t0, int & ttOut);

	PtrJob branch(int & ttOut);

	// Double branching. List reducing is also inside.
	PtrJob branch2(int & ttOut);

	// Remove already used positions from reduced list, mainly used when creating kids
	void cleanRList()
	{
		auto end = remove(rList, rList + N_rlist, -1);
		N_rlist = end - rList;
	}

	// Do split. Return the opt sequence or NULL if split cannot be performed
	PtrJob split(short k, int & ttOut);

	void abandon(){ N = -1; };	// Current node will not lead to an optimal solution. It can be decided for example by sigma test on bd2.
	bool isAbandoned(){ return N == -1; }

	// Check sigma in the second db, return true if strictly dominated. False if not. Update db2 if necessary
	bool checkSigma(){
		//!TAG0216 return false;
		if (!Config::ENABLE_DB_SIGMA || N < Config::MEM_N_LB)return false;	// Skip small pbs, to avoid split
		bool res = false;
#ifdef ENABLE_DB_SIGMA_CODE
		auto tt_sump = TT(sigma, N_sigma);
		auto searchRes = Problem::pbMemo.findSol(N_sigma, 0, tt_sump.second, sigma);
		if (searchRes == NULL)
		{
			// Check in db2
			searchRes = Problem::sigmaMemo.findSol(N_sigma, 0, tt_sump.second, sigma);
			if (searchRes == NULL){// Add sigma
				Problem::sigmaMemo.addSol(N_sigma, 0, tt_sump.second, tt_sump.first, sigma, -1);
			}
			else{
#ifdef USE_MEM2
				int & tt = searchRes->t0_tt[0];
#else
				int & tt = searchRes->tt;
#endif
				if (tt < tt_sump.first)res = true;
				else if (tt > tt_sump.first){ // Update db2
					tt = tt_sump.first;
					mcopy(sigma, searchRes->sol, N_sigma);
					res = false;
				}
			}
		}
#ifdef USE_MEM2
		else if (searchRes->t0_tt[0] < tt_sump.first)
#else
		else if (searchRes->tt < tt_sump.first)
#endif
			res = true;
		if (res) {
			++ Stat::ctrCutDb2;
			Stat::szCutDb2 += N;
		}
#endif
		return res;
	}

	// Tried on branching jobs, no effect, possibly covered by posElim.
	bool testAdjacentOrder(int t0, const Job& j1, const Job& j2, int w1=1, int w2=1){
		int lhs = w1*j1.T(t0 + j1.p) + w2*j2.T(t0 + j1.p + j2.p);//w1*(j2.p - max(0, j1.d - t0 - j1.p));
		int rhs = w2*j2.T(t0 + j2.p) + w1*j1.T(t0 + j1.p + j2.p);//w2*(j1.p - max(0, j2.d - t0 - j2.p));
		if ( lhs>rhs ) 
			return false;
		//else if (lhs == rhs && j1.p < j2.p) 
			//return false;
		return true;
	}

	template<typename T>
	Job* cloneFromIds(const T srcids, int n)
	{
		if ((srcids) != NULL && n > 0)
		{
			Job* dest = ALLOC<Job>(n);
			FOR_E(i, n)dest[i] = currIns[srcids[i]-1];
			return dest;
		}
		return NULL;
	}

	// Return total tardiness and completion time    
	pair<int, int> TT0(PtrJobConst beg, int n, int C = 0)
	{
		int tt = 0;
		FOR_E(it, n){
			C += beg[it].p;
			tt += beg[it].T(C);
		}
		return make_pair(tt, C);
	}
	template<typename JOBIDS>
	static pair<int, int> TT(JOBIDS beg, int n, int C = 0)
	{
		int tt = 0;
		FOR_E(it, n){
			C += currIns[ID(beg[it]) - 1].p;
			tt += currIns[ID(beg[it]) - 1].T(C);
		}
		return make_pair(tt, C);
	}
	
	static CtrPack saveRestoreCtr(CtrPack toRestore = CtrPack()){
		if (toRestore.ctrAllNodes == -2){// save
			toRestore.ctrAllNodes = ctrAllNodes;
			toRestore.ctrMovL = ctrMovL;
			toRestore.ctrMerL = ctrMerL;
			toRestore.ctrSzMerL = ctrSzMerL;
			toRestore.ctrMovR = ctrMovR;
			toRestore.ctrMerR = ctrMerR;
			toRestore.ctrSzMerR = ctrSzMerR;
			toRestore.ctrJobsMem = ctrJobsMemNotUsed;
			toRestore.ctrNodeSingleKid = ctrNodeSingleKid;
			toRestore.ctrNodeMultiKid = ctrNodeMultiKid;
			toRestore.ctrKidsNotSingle = ctrKidsNotSingle;
		}
		else{//restore
			ctrAllNodes = toRestore.ctrAllNodes;
			ctrMovL = toRestore.ctrMovL;
			ctrMerL = toRestore.ctrMerL;
			ctrSzMerL = toRestore.ctrSzMerL;
			ctrMovR = toRestore.ctrMovR;
			ctrMerR = toRestore.ctrMerR;
			ctrSzMerR = toRestore.ctrSzMerR;
			ctrJobsMemNotUsed = toRestore.ctrJobsMem;
			ctrNodeSingleKid = toRestore.ctrNodeSingleKid;
			ctrNodeMultiKid = toRestore.ctrNodeMultiKid;
			ctrKidsNotSingle = toRestore.ctrKidsNotSingle;
		}
		return toRestore;
	}

	// --------------------------------------------------------------
	// Experimental functions
	// --------------------------------------------------------------
#pragma region TEST
	// Test whether early jobs are in edd
	bool testEarlyEdd(PtrJob sol, int n, int t0){
		FAIL_ON(n < 1 || sol == NULL);
		int lastD = sol[0].d;
		FOR_E(i, n){
			t0 += sol[i].p;
			if (t0 > sol[i].d){ lastD = sol[i].d; continue; } // Reaching the late job
			if (i != 0 && lastD > sol[i].d) return false;
			lastD = sol[i].d;
		}
		return true;
	}

	// Test whether optimal sequence respect Emmons, consider E=L=C
	bool testOptEmmons(PtrJob sol, int n, const vector<int> & Cs){
		FOR_E(i, n - 1){
			int maxed = Cs[i]>sol[i].d ? Cs[i] : sol[i].d;
			FOR_BE(j, i + 1, n){
				if (cmpJobP0(sol[i], sol[j]) && sol[i].d >= Cs[j] && sol[i].d >= sol[j].d){
					cerr << "First cond not verified i=" << i << "; j=" << j << endl;
					return false;
				}
				if (!cmpJobP0(sol[i], sol[j]) && (sol[j].d < maxed || sol[j].d + sol[j].p < Cs[i])){
					cerr << "Second cond not verified i=" << i << "; j=" << j << endl;
					return false;
				}
			}
		}
		return true;
	}

	// Compute the time window [t',t"] for t0 in which sol stays optimal
	pair<int, int> findOptWin(PtrJob sol, int n, int t0){
		vector<int> Cs(n, 0);
		Cs[0] = t0;
		FOR_E(i, n){
			if (i > 0)Cs[i] = Cs[i - 1];
			Cs[i] += sol[i].p;
		}
		if (!testOptEmmons(sol, n, Cs))
		{
			cerr << "Emmons test failed!\n";
			cin >> t0;
		}
		int deltaMinus = t0;// , deltaPlus = INT_MAX;
		FOR_E(i, n-1){
			if (deltaMinus == 0) break; //no need to continue;
			FOR_BE(j, i + 1, n){
				if (cmpJobP0(sol[i], sol[j])){
					// ! We depends on the modified Emmons' rules
					if (sol[i].d < sol[j].d)continue;
					if (sol[i].d < Cs[j] ) {//update deltaminus
						if (deltaMinus> Cs[j] - sol[i].d - 1){
							deltaMinus = Cs[j] - sol[i].d - 1;
							if (deltaMinus == 0) break; //no need to continue;
						}
					}
					else 
						cerr << "Unknown situation;";
				}
			}
		}
		
		//return make_pair(t0, t0 + deltaPlus);
		return make_pair(t0 - deltaMinus, t0);
	}

	bool testOptWin(PtrJob sol, int n, int t0){
		//if (n > 10)return true;
		auto win = findOptWin(sol, n, t0);
		if (win.first == win.second)return true;
		//cout << "\nPb size = " << n << "; Window = [" << win.first << "," <<win.second <<"]"<< endl;
		Problem p(sol, n, t0);
		int tt;
		for (int i = win.second - 1; i >= win.first; i--){
			if (i==t0) continue;
			p.startT = i;
			//p._isEiLiComputed = false;
			p._isRLComputed = false;
			//p._isSeiSorted = false;
			p._isEddSorted = false;
			FREE<Pos>(&p.rList, p.N);
			auto res = p.solve(tt);
			FOR_E(j, n)
				if (res[j] != sol[j]){
					int tt0 = TT(sol, n, i).first;
					tt = TT(res, n, i).first;
					if (tt0 > tt){
					cout << "\nPb size = " << n << "; Window size=" << win.second - win.first + 1 << endl;
					cout << "Not matching! id = " << p.id << endl;
					cout << "tt sol = " << tt0;
					cout << "\ntt res = " << tt;
						getchar();
						getchar();
					}
					FREE<Job>(&res, p.N);
					return false;
				}
			FREE<Job>(&res,p.N);
		}
		return true;
	}

	void testAllEarlyLate(PtrJob sol, int n, int t0, int tt){
		if (tt == 0){
			cout << id << ":" << N << ":All early\n";
			return;
		}
		FOR_E(i, n){
			t0 += sol[i].p;
			if (t0 <= sol[i].d)return;
		}
		cout << id << ":" << N << ":All late\n";
	}

	// Test if the initial instance can be simplified: some jobs are never tardy
	// ! Does not work for [0.2,0.6] since d is in [P(1-T-R/2), P(1-T+R/2)]=[0.3P, 0.5P]
	void kernelize(){
		int sump = P;
		int ctrEarly = 0;
		for (int i = N - 1; i >= 0; --i){
			if (EDD(i).d >= sump){
				sump -= EDD(i).p;
				++ ctrEarly;
			}
		}
		cout<<id<<" : " <<ctrEarly <<"/"<<N<< " jobs always early. \n";
	}

#pragma endregion
};

//////////////////////////////////////////////
//
// Inline functions have to be defined here
// Only inline basic functions, i.e. without recursive call
//////////////////////////////////////////////
INLINE PtrPos Problem::reducedList()
{
	if (_isRLComputed) return rList;
	if (Config::DISABLE_REDUCTION){
		rList = ALLOC<Pos>(N - indN_EDD);
		N_rlist = 0;
		for (int i = indN_EDD; i < N; ++i)rList[N_rlist++] = i;
		_isRLComputed = true;
		return rList;
	}
	bool rule1, rule2, rule3, rule4;
	rule1 = rule2 = false;
	rule3 = rule4 = true;	// To be consistent with Andrea's version
	N_rlist = 0;
	FAIL_ON(rList != NULL);
	if (N == 0) return rList;

	rList = ALLOC<Pos>(N-indN_EDD);

	int* Cn = ALLOC<int>(N - indN_EDD + 1);                                   // Cn[0] = Cn(s) in the paper
	int maxDP = 0;                                                            // Used for rule 3: max d_[i]+p_[i] for s+1<=i<=r-1
	PtrJobIndex edd2 = edd;
	int indN_EDD2 = indN_EDD;
	//Tni(r-indN_EDD2) is the tt of edd2 with the longest job moved to position r
	int* Tni = computeTni(edd2, N, indN_EDD2, startT);				  

	int minTT = INT_MAX;													  // Used for rule 4.3
	Cn[0] = startT;
	for (int i = 0; i <= indN_EDD; i++) Cn[0] += EDD(i).p;                    // Cn[0] = Cn(s) now
	for (int r = indN_EDD; r < N; r++){// Test each potential position

		if (r > indN_EDD){
			Cn[r - indN_EDD] = Cn[r - indN_EDD - 1] + EDD(r).p;               // Cn[r-indN_EDD] = Cn(r) 
		}
		if (r >= indN_EDD + 2)
			maxDP = max(maxDP, EDD(r - 1).d + EDD(r - 1).p);                   // For Rule 3.
		// Get the min T(n,i) for all i<=r-1
		if (r >= indN_EDD2 + 1){
			int tt_r_1 = Tni[r - 1 - indN_EDD2];
			if (minTT > tt_r_1)
				minTT = tt_r_1;
		}
		bool skip1 = false, skip2 = false;
		// Rule 1
		if (rule1)
			if (r<N - 1 && Cn[r - indN_EDD] > EDD(r + 1).d)
				skip1 = true;//continue;      // Skip this position
		// Rule 2 (may be integrated to rule 3)
		if (rule2)
			if (r > indN_EDD && Cn[r - indN_EDD] <= EDD(r).d + EDD(r).p)
				skip2 = true;//continue;

		// Rule 3
		if (rule3)
			if (r >= indN_EDD + 2 && Cn[r - indN_EDD] <= maxDP)
				continue;
		// Rule 4
		if (rule4){
			if (indN_EDD2 > r) continue;
			int t_n_r = Tni[r - indN_EDD2];
			if ((r >= indN_EDD2 + 1 && t_n_r >= minTT)
				|| (r<N - 1 && t_n_r > Tni[r + 1 - indN_EDD2]))
				continue;
			if (skip1 || skip2)
				cerr << "Not ok! " << id;

		}
		// Add r to the candidate list
		rList[N_rlist++] = r;
	}

	_isRLComputed = true;
	if (N_rlist == 0)
	{
		cout << "Empty reduce list! Suspended..." << endl;
		cout << id << " " << N << " " << startT << " " << *jL << endl;
		FOR_E(i, N)cout << EDD(i) << endl;
		getchar();
	}
	FREE<int>(&Tni, N-indN_EDD);
	FREE<int>(&Cn, N-indN_EDD+1);
	LOG("RedPos", strPoList(rList, N_rlist));
	return rList;
}

/// Compute T(n,i) from the paper: the tt of the edd sequence with the indN-th job switched to ind
INLINE int* Problem::computeTni(PtrJobIndexConst edd2, int n, short indN, int startt)
{
	int * res = ALLOC<int>(n - indN);
    int tt = 0, startTN = 0;
    for (int i=0; i < n; i++){
		if (i == indN) startTN = startt;
        startt += jobs[edd2[i]].p;
        tt += jobs[edd2[i]].T(startt);
    }
	res[0] = tt;
	Job & longest = jobs[edd2[indN]];
	// deduce res[1...N-indN]
    for (int i = 1; i <(n-indN); i++){
		// res[i-1]: jL is at indN+(i-1), the job after is indN+i. To get res[i]: we switch these 2 jobs
		Job & after = jobs[edd2[indN + i]];
		res[i] = res[i - 1]
			- (longest.T(startTN + longest.p) + after.T(startTN + longest.p + after.p))
			+ after.T(startTN + after.p) + longest.T(startTN + after.p + longest.p);
		startTN += after.p;
    }
    return res;
}

#endif
