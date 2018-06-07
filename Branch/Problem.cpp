#include "Problem.h"
#include "Job.h"
#ifdef INLINE
#undef INLINE
#endif
#define INLINE					//change this line to enable inline. Not so useful for recursive functions
#pragma warning(disable:4127)	// For FAIL_ON(true): l'expression conditionnelle est une constante

vector<Job> Problem::currIns;
vector<Job>* SolvedPb::ptrCurrIns=NULL;
Mem Problem::pbMemo, Problem::sigmaMemo;
unsigned int Problem::idNext=0;
long long Problem::ctrAllNodes=0;
long long Problem::ctrMovL = 0;
long long Problem::ctrMerL = 0;
long long  Problem::ctrSzMerL=0;
long long Problem::ctrMovR = 0;
long long Problem::ctrMerR = 0;
long long Problem::ctrSzMerR=0;
long long  Problem::ctrJobsMemNotUsed=0;
long long Problem::ctrNodeSingleKid = 0;
long long Problem::ctrNodeMultiKid = 0;
long long Problem::ctrKidsNotSingle = 0;

// This constructor allows to add jobs later. Then init() should be called manually
Problem::Problem(short size, int startTime) :
_isRLComputed(false),
_isEddSorted(false),
//mmToldMe(0),
jL(NULL),
edd(NULL),
rList(NULL),
#ifdef ENABLE_DB_SIGMA_CODE
	sigma(NULL),
	N_sigma(0),
#endif
N(size),
N_rlist(0),
P(0)//, tt(0)
{
    jobs = ALLOC<Job>(N);
    startT = startTime;
}

Problem::Problem(PtrJobConst beg,
    int n,
    int startTime, 
    bool reuseel, 
    bool eddSorted) : Problem(n, startTime){
    mcopy( beg,jobs, N);
    _isEddSorted = eddSorted;
    init(startTime);
    ++ctrAllNodes;
}

Problem::Problem(PtrJobConst jobsSrc, 
    PtrJobIndexConst beg,
    int n,
    int startTime,
    bool reuseEL,
    bool eddSorted) : Problem(n, startTime){
	//! ReuseEL no more supported since split is removed
    FOR_E(i, n)
    {
        jobs[i] = jobsSrc[beg[i]];
    }
    _isEddSorted = eddSorted;
    init(startTime);
    ++ctrAllNodes;
}

// This is called when all jobs are added
void Problem::init(int startTime ){
    if (startTime > 0) startT = startTime;
    if (N <= 0) return;
	id = idNext++;

    P = 0;
    edd = ALLOC<JobIndex>(N);

    jL = &jobs[0];
	int indLj = 0;
    for (int i = 0; i < N; i++){
        edd[i]=i;		
        P += jobs[i].p;
		if (cmpJobP0(*jL, jobs[i])){
			jL = &jobs[i];
			indLj = i;
		}
    }
    if (!_isEddSorted)
        msort(jobs, edd, N, cmpJobDIndex);
	
	// init id edd
	// FOR_E(i, N)EDD(i).idEdd = i;

    // Find the longest job in edd
    auto ind = find(edd, edd+N, indLj);
    assert(ind != edd+N);
    indN_EDD = short(ind - edd);
	LOG("New Node", *this);
}


INLINE PtrJob Problem::solve(int & ttOut, short paraK)
{
	LOG("Solve Beg", id);
	if (paraK < 2) paraK = 0;
	if (N == 0){ FAIL_ON(N==0); ttOut = 0; return ALLOC<Job>(0); }
	if (N == 1 )  {										// BUG 0721: a pb of size 1 can have a 0-size lkid
		ttOut = jobs[0].T(startT + jobs[0].p);
		return CLONE(jobs, 1);
	}
#ifdef BRUTE
	if (N <= BRUTE_MAX) {
		auto res = solveBrute(jobs, N, startT, ttOut);
		return res;
	}
#endif

	// Search first in bd
	if (Config::ENABLE_MEM_PB ){
		auto searchRes = Problem::pbMemo.findSol(N, startT, P, jobs);
		if (searchRes != NULL)
		{
			LOG("solve()", id, N, "Found in memory!");
			//cout << id << "," << N << endl;
#if defined(USE_MEM2) || defined(USE_MEM4)
			ttOut = searchRes->t0_tt[startT];
#else
			ttOut = searchRes->tt;
#endif
			//!TAG0216  if (searchRes->id != id)Stat::ctrMemReplace++; 
			LOGEVENT(Event::Action::CUT_MEMO, id, searchRes->id);
			LOGEVENT(Event::Action::USED, searchRes->id, id);
			return cloneFromIds(searchRes->sol, N);
		}
	}
    auto res = Config::DOUBLE_DEC ? branch2(ttOut) : branch(ttOut);
	LOG("Solve End", id, "Solved by branch(", paraK);
	return res;
}

INLINE PtrJob Problem::branch(int & bestTT)
{
	if (N <= 0) return NULL;
	LOG("Branch() Beg", *this);

	PtrPos posList = reducedList();
	if (N_rlist <= 0){
		// Current pb has no optimal solution. (This can happen when the Ei values comes from parent pbs)
		// Or when the only pos in rlist is consumed by makeLKids
		LOG("branch():N_rList<=0", N_rlist, " Return null");
		cerr << "Rlist empty. id,N,startT=" << id << "," << N << "," << startT << endl;
		getchar();
		return NULL;
	}
#ifdef TEST_NB_KIDS
	else if (N_rlist == 1)
	{
		if (rList[0] > 0 && rList[0] < N - 1)
		{
			ctrKidsNotSingle += 2;
			++ctrNodeMultiKid;
		}
		else	++ctrNodeSingleKid;
	}
	else {
		++ctrNodeMultiKid;
		FOR_E(i, N_rlist){
			if (rList[i] > 0 && rList[i] < N - 1)
			{
				ctrKidsNotSingle += 2;
			}
			else ++ctrKidsNotSingle;
		}
	}
#endif
	int tmpTT = 0;
	PtrJob sol = ALLOC<Job>(N);
	PtrJobIndex edd2 = ALLOC<JobIndex>(posList[N_rlist - 1]);// edd2 contains left jobs in edd order
	for (int iedd = 0, iedd2 = 0; iedd2 < posList[N_rlist - 1]; ++iedd, ++iedd2)
		if (iedd == indN_EDD) --iedd2; else edd2[iedd2] = edd[iedd];
	bestTT = INT_MAX;

	FOR_E(i, N_rlist){
		short p = posList[i];
		LOG("Branch()", id, jL->id, "->", p);
		PtrJob solTmpL = NULL, solTmpR = NULL;
		int tt = 0;
		int startTRight = startT;
		// Avoid p<= in case both flags be true
		if (p > 0){
			// Solve left
			Problem pleft(jobs, edd2, p, startT, Config::DEDUCE_EL, true);
			LOGEVENT(Event::Action::C_SPLIT, pleft.id, id, p, 'L');
#ifdef ENABLE_DB_SIGMA_CODE
			if(Config::ENABLE_DB_SIGMA && N_sigma > 0){
				pleft.sigma = CLONE(sigma, N_sigma);
				pleft.N_sigma = N_sigma;
			}
#endif
			// The sigma must be valid otherwise the father node should already be cut
			
			solTmpL = pleft.solve(tt);
			if (solTmpL == NULL){// Current node will not be optimal: bad sigma reached in one of its subnodes
				FAIL_ON (!pleft.isAbandoned());
				continue;
			}
			startTRight = pleft.P + startT;
		}
		startTRight += jL->p;
		tt += jL->T(startTRight);
		if (p < N - 1){							// Solve right
			Problem pright(jobs, edd + p + 1, N - p - 1, startTRight, Config::DEDUCE_EL, true);
			LOGEVENT(Event::Action::C_SPLIT, pright.id, id, p, 'R');
#ifdef ENABLE_DB_SIGMA_CODE
			if (Config::ENABLE_DB_SIGMA){
				// Set sigma
				pright.N_sigma = N_sigma + p + 1;
				pright.sigma = ALLOC<short>(pright.N_sigma);
				mcopy(sigma, pright.sigma, N_sigma);
				FOR_E(i, p)pright.sigma[N_sigma + i] = solTmpL[i].id;
				pright.sigma[pright.N_sigma - 1] = jL->id;
				//check sigma dominance
				if (pright.checkSigma())continue;
			}
#endif
			solTmpR = pright.solve(tmpTT);
			if (solTmpR == NULL){
				FAIL_ON( !pright.isAbandoned());
				continue;
			}
			tt += tmpTT;
			//if (Config::ENABLE_MEM_PB == 3){
			//	int n = pright.N;
			//	int sumP = pright.P;
			//	int t0 = startTRight;
			//	auto solmem = solTmpR;
			//	while (n > 10){
			//		auto searchRes = Problem::pbMemo.findSol(n, t0, sumP, solmem);
			//		if (searchRes.first == NULL)
			//		{// Add
			//			pbMemo.addSol(n, t0, sumP, tmpTT, solmem);
			//		}
			//		--n;
			//		// Mind the order
			//		sumP -= solmem->p;
			//		t0 += solmem->p;
			//		tmpTT -= solmem->T(t0);
			//		solmem += 1;
			//	}
			//}
		}
		LOG("Branch()", id, jL->id, "->", p, tt);
		// Keep best
		if (bestTT > tt){
			bestTT = tt;
			mcopy(solTmpL, sol, p);
			sol[p] = *jL;
			mcopy(solTmpR, sol + p + 1, N - p - 1);
		}
		FREE<Job>(&solTmpL, p);	// TODO: avoid allocate every time. Use a global double buffer
		FREE<Job>(&solTmpR, N-p-1);
		if (bestTT == 0){ //!AG
			//if(N_rlist - i - 1>0) cout << "Branch Skipping " << N_rlist - i - 1 << endl;
			break;	// The case where sol has TT=0, no need to continue.
		}
	}
	FREE<JobIndex>(&edd2, posList[N_rlist-1]);
	if (bestTT == INT_MAX){
		abandon();
		//cerr << "BestTT = INT_MAX"<<endl;
		return NULL;
	}

	if (Config::ENABLE_MEM_PB) {  //!AG && N_rlist > 1. Add condition of Andrea: do not addSol if the branching is unique
		if ((Config::ENABLE_MEM_PB != 4 && Config::ENABLE_MEM_PB != 6 && Config::ENABLE_MEM_PB != 8) //|| (mmToldMe != 1 && mmToldMe != 2)
			){
			pbMemo.addSol(N, startT, P, bestTT, sol,id);
			//pbMemo.addSol(N, _isLeft ? -startT : startT, _isRight ? -P : P, bestTT, sol);
			LOG("Branch()", id, "Add to memory!");
		}
	}
	LOG("Branch() End", id, strIdList(sol, N), bestTT);
	if (sol->id < 0){ cerr << id<<": Uninitiated sol." << endl; getchar(); }
	return sol;
}

// double decomposition
INLINE PtrJob Problem::branch2(int & bestTT)
{
	if (N <= 0) return NULL;
	// Single branching
	if (indN_EDD == 0) return branch(bestTT);
	LOG("Branch2() Beg", *this);
	FAIL_ON(_isRLComputed);
	if (Config::DPREM_VERSION){ cerr << "D' version is not supported by double branching!\n"; getchar(); }
	
	PtrJobIndex edd2 = ALLOC<JobIndex>(N);
	int nbJobs = 0;
	FOR_E(i, N)
		if (!cmpJobP0(EDD(0), EDD(i)))
			edd2[nbJobs++] = edd[i];

	vector<short> lawlerPos, dcPos;
	lawlerPos.reserve(N - indN_EDD);
	dcPos.reserve(nbJobs + 1);
	if (Config::DISABLE_REDUCTION){
		// Return all positions
		for (int i = indN_EDD; i < N; ++i)lawlerPos.push_back(i);
		FOR_E(i, nbJobs)dcPos.push_back(i);
	}
	else{
		bool rule3, rule4;
		rule3 = rule4 = true;	// To be consistent with Andrea's version
		// Compute positions for lawler
		int minTT = INT_MAX;													  // Used for rule 4.3, for d' version
		int maxDP = 0;                                                            // Used for rule 3: max d_[i]+p_[i] for s+1<=i<=r-1
		int* Cn = ALLOC<int>(N);												  // Cn2[0]= like Cn but the job switched is the smallest due date job
		int* Tni = computeTni(edd, N, indN_EDD, startT);					  // Tni(r-indN_EDD2) is the tt of edd2 with the longest job moved to position r
		Cn[0] = startT;
		for (int i = 0; i <= indN_EDD; i++) Cn[0] += EDD(i).p;                    // Cn[0] = Cn(s) now

		for (int r = indN_EDD; r < N; r++){// Test each potential position
			if (r > indN_EDD)
				Cn[r - indN_EDD] = Cn[r - indN_EDD - 1] + EDD(r).p;               // Cn[r-indN_EDD] = Cn(r) 
			if (r >= indN_EDD + 2)												  // !AG Changed for AG: should be +2, EDD(r-1)
				maxDP = max(maxDP, EDD(r - 1).d + EDD(r - 1).p);                  // For Rule 3.
			// Get the min T(n,i) for all i<=r-1
			if (r >= indN_EDD + 1){
				int tt_r_1 = Tni[r - 1 - indN_EDD];
				if (minTT > tt_r_1)
					minTT = tt_r_1;
			}
			// Rule 3
			if (rule3 && r >= indN_EDD + 2 && Cn[r - indN_EDD] <= maxDP) continue;
			// Rule 4
			if (rule4 && indN_EDD <= r){
				int t_n_r = Tni[r - indN_EDD];
				if ((r >= indN_EDD + 1 && t_n_r >= minTT)
					|| (r<N - 1 && t_n_r > Tni[r + 1 - indN_EDD]))
					continue;
			}
			// Add r to the candidate list
			lawlerPos.push_back(r);
		}

		// Compute positions for Decomposition II
		maxDP = 0;
		minTT = INT_MAX;
		FREE<int>(&Tni, N - indN_EDD);
		Tni = computeTni(edd2, nbJobs, 0, startT);
		Cn[0] = startT + EDD(0).p;					// C of edd with job k (smallest due date) switched to 0
		for (int r = 0; r < nbJobs; r++){			// Test each potential position
			if (r > 0)Cn[r] = Cn[r - 1] + EDD(r).p;
			if (r >= 2)															   //!AG Changed for AG: should be 2, EDD(r-1)
				maxDP = max(maxDP, EDD(r - 1).d + EDD(r - 1).p);                   // For Rule 3.
			// Get the min T(n,i) for all i<=r-1
			if (r >= 1){
				int tt_r_1 = Tni[r - 1];
				if (minTT > tt_r_1)
					minTT = tt_r_1;
			}
			// Rule 3
			if (rule3 && r >= 2 && Cn[r] <= maxDP) continue;
			// Rule 4
			if (rule4)
				if ((r >= 1 && Tni[r] >= minTT)
					|| (r<nbJobs - 1 && Tni[r] > Tni[r + 1]))
					continue;
			// Add r to the candidate list
			dcPos.push_back(r);
		}
		FREE<int>(&Cn,N);
		FREE<int>(&Tni, nbJobs);
	}

	if (lawlerPos.empty() || dcPos.empty()){
		LOG("branch2():double rList empty", lawlerPos.size(), dcPos.size(), " Return null");
		cerr << "Double dec, rlist empty, id,N,startT,posLawlerSize, posDecSize=" << id << "," << N << "," << startT << "," << lawlerPos.size() << "," << dcPos.size() << endl;
        return NULL;
	}

	bool isUnique = false;
	if (lawlerPos.size() == 1 && dcPos.size() == 1)
		isUnique = true;
	_isRLComputed = true;

	// For each compatible pos pair, decompose
	int tmpTT = 0;
    PtrJob sol = ALLOC<Job>(N);    
	PtrJobIndex eddMid = ALLOC<JobIndex>(lawlerPos.back());		// eddmid contains mid jobs in edd order
	bestTT = INT_MAX;
	int lastDcPosComp = 0;										// Max dc pos that is compatible with last lawer pos: all dcpos < that nb must be compatible with current lawer pos
	for (auto p : lawlerPos){
		for (auto p2 : dcPos){
			if (p2 > lastDcPosComp){
				// Test compatibility: maxd of dc's left must < min d of lawler's right
				if (p < N - 1 && !cmpJobD0(jobs[edd2[p2]], EDD(p + 1)))	continue;
				else lastDcPosComp = p2;
			}

			// Branching
			PtrJob solTmpL = NULL, solTmpM = NULL, solTmpR = NULL;
			int tt = 0;
			int startTMR = startT;

			if (p2 > 0){								
				// Solve left
				Problem pleft(jobs, edd2+1, p2, startT, Config::DEDUCE_EL, true);
				LOGEVENT(Event::Action::C_BRANCH, pleft.id, id, p2, 'L');
#ifdef ENABLE_DB_SIGMA_CODE
				if (Config::ENABLE_DB_SIGMA && N_sigma > 0){
					pleft.sigma = CLONE(sigma, N_sigma);
					pleft.N_sigma = N_sigma;
				}
#endif
				solTmpL = pleft.solve(tt);
				if (solTmpL == NULL){
					//cerr << __LINE__ << ": pleft.solve returns NULL. N,p2,p,pleft.N,pleft.startT=" << N << "," << p2 << "," << p << "," << pleft.N << "," << pleft.startT << endl;
					//cerr << strIdList(pleft.jobs, pleft.N)<<endl;
					FAIL_ON( !pleft.isAbandoned());
					continue;
				}
				startTMR += pleft.P;
			}
			startTMR += EDD(0).p;
			tt += EDD(0).T(startTMR);
			// Middle
			if (p2 + 1 < p){
				int nMid = 0;
				short nextEddLeft = 0;
				FOR_E(i, p+1){
					if (i == indN_EDD)continue;
					if (nextEddLeft<=p2 && edd[i] == edd2[nextEddLeft]){
						++nextEddLeft;
						continue;
					}
					eddMid[nMid++] = edd[i];
				}
				if (nMid != (p - p2 - 1))
					cerr << id<< " "<<nMid << " " << p << " " << p2 << endl;
				FAIL_ON(nMid != (p - p2 - 1));
				Problem pmid(jobs, eddMid, nMid, startTMR, Config::DEDUCE_EL, true);
				LOGEVENT(Event::Action::C_BRANCH, pmid.id, id, p2+1, 'M');
				//if (p2 == 0 && p == N - 1)pmid.mmToldMe = 23;
				//else if (p2 == 0)pmid.mmToldMe = 21;
				//else if (p == N - 1)pmid.mmToldMe = 22;
#ifdef ENABLE_DB_SIGMA_CODE
				if (Config::ENABLE_DB_SIGMA){
					// Set sigma
					pmid.N_sigma = N_sigma + p2 + 1;
					pmid.sigma = ALLOC<short>(pmid.N_sigma);
					mcopy(sigma, pmid.sigma, N_sigma);
					FOR_E(i, p2)pmid.sigma[N_sigma + i] = solTmpL[i].id;
					pmid.sigma[pmid.N_sigma - 1] = EDD(0).id;
					// Check sigma dominance
					if (pmid.checkSigma()) continue;
				}
#endif
				solTmpM = pmid.solve(tmpTT);
				if (solTmpM == NULL){
					//cerr << __LINE__ << ": pmid.solve returns NULL. N,p2,p,mid.N,mid.startT=" << N << "," << p2 << "," << p << "," << pmid.N << "," << pmid.startT << endl;
					//cerr << strIdList(pmid.jobs, pmid.N) << endl;
					FAIL_ON(!pmid.isAbandoned());
					continue;
				}
				tt += tmpTT;
				startTMR += pmid.P;
			}
			startTMR += jL->p;
			tt += jL->T(startTMR);
			// Right
			if (p < N - 1){							
				Problem pright(jobs, edd + p + 1, N - p - 1, startTMR, Config::DEDUCE_EL, true);
				LOGEVENT(Event::Action::C_BRANCH, pright.id, id, p, 'R');
#ifdef ENABLE_DB_SIGMA_CODE
				if (Config::ENABLE_DB_SIGMA){
					// Set sigma
					pright.N_sigma = N_sigma + p + 1;
					pright.sigma = ALLOC<short>(pright.N_sigma);
					mcopy(sigma, pright.sigma, N_sigma);
					FOR_E(i, p2)pright.sigma[N_sigma + i] = solTmpL[i].id;
					pright.sigma[N_sigma + p2] = EDD(0).id;
					FOR_E(i, p - p2 - 1)pright.sigma[N_sigma + p2 + i + 1] = solTmpM[i].id;
					pright.sigma[pright.N_sigma - 1] = jL->id;
					// Check sigma dominance
					if (pright.checkSigma()) continue;
				}
#endif
				solTmpR = pright.solve(tmpTT);
				if (solTmpR == NULL){
					//cerr << __LINE__ << ": pright.solve returns NULL. N,p2,p,r.N,r.startT=" << N << "," << p2 << "," << p << "," << pright.N << "," << pright.startT << endl;
					//cerr << strIdList(pright.jobs, pright.N) << endl;
					FAIL_ON(!pright.isAbandoned());
					continue;
				}
				tt += tmpTT;
			}
			// Keep best
			if (bestTT > tt){
				bestTT = tt;
				mcopy(solTmpL, sol, p2);
				sol[p2] = EDD(0);
				mcopy(solTmpM, sol+p2+1, p-p2-1);
				sol[p] = *jL;
				mcopy(solTmpR, sol + p + 1, N - p - 1);
			}
			FREE<Job>(&solTmpL, p2);	  // TODO: avoid allocate every time. Use a global double buffer
			FREE<Job>(&solTmpM, p - p2 - 1);
			FREE<Job>(&solTmpR, N - p - 1);
			if (bestTT == 0){ //!AG 
				break;		  // The case where sol has TT=0, no need to continue.
			}
		}// for each dc pos
    }// for each lawler pos
	
	FREE<JobIndex>(&edd2,N);
	FREE<JobIndex>(&eddMid,lawlerPos.back());
	if (bestTT == INT_MAX){
		abandon(); return NULL;
	}

	if (Config::ENABLE_MEM_PB )	// && !isUnique. Add condition of Andrea: do not addSol if the branching is unique
	{
		if ( (Config::ENABLE_MEM_PB != 4 && Config::ENABLE_MEM_PB!=6 && Config::ENABLE_MEM_PB!=8) //|| (mmToldMe != 1 && mmToldMe != 2)
			){
			pbMemo.addSol(N, startT, P, bestTT, sol,id);
			LOG("Branch()", id, "Add to memory!");
		}
	}
	LOG("Branch2() End",id, strIdList(sol,N),bestTT );
	return sol;
}

// Call this when n<=4
INLINE PtrJob Problem::solveBrute(PtrJob jobs, short n, int t0, int & ttOut){
	FAIL_ON(n > 5);
	if (n <= 0){
		ttOut = 0;
		return NULL;
	}
	else if (n == 1) {
		ttOut = jobs[0].T(t0 + jobs->p);
		return CLONE(jobs, 1);
	}

	short index[5] = {0,1,2,3,4};
	PtrJob sol = CLONE(jobs, n);
	ttOut = 0;
	int t02 = t0;
	FOR_E(it, n){
		t02 += jobs[it].p;
		ttOut += jobs[it].T(t02);
	}

	while (next_permutation(index, index+n)){
		int tt = 0;
		t02 = t0;
		FOR_E(it, n){
			t02 += jobs[index[it]].p;
			tt += jobs[index[it]].T(t02);
		}
		if (tt < ttOut){
			ttOut = tt;
			FOR_E(it, n){
				sol[it] = jobs[index[it]];
			}
		}
	}
	return sol;
}