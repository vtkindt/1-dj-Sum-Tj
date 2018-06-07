#ifndef H_H
#define H_H
/**
* This file defines several auxilary classes:
* Config: used to read parameters from config.ini
* Result: used to store results
* Stat: used to make statistics
* Alloc: used to manage memory allocation
*/

#include <cassert>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <direct.h>
#include <ctime>
#include <tuple>
#include <list>
#include <numeric>
#include "windows.h"
#include "psapi.h"	//for Windows' getProcess...
using namespace std;

///////////////////////
// Triggers
///////////////////////
//#define VLD_FORCE_ENABLE		// Comment this line to disable leak detection
//#include "vld.h"				// Comment this line to disable leak detection
#define BRUTE					// Solve small nodes by brute force. !! Disable this if we want early jobs in edd
#define BRUTE_MAX 5
#define MAX_N_JOBS 1300
#define MAX_P_SUM 130000
#define _1M 1000000
//#define VRIFY_RES				// Enable the code for verifying opt value

// MEM--------------------
#define UNIT_JOB_SEQ_LEN	64	// Used in Seq class: when requiring a job sequence, we allocate by this unit size (nb of jobs)
#define MEM_STRA_RATIO_LB	0	// We start by a small value to see what can we get. (0.3). This should not be important. (obselete)
#define MEM_STRA_NEED_FACTOR (10) // To store a new sol of 9 bytes, we clear out 9*FACTOR bytes!
// According to the test on 700.30 Lim=2G. Stra=5:
// - malloc/free		ram=2.6G
// - dlmalloc/dlfree	ram=2.1G (prefered)
// - USE_VIRTUAL_ALLOC	ram>3G,  and much slower(perhaps cause of 0-initialization). should not be enabled.
//							Note that VirtualAlloc/Free get/release mem from/to the system directly. dlmalloc is based on this by adding
//							a system keeping freed blocks for further use (best-fit)
#define MALLOC_FUNC dlmalloc    // can be dlmalloc
#define FREE_FUNC dlfree
//#define USE_VIRTUAL_ALLOC		// After testing on 700.30 with Lim=2G, this should not be used: much longer and much more ram, even though the ram decrease frequently
//#define ENALBE_ALLOC			// Use Alloc class for mem pool management: this does not concern mem.h in which native allocation by malloc is preferred
								// TODO test this on a big instance with swapping happen

#define USE_MEM0				// The database using single linked list
//#define MEM_FIELD_ID			// Add a field to SolvedPb. Enable this when log events
#define MEM_FIELD_HIT			// Add a field to SolvedPb. Enable this for stra 9

//#define ENABLE_DB_SIGMA_CODE	// Enable code concerned by db_sigma

//#define USE_MEM2				// The database with same sequences stored once for different t0's
//#define USE_MEM3				// The database without list
//#define USE_MEM4				// Mem2 without list
//#define MEM_LOG				// Log memory: for each size, avg chains, avg hits...
//#define MEM_LOG_2				// Created to test whether there exist many pbs with diff t0 share the same sol
//#define MEM_LOG_3				// Created to log sol that contains optimal subsol
//#define MEM_LOG_4				// For each N, log nbItem, nbItemUseful, nbBucket
//#define MEM_LOG_ML				// Log all sets and their hits for machine learning training

#define MAX_SIZE_UNIT_ALLOC 8  // The max sizeof(T) for T in Alloc<T>. Default to sizeof Job
//#define TEST_NB_KIDS
//#define LOG_ENABLE			// Enable logging to log.txt
//#define PRINT_SOL 
//#define STAT_GEN_MODE			// Enable this to generate stat from res files. Set start id to 201
//#define USE_JOBSET_BIS		// Enable another instance set

#ifdef LOG_ENABLE
#define LOG(...) log(__VA_ARGS__);
#else
#define LOG(...) 
#endif
#define FOR_E(i,n) for(int i=0; i<n; ++i)
#define FOR_BE(i, beg, end) for(int i=beg; i<end; ++i)
#define CLONE(arrptr, n)  ((arrptr==NULL)? NULL : mcopy(arrptr,ALLOC<remove_reference<decltype(*arrptr)>::type>(n),n) )
#define FAIL_ON(cond) {if(cond){ cerr<<"Failed at line "<< __LINE__ <<"; file "<< __FILE__<<", ("<<#cond<<") should not be true. Pause.\n"; getchar();}}
//#define FREE(mem)   if(mem!=NULL){free(mem);mem = NULL;}
//#define FAIL_ON(cond,msg) {if(cond){ cerr<<"Failed at line "<< __LINE__ <<"; file "<< __FILE__<<", ("<<#cond<<") should not be true. <<msg<<\n"; throw("");}}
#define MAX(a,b) ((a)>(b)?(a):(b))

extern "C" void* dlmalloc(size_t);
extern "C" void dlfree(void *p);
class Problem;
extern ofstream logfile;

template<typename T>
inline T* ALLOC(int n)
{
	if (n <= 0)return NULL;
#ifdef ENALBE_ALLOC
	return Alloc::alloc<T>(n);
#else
	return ALLOC0<T>(n);
#endif
}

template<typename T>
inline T* ALLOC0(int n)
{
	if (n <= 0)return NULL;
#ifdef USE_VIRTUAL_ALLOC
	T* res = (T*)VirtualAlloc(NULL, sizeof(T)*(n), MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
#else
	T* res = (T*)MALLOC_FUNC(sizeof(T)*(n));
#endif

	if (res == NULL){
		cerr << "ALLOC0: Fail on memory allocation ! sizeof=" << sizeof(T) << "; n=" << n << "; membytes = " << Problem::pbMemo.bytesMem << "; nbentry = " << Problem::pbMemo.nbentry << "; ram = " << get_ram_usage() << ". Deleting mem items...\n";
		Config::MALLOC_FAILED = true;
		Mem& memo = Problem::pbMemo;
		long long lim = memo.DB_LIM;
		memo.DB_LIM = memo.bytesMem - 10*_1M;		// At a step of 10M
		memo.removeOldDownToLim(false);
#ifdef USE_VIRTUAL_ALLOC
		while ((res = (T*)VirtualAlloc(NULL, sizeof(T)*(n), MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE)) == NULL){
#else
		while ((res = (T*)MALLOC_FUNC(sizeof(T)*(n))) == NULL){
#endif
			memo.DB_LIM -= _1M;		// At a step of 1M
			if (memo.DB_LIM < _1M){
				cerr << "No more space even we remove most sols.Pause." << endl;
				cerr << "ALLOC0: Fail on memory allocation ! " << " Memory size = " << memo.bytesMem << " bytes; Nb of sequences = " << memo.nbentry << "\n";
				logfile << "ALLOC0: Fail on memory allocation ! " << " Memory size = " << memo.bytesMem << " bytes; Nb of sequences = " << memo.nbentry << "\n";
				getchar();
			}
			memo.removeOldDownToLim(false);
		}
		// Remove another 10M for not supervised malloc
		memo.DB_LIM -= _1M;		// At a step of 1M
		if (memo.DB_LIM > _1M) memo.removeOldDownToLim(false);
		else cerr << "ALLOC0: Bytemem now < 1M !";
		memo.DB_LIM = lim;
	}
	return res;
}

template <typename T>
inline void FREE(T ** ptr, int n=-1){
	if (ptr==NULL || *ptr == NULL) return;
	if(n==-1){
		FREE_FUNC(*ptr);*ptr=NULL;return;	
	}
#ifdef ENALBE_ALLOC
	Alloc::release<T>(*ptr, n);
#else
	FREE0<T>(ptr);
#endif
}

template <typename T>
inline void FREE0(T ** ptr){
	if (ptr == NULL || *ptr == NULL) return;
#ifdef USE_VIRTUAL_ALLOC
	VirtualFree(*ptr, 0, MEM_RELEASE);
#else	
	FREE_FUNC(*ptr);
#endif
	*ptr = NULL;
}

template<typename T>
inline T* mcopy(const T* src, T* dest, int n) 
{
	if (src != NULL && n > 0)						// Without this test, behavior undefined 
		memcpy(dest, src, n*sizeof(T));
	return dest;
}

struct Config
{
    static int SIZE_START;		// =100 Program starts by solving instance of this size. 
	static int SIZE_STEP;
	static int ID_START;
	static bool SOLVE_ONE;		// Whether only solve one instance

	// 1: original memo; 
	// 2: memcut left+right ; (no more supported)
	// 3:mem_ass: add sub sols (different strategies are tested); 
	// 22: memcut2 + add subsols
	// 4: exclude first/last node from adding to bd
	// 5: search bd before split in solve. This changes nothing if we only add branched sol into bd
	// 6: 4+5
	// 7: add sol after split and search bd before split (5)
	// 8: = 4+7
    static int ENABLE_MEM_PB;
	static int MEM_N_LB;			   // Min size of sol to add into memory
	static int MEM_N_UB;
	static double MEM_N_UB_RATIO;		   // Do not add to memo if n>N/ratio. Default = 1.
    static long long RAM_LIM_MO;	   // The max size of ram usage. In MegaBytes. The result will be approx
	// Strategies when the memory is full
	// 0. [default] Add no more new sol.
	// 1. Replace LastUnused Or Last. Return if the chain is empty
	// 2. Replace LastUnused Or Last, otherwise remove global old entries until got enough space.
	// 3. (Obsolete, since clearly worse than 2) Directly remove old entries until got enough space.
	// 4. Replace ou remove all sol > N/RATIO
	// 5. 4 with dynamic ratio (each time starting with highest ratio)
	// 6. Replace ou remove big sols
	// 7. Consider the current global list when limit reached, remove unused among old, then old.
	// 8. Like stra5 but each time remove sequences on a range of sizes (defined by ratio step)
	// 9. Remove all unused and at the same time decrease 'hit' of remain sequences!
	//!!! TO BE CLEAR : JUST USE 9 !
	static int MEM_LIM_STRATEGY; 
	static double MEM_STRA_RATIO;	//Used by stra4. Default 1.
	static double MEM_STRA_RATIO_STEP;

	static bool DOUBLE_DEC;		// Enable double decomposition
    static bool ONLY_HARDEST;	// Only solve 0.2 0.6
	static bool USE_SLAVE;		// Whether use another program to solve
	static bool IS_SLAVE;		// Whether this is a slave program
	static bool SELF_SLAVE;
	static bool MEM_ALLOCATOR;	// Whether use my allocator for mem block allocation
	static string CONFIG_FILE;
	static string PROG_NAME;
	static string EXE;			// Program to call when TT2001 is enabled
	static string NOTIF;
	// Split is now always disabled for N>=MEM_N_LB. The flag only concerns small pbs
    static bool DISABLE_SPLIT;
    static int CURR_SIZE_INS;
	static bool MALLOC_FAILED;
	static bool ENABLE_DB_SIGMA;

    static bool readConfig(string fname)
    {
        bool isdefault = true;
        ifstream ifs (fname);
        if (!ifs)
        {
            cout << "\n Config file not found." << endl;
            return false;
        } else
        {
            char* buffer = _getcwd(NULL, 0);
			cout << "Config file: " << buffer << "\\" << fname << endl;
            free(buffer);
        }
        // Parse file
        string line;
        while (!ifs.eof()){
            getline(ifs, line);
            if (line.empty())continue;
            while (isspace(*line.rbegin()) && line.end() != line.begin())
                line.erase(line.end()-1);
            int pequal = int(line.find("="));
            if (pequal == string::npos || *line.rbegin() == '=') continue;
            string name = line.substr(0, pequal);
			string valstr = line.substr(pequal + 1);
			if (name == "NOTIF"){
				NOTIF = valstr;
				isdefault = false;
				continue;
			}
			else if (name == "EXE"){
				EXE = valstr;
				isdefault = false;
				continue;
			}
			else if (name == "MEM_STRA_RATIO"){
				stringstream(valstr) >> MEM_STRA_RATIO; continue;
			}
			else if (name == "MEM_N_UB_RATIO"){
				stringstream(valstr) >> MEM_N_UB_RATIO; continue;
			}
			else if (name == "MEM_STRA_RATIO_STEP"){
				stringstream(valstr) >> MEM_STRA_RATIO_STEP;  continue;
			}
			int value;
            stringstream(valstr)>>value;
            isdefault = false;
            if (name == "SIZE_START")SIZE_START = value;
			else if (name == "ID_START")ID_START = value;
			else if (name == "SIZE_STEP")SIZE_STEP = value;
            else if (name == "K")K = value;
			else if (name == "K_MIN")K_MIN = value;
			else if (name == "KID_MIN")KID_MIN = value;
            else if (name == "BM_MODE")BM_MODE = value;
			else if (name == "RAM_LIM_MO")RAM_LIM_MO = (value);
			else if (name == "MEM_LIM_STRATEGY")MEM_LIM_STRATEGY = (value);
			//else if (name == "MEM_STRA_RATIO")MEM_STRA_RATIO = (value);
			else if (name == "MEM_N_LB")MEM_N_LB = (value);
			else if (name == "MEM_N_UB")MEM_N_UB = (value);
			//else if (name == "MEM_N_UB_RATIO")MEM_N_UB_RATIO = (value);
            else if (name == "SOLVE_ONE")SOLVE_ONE = (value != 0);
			else if (name == "DOUBLE_DEC")DOUBLE_DEC = (value != 0);
            else if (name == "ENABLE_MEM_PB")ENABLE_MEM_PB = value;
            else if (name == "DEDUCE_EL")DEDUCE_EL = (value != 0);
			else if (name == "ONLY_HARDEST")ONLY_HARDEST = (value != 0);
			else if (name == "DISABLE_REDUCTION")DISABLE_REDUCTION = (value != 0);
			else if (name == "MEM_ALLOCATOR")MEM_ALLOCATOR = (value != 0);
			else if (name == "IS_SLAVE")IS_SLAVE = (value != 0);
			else if (name == "ENABLE_DB_SIGMA"){ 
#ifndef ENABLE_DB_SIGMA_CODE
				cerr << "DB_SIGMA should be enabled in h.h. Pause" << endl; getchar(); 
#endif
				ENABLE_DB_SIGMA = (value != 0);
			}
			else if ((name == "SIZE_MIN_CUT")|| (name == "MEMCUT_FACTOR") || name == "VERIFY_PREC" || name == "DPREM_VERSION" || (name == "RULE4_PREM") || (name == "DISABLE_SPLIT")){
				cerr << name <<" is no more supported. Pause\n"; getchar();
			}
        }
        if (K < 2)
        {
            BM_MODE=PREFER_RM=PREFER_LM = ENABLE_LM = ENABLE_RM = 0;
        }
        else
        {
            switch (BM_MODE)
            {
            case 1: ENABLE_LM = true; break;
            case 2: ENABLE_RM = true; break;
            case 3: ENABLE_LM = ENABLE_RM = true; break;
            case 4: PREFER_LM = ENABLE_LM = ENABLE_RM = true; break;
            case 5: PREFER_RM = ENABLE_LM = ENABLE_RM = true; break;
            default: ENABLE_LM = true;
            }
        }
		if (DPREM_VERSION) RULE4_PREM = false;
		if (MEMCUT_FACTOR < 2)MEMCUT_FACTOR = 2;
		if (ENABLE_MEM_PB == 2){ cerr << "MEM_CUT is no more supported.\n"; getchar(); }
		FAIL_ON(MEM_LIM_STRATEGY==3);
		if (!EXE.empty()){
			USE_SLAVE = true;
			if (EXE == "self.exe")SELF_SLAVE = true;
		}
        return isdefault;
    }
    static void printConfig()
    {
        cout << "SIZE_START      =" << SIZE_START <<
			"\nSIZE_STEP        =" << SIZE_STEP <<
			"\nID_START        =" << ID_START <<
            //"\nK               =" << K <<
			//"\nK_MIN           =" << K_MIN <<
			//"\nKID_MIN         =" << KID_MIN <<
            //"\nENABLE_LM       =" << ENABLE_LM <<
            //"\nENABLE_RM       =" << ENABLE_RM <<
            //"\nPREFER_LM       =" << PREFER_LM <<
            //"\nPREFER_RM       =" << PREFER_RM <<
            "\nSOLVE_ONE       =" << SOLVE_ONE <<
            //"\nENABLE_MOVE     =" << ENABLE_MOVE <<
            //"\nSIZE_MIN_CUT  =" << SIZE_MIN_CUT <<
			"\nMEM_LIM_STRATEGY=" << MEM_LIM_STRATEGY <<
			"\nMEM_STRA_RATIO=" << MEM_STRA_RATIO <<
			"\nMEM_STRA_RATIO_STEP=" << MEM_STRA_RATIO_STEP <<
			"\nRAM_LIM_MO=" << RAM_LIM_MO <<
			"\nMEM_N_LB=" << MEM_N_LB <<
			"\nMEM_N_UB=" << MEM_N_UB <<
			"\nMEM_N_UB_RATIO=" << MEM_N_UB_RATIO <<
			"\nENABLE_MEM_PB=" << ENABLE_MEM_PB <<
			"\nONLY_HARDEST=" << ONLY_HARDEST<<
			"\nDOUBLE_DEC=" << DOUBLE_DEC <<
			"\nMEM_ALLOCATOR=" << MEM_ALLOCATOR <<
			//"\nDISABLE_SPLIT=" << DISABLE_SPLIT << 
			"\nDISABLE_REDUCTION=" << DISABLE_REDUCTION << 
			"\nDB_SIGMA=" << ENABLE_DB_SIGMA <<
			"\nNOTIF=" << NOTIF <<
			endl; 
		if (!EXE.empty()){ cout << "\nUsing EXE =" << EXE << endl; }
    }
	static void notify(string msg=""){
		if (msg.empty() && !Config::NOTIF.empty())
			msg = Config::NOTIF;
		if (!msg.empty())
		{
			//cerr << "Notification disabled by lei..." << endl;
			system(("echo Subject: " + msg + " > mail.txt").c_str());
			system("curl.exe smtp://smtp.gmail.com:587 -v --mail-from \"yourmail@gmail.com\" --mail-rcpt \"yourmail@gmail.com\" --ssl -u yourmail@gmail.com:passwork -T mail.txt -k --anyauth");
			system("del mail.txt");
		} 

	}
	
	// -------------------------------------
	// Not used
	// -------------------------------------
	static int K;				// =0 (i.e. no merge)
	static int K_MIN;
	static int KID_MIN;
	static int BM_MODE;
	static bool ENABLE_LM;
	static bool ENABLE_RM;
	static bool PREFER_LM;
	static bool PREFER_RM;
	static int SIZE_MIN_CUT;		// apply mem cut when the larger pb is greater than this
	static int MEMCUT_FACTOR;	// apply mem cut when the smaller pb is smaller than N/factor
	static bool VERIFY_PREC;	// When computing precedence test to avoid duplication
	static bool DPREM_VERSION;	
	static bool DEDUCE_EL;
	static bool RULE4_PREM;
	static bool DISABLE_REDUCTION;
	static map<short, vector<tuple<int,int,int, int>>> TmpMap;
};


struct Result
{
	Result() :id(-1), tt(-1),t_cpu(-1), t_wall(-1), nbmovL(-1), nbmerL(-1), szmerL(-1), nbmovR(-1), 
		nbmerR(-1), szmerR(-1), nbnode(-1), hit(-1), miss(-1), nbentryMem(-1), nbJobsMem(-1), nbBytesMem(-1)
		, hit2(-1), miss2(-1), nbentryMem2(-1), nbBytesMem2(-1), ram(-1), rmFail(-1), minRmRatio(-1){}
	int id;
	int tt;
	double t_cpu, t_wall;

    long long nbmovL, nbmerL, szmerL;// , nbmerLHalf;
    long long nbmovR, nbmerR, szmerR;
    long long nbnode;
	long long hit, miss, nbentryMem, nbJobsMem, nbBytesMem;
	long long hit2, miss2, nbentryMem2, nbBytesMem2;
	long long ram, rmFail;
	double minRmRatio;
	string info;
};

struct Stat
{
	Stat() :R(-1), T(-1),tmin(-1), tavg(-1), tmax(-1),nbmovAvg(-1), nbmerAvg(-1), 
		sizemerAvg(-1), nbnodeAvg(-1), hitAvg(-1), missAvg(-1), nbentryMemAvg(-1), nbJobsMemAvg(-1), nbBytesAvg(-1){}
    float R, T;
    float tmin, tavg, tmax;
    long long nbmovAvg, nbmerAvg, sizemerAvg;
    long long nbnodeAvg;
	long long hitAvg, missAvg, nbentryMemAvg;
	long long nbJobsMemAvg, nbBytesAvg;

	static vector<pair<int,int>> rangeCtr1;
	static vector<pair<int, int>> rangeCtr2;
	static long long ctrMemReplace; // Nb of replacement happened
	static long long ctrMemRemove;	// Nb of mem cleaning happened
	static long long ctrMemRemoveFail;	// Nb of mem cleaning happened
	static long long szMemRemove;	// Total bytes removed
	static long long ctrCutDb2;
	static long long szCutDb2;
	static double minRmRatio;
	static void resetStatic(){
		szCutDb2 = ctrCutDb2 = szMemRemove = ctrMemRemove = ctrMemRemoveFail = ctrMemReplace = 0;
		minRmRatio = 1;
	}

};

// ------------------------------------
// Memory pool for efficient allocation (avoid memory fragmentation) 
// ------------------------------------
class Job;
class Alloc{
public:
	// One pool contains a list of free blocks of the same size
	struct Pool 
	{
	private:
		Pool& operator=(Pool&p) = default;
	public:
		explicit Pool(Pool&p){ cerr << "Copy constructor of Pool should not be called.\n"; };
		list<void*> blockList;
		Pool(){}
		int clean(){ // Release unused blocks. Called when memory not sufficient.
			if (blockList.empty())return 0;
			int nbRemove = int(blockList.size());
			for (auto & block : blockList)FREE_FUNC(block);
			blockList.clear();
			return nbRemove;
		}
		void * get(int n){
			if (!blockList.empty()){
				void* res = blockList.front();
				blockList.erase(blockList.begin());
				return res;
			}
			else{
				return (void*)ALLOC0<char>(n);
			}
		} 
		void put(void * ptr){
			blockList.push_back(ptr);
		}
		~Pool(){
			for (auto & block : blockList){
				FREE_FUNC(block);
			}
		}
	};

	static map<int, Pool> poolsBySize;		//poolsBySize[i] = a pool of blocks of size 2*i
	static long long memSize;				// In bytes
	static int maxIUsed;					//max effective index in poolsBySize
	static void init(){
		maxIUsed = memSize = 0;
		poolsBySize.clear();
	}
	template<typename T>
	static T* alloc(int n){
		n *= sizeof(T);
		auto & pool = poolsBySize[n];
		auto oldSize = pool.blockList.size();
		void * res = pool.get(n);
		if (pool.blockList.size()>oldSize) memSize += n;
		return (T*)res;
	}
	template<typename T>
	static void release(T * ptr, int n){
		n *= sizeof(T);
		poolsBySize[n].put(ptr);
	}
	static void clear(){
		//long long nbBytesRemoved = 0;
		//for (auto &it : poolsBySize)it.second.clean();
		poolsBySize.clear();
	}
	static void logAlloc(){
		for(auto &it : poolsBySize)logfile << it.first<< ", " << it.second.blockList.size()<<endl;
	}
};

// An event corresponds to one line in the output file
class Event{
public:
	enum Action{
		C_SPLIT, C_BRANCH, SOVLED, CUT_LB, CUT_MEMO, USED
	};
	Action act;
	unsigned int id;
	unsigned int pid;
	short brPos;
	char lmr;
	char solvedBy; //1 memo, 2 branching

	Event(Action actt, unsigned int idd, unsigned int pidd, short brPoss, char isLeftt, char solvedByy) : 
		act(actt), id(idd), pid(pidd), brPos(brPoss), lmr(isLeftt), solvedBy(solvedByy){}
	// Called for CREATE
	Event(Action actt, unsigned int idd, unsigned int pidd, short brPoss, char lmrr) :Event(actt, idd, pidd, brPoss, lmrr, 0){}
	// Called for CUT_?
	Event(Action actt, unsigned int idd,  unsigned int pid) :Event(actt, idd, pid, -1, false, -1){}

	static void logEvent(const Event & e){
		static string types[6] = {"C1","C2", "S","K1","K2","U"};
		logfile << types[e.act] << "," << e.id << "," << e.pid << ",";
		if (e.act == Action::CUT_MEMO)logfile << "Cut by " << e.pid;
		logfile << ",";
		if (e.act == Action::C_BRANCH)
			logfile << e.brPos << char(e.lmr); //link info
		else if (e.act == Action::C_SPLIT)
			logfile << "s";
		logfile << endl;
	}
};
#define LOGEVENT(...) 
//#define LOGEVENT(...) Event::logEvent(Event(__VA_ARGS__))

double get_wall_time();
double get_cpu_time();
unsigned long long get_cpu_cycle();
long long get_ram_usage();
#endif