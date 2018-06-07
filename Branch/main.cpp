#include <thread>
#include <iomanip>
#include <exception>
#include <cstdio>
#include "h.h"
#include "Problem.h"
//using namespace std::chrono;

// Here are default values of parameters. They may be overwritten by the config file
int Config::SIZE_START = 100;   
int Config::SIZE_STEP = 0;   
int  Config::ID_START = 1;      // 1-200
int Config::K = 0;
int Config::CURR_SIZE_INS = -1;
int Config::K_MIN = 4;          // When would br be called. Should > 2
int Config::KID_MIN = 4;        // makeKids do not generate nodes less than this
int Config::BM_MODE = 0;
bool Config::ENABLE_LM = false;
bool Config::ENABLE_RM = false;
bool Config::PREFER_LM = false;
bool Config::PREFER_RM = false;
bool Config::SOLVE_ONE = false;
bool Config::DOUBLE_DEC = false;
int Config::ENABLE_MEM_PB = 0;
long long Config::RAM_LIM_MO = 999999;
int Config::MEM_LIM_STRATEGY = 0;
double Config::MEM_STRA_RATIO = 1;
double Config::MEM_STRA_RATIO_STEP = 0.02;
int Config::MEM_N_LB = 10;
int Config::MEM_N_UB = MAX_N_JOBS - 1;
double Config::MEM_N_UB_RATIO = 1;
bool Config::ONLY_HARDEST = false;
bool Config::DISABLE_REDUCTION = false;
bool Config::USE_SLAVE = false;
bool Config::IS_SLAVE = false;
bool Config::SELF_SLAVE = false;
bool Config::MEM_ALLOCATOR = false;
bool Config::MALLOC_FAILED = false;
bool Config::ENABLE_DB_SIGMA = false;
string Config::CONFIG_FILE = "";
string Config::PROG_NAME = "";
string Config::EXE = "";
string Config::NOTIF = "";
map<short, vector<tuple<int, int, int, int>>> Config::TmpMap; //pid id t0 tt

// Obsolete
int Config::SIZE_MIN_CUT = 0;				// apply mem cut when the larger pb is greater than this
int Config::MEMCUT_FACTOR = 2;				// apply mem cut when the smaller pb is smaller than N/factor
bool Config::VERIFY_PREC = false;
bool Config::DPREM_VERSION = false;
bool Config::DEDUCE_EL = false;
bool Config::RULE4_PREM = false;			//Default use rule 4 (instead of 4')
bool Config::DISABLE_SPLIT = false;

// Static variables for class Stat
vector<pair<int, int>> Stat::rangeCtr1;
vector<pair<int, int>> Stat::rangeCtr2;
long long Stat::ctrMemReplace;
long long Stat::ctrMemRemove;
long long Stat::ctrMemRemoveFail;
long long Stat::szMemRemove;
long long Stat::ctrCutDb2;
long long Stat::szCutDb2;
double Stat::minRmRatio = 1;

// Static variables for class Alloc
map<int, Alloc::Pool> Alloc::poolsBySize;
long long Alloc::memSize = 0;
int Alloc::maxIUsed = 0;


// Var/functions declarations
string fileSuffix = ".txt"; //"_turk.txt";
string outFolder = "out\\";
string dataFolder = "data\\";
ofstream logfile;
int run();
void eatMemoryBy1G();
void getHardInfo();
void setMaxMemLimit(long long);
void resetStatic();
void writeResHeader(string filename);

// Main function
int main(int argc, char* argv[]){
    std::ios_base::sync_with_stdio(false);	// To be faster
	Config::PROG_NAME = string(argv[0]);
	Config::CONFIG_FILE = argc > 1 ? argv[1] : "config.ini";
	Config::readConfig(Config::CONFIG_FILE);
	if (argc == 4)
	{
		// If called in solve_one mode instead of tester mode
		Config::IS_SLAVE = true;
		Config::SELF_SLAVE = false;
		Config::USE_SLAVE = false;
		Config::SIZE_START = atoi(argv[2]);
		Config::ID_START = atoi(argv[3]);
		Config::SOLVE_ONE = 1;
		Config::EXE = "";
		Config::NOTIF = "";
		cout << "Slave: ";
	}else	Config::printConfig();
	FAIL_ON(sizeof(Job) != 8);	// See max_size_unit_block
	if (!Config::USE_SLAVE){
		//cout << "Set hard mem limit to " << Config::RAM_LIM_MO + 500 << endl;
		//setMaxMemLimit((Config::RAM_LIM_MO + 500) * _1M);
		//cin.get();
	}
    run();
	if (!Config::IS_SLAVE){
		//Config::notify();
		cout << "\nProgram finished. Press any key to quit.\n";
		//cin.get();
	}
    return 0;
}

void saveCurrentRes(string filename, const Result & res);
void makeStat(string filename, Result* result, int nbIns = 200, int nbInsPerRT = 10);
int run()
{
    stringstream ss;
    int step = 100;
	if (Config::SIZE_STEP >0) step = Config::SIZE_STEP;
	//system(("MD " + outFolder).c_str());
	for (int sIns = Config::SIZE_START; sIns <= 1500; sIns += step) {
		//if (sIns == 900)Config::notify("800 finished");
		//else if (sIns == 1000)Config::notify("900 finished");
		//else if (sIns == 1100)Config::notify("1100 finished");

#ifdef MEM_LOG_ML
		ss.str("ml\\200_hard.txt");
#endif
		int sizeIns = sIns, nbIns = 200, nbInsPerRT = 10;

		if (!Config::IS_SLAVE) {
			cout << "\n ===== Size: " << sizeIns << " Number: " << nbIns << " =====" << endl;
			ss.str("");
#ifdef MEM_LOG_ML
			ss << "ml\\";
#endif
			ss << outFolder << "res" << sizeIns << fileSuffix;
			writeResHeader(ss.str());
		}

		// Results to write
		Result result[200];

		// Read existing results from res. Used when not starting from the first instance
		if (
#ifndef STAT_GEN_MODE
			sIns == Config::SIZE_START &&
#endif
			!Config::SOLVE_ONE && Config::ID_START > 1)
		{
			ss.str("");
#ifdef MEM_LOG_ML
			ss << "ml\\";
#endif
			ss << outFolder << "res" << sIns << fileSuffix;
			ifstream fres(ss.str());
			char line[150];
			stringstream ssLine;
			for (int indId = 0; indId < Config::ID_START - 1; indId++)
			{
				fres.getline(line, 150);
				if (fres.fail()) {
					cerr << "Read res file failed at id = " << indId + 1 << ". Check the res file. \n";	// If failed check the res file
					getchar();
				}
				ssLine.str(line);
				ssLine.seekg(0);
				ssLine >>
					result[indId].id >>
					result[indId].tt >>
					result[indId].t_cpu >>
					result[indId].nbmovL >>
					result[indId].nbmerL >>
					result[indId].szmerL >>
					result[indId].nbmovR >>
					result[indId].nbmerR >>
					result[indId].szmerR >>
					result[indId].nbnode >>
					result[indId].hit >>
					result[indId].miss >>
					result[indId].nbentryMem >>
					result[indId].nbJobsMem;
				// Make sure the file is in good format: we changed the nbCol au fil de temps
				if (indId > 1) FAIL_ON(result[indId].id != result[indId - 1].id + 1);
			}
			fres.close();
		}

		vector<int> tts;                        // TT values which will be compared to old values in the file, for test purpose
		tts.reserve(nbIns);
		for (double R = 0.2; R <= 1.0; R += 0.2) {
			if ((R / 0.2) * 4 * nbInsPerRT < Config::ID_START) continue;
			for (double T = 0.2; T <= 0.8; T += 0.2) {
				if ((R / 0.2-1) * 4 * nbInsPerRT + (T/0.2)*nbInsPerRT < Config::ID_START) continue;
				if (Config::ONLY_HARDEST && (int(R*10) != 2 || int(T*10) != 6)) {
						continue;
				}
				for (int j = 0; j < nbInsPerRT; j++)// For each instance
				{
					int indId = int((R/0.2-1)*4*nbInsPerRT + (T/0.2-1)*nbInsPerRT + j);
					if (Config::SOLVE_ONE && sIns == Config::SIZE_START &&
						(indId + 1 < (Config::ID_START)))
						continue;
					// Here we read instances, solve them and measure time
					ss.str("");
					ss << dataFolder << sIns << "\\SDT_"<<sIns<<"_"<<setprecision(1)<<fixed<<R<<"_"<<T<<"_"<<(j+1)<< fileSuffix;
					ifstream file(ss.str());
					if (!file) {
						cout << ss.str() << " not found. Exiting..." << endl; return 0;
					}
					// Ignore comments
					while (file.peek() == '#')file.ignore(1000, '\n');
					short p;
					int d;
					Problem::currIns.clear();
					for (short k = 0; k < sizeIns; k++) {
						file >> p >> d;
						Problem::currIns.emplace_back(k + 1, p, d);
					}

					if (!file.good()) {
						cerr << "Problem when reading instance "<<ss.str() << endl;
						cerr << "Pause..." << endl;
						getchar();
					}

					
#ifdef STAT_GEN_MODE
					if (false)
#else
					if (Config::SOLVE_ONE || sIns > Config::SIZE_START || (sIns == Config::SIZE_START && indId + 1 >= Config::ID_START))
#endif
					{
						/////////////////////////
						// Solving...
						/////////////////////////
						
						// Call another program to solve each instance
						if (Config::USE_SLAVE) {
							cout << "==========\nSolving Instance " << ss.str() << "...";
							int pos = Config::PROG_NAME.find_last_of('\\');
							string currDir = (pos == Config::PROG_NAME.npos) ? "" : Config::PROG_NAME.substr(0, pos+1);
							if (Config::EXE == "self.exe") {
								ss.str("");
								ss << "copy " << Config::PROG_NAME << " " << currDir << "self.exe /Y";
								//cout << ss.str() << endl;
								system(ss.str().c_str());
							}
							ss.str("");
							ss << currDir << Config::EXE << " " << Config::CONFIG_FILE << " " << sizeIns << " " << indId + 1;
							//cout << ss.str() << endl;
							system(ss.str().c_str());
							if (Config::SOLVE_ONE)return 0;
							continue;
						}

						bool resolve = false;
						// Solve using the current code
						// Init before solving
						SolvedPb::ptrCurrIns = &Problem::currIns;
						Config::CURR_SIZE_INS = sizeIns;
						Config::MALLOC_FAILED = false;
						resetStatic();

						Alloc::init();
						ss.str(""); ss << "log" << sizeIns << "." << (indId + 1) << ".txt";
						auto logFileName = ss.str();
						logfile.open(logFileName, ios::out);

						Problem pb(sizeIns);
						FOR_E(ijob, sizeIns)pb.jobs[ijob] = Problem::currIns[ijob];
						pb.init();
					
						//!reuse the memory
						if (resolve) {
							cout << "Processing the memory for resolving...Current nbentry=" << Problem::pbMemo.nbentry << endl;
							Problem::pbMemo.removeUnusedAndUpdateHit();
							Problem::pbMemo.hit = Problem::pbMemo.miss = 0;
							cout << "After cleaning, useful nbentry=" << Problem::pbMemo.nbentry << endl;
							Config::MEM_LIM_STRATEGY = 0;
							Config::RAM_LIM_MO = 0;
							Config::CURR_SIZE_INS = sizeIns;
							Config::MALLOC_FAILED = false;
							Alloc::init();
							pb._isRLComputed = pb._isEddSorted = pb.N_rlist = 0;
							Problem::ctrAllNodes = 1;
							Problem::idNext = Problem::ctrMovL = Problem::ctrMerL = Problem::ctrSzMerL =
								Problem::ctrMovR = Problem::ctrMerR = Problem::ctrSzMerR = Problem::ctrJobsMemNotUsed =
								Problem::ctrKidsNotSingle = Problem::ctrNodeSingleKid = Problem::ctrNodeMultiKid = 0;
							Stat::resetStatic();
						}
						//pb.mmToldMe = 1;
						//cout << "Go...";
						int ttRes = 0;
						//clock_t time = clock();
						double timew = get_wall_time();
						double timec = get_cpu_time();
						//auto cyc = get_cpu_cycle();
						auto res = pb.solve(ttRes, short(Config::K));
						//double timeElapsed = double(clock() - time) / CLOCKS_PER_SEC;
						double timeElapsedW = get_wall_time() - timew;
						double timeElapsedC = get_cpu_time() - timec;
						//if(sizeIns>=900)Config::notify(to_string(sizeIns)+"."+to_string(indId+1)+"Finished. "+to_string(timeElapsedC)+"/"+to_string(timeElapsedW));
						//pb.pbMemo.DB_LIM = 1;
						//pb.pbMemo.removeOldDownToLim(false);
						//cout << "clock() gives "<<timeElapsed<<"; get_cpu_time = "<<timeElapsed2<<"; cycles = "<<get_cpu_cycle()-cyc<<endl;
						//logfile << "\n\nFinish\n\n";
						//Alloc::logAlloc();
						//cout << "Before reset.\n";
						//getchar();
						//Problem::pbMemo = Mem();
						//FOR_BE(i, Config::MEM_N_LB, Config::CURR_SIZE_INS)logfile << Problem::pbMemo.mapNbJobs[i].max_load_factor() << "\t" << Problem::pbMemo.mapNbJobs[i].load_factor() << "\t" << Problem::pbMemo.mapNbJobs[i].bucket_count()<< endl;
						//cout << "After reset\n";
						//getchar();

						//Log extrem nodes
						//if(logfile.is_open())logfile.close();
						//ss.str(""); ss << "log" << sizeIns << "." << (indId + 1) << ".txt";
						//logfile.open(ss.str(), ios::out);
						//for (auto it = Config::TmpMap.rbegin(); it != Config::TmpMap.rend(); ++it){
						//	logfile<<"\n"<< it->first<< "=======================\n";
						//	//sort(it->second.begin(), it->second.end());
						//	for (auto & p : it->second)
						//		logfile << get<0>(p) << "\t" << get<1>(p) << "\t" << get<2>(p) << "\t" << get<3>(p) << "\n";
						//	// Test dominance
						//	FOR_E(i, it->second.size())FOR_BE(j, i + 1, it->second.size()){
						//		if (get<2>(it->second[i]) <= get<2>(it->second[j]) && get<3>(it->second[i]) <= get<3>(it->second[j]))
						//			logfile << get<1>(it->second[i]) << " -> " << get<1>(it->second[j]) << endl;
						//	}
						//	logfile << endl << endl;
						//}

						// Log split count
						//for (int i = MAX_N_JOBS; i > 0; --i)
						//	if (Config::ctrNodes[i]>0)logfile << i << "\t" 
						//		<< Config::ctrSplit[i] << "\t" 
						//		<< Config::ctrNodes[i] << "\t"
						//		<< Config::ctrSplit[i]*100.0/Config::ctrNodes[i] << endl;
						
						// Remove the log file if it's empty
						bool toRemove = (logfile.tellp() <= 1);
						logfile.close();
						if (toRemove) { ::remove(logFileName.c_str()); }

#if defined(MEM_LOG_ML)
						ss.str(""); ss << "ml\\memlog_" << sizeIns << "_" << idIns << ".txt";
						ofstream logmem(ss.str());
						Problem::pbMemo.log(logmem);
#endif
						ttRes = Problem::TT(res, pb.N).first;
						//assert(ttRes == Problem::TTById(res, pb.N));
						tts.push_back(ttRes);
						result[indId].id = indId + 1;
						result[indId].tt = ttRes;
						result[indId].t_cpu = int(timeElapsedC);
						result[indId].t_wall = int(timeElapsedW);
						result[indId].nbmovL = Stat::ctrMemReplace;//Problem::ctrMovL;
						result[indId].nbmerL = Stat::ctrMemRemove; //Problem::ctrMerL;
						result[indId].szmerL = Stat::ctrMemRemove == 0 ? 0 : (Stat::szMemRemove / Stat::ctrMemRemove);// Problem::ctrMerL == 0 ? 0 : Problem::ctrSzMerL / Problem::ctrMerL;
						result[indId].nbmovR = Problem::ctrMovR;
						result[indId].nbmerR = Stat::ctrCutDb2;// Problem::ctrMerR;
						result[indId].szmerR = Stat::ctrCutDb2 == 0 ? 0 : Stat::szCutDb2 / Stat::ctrCutDb2;// Problem::ctrMerR == 0 ? 0 : Problem::ctrSzMerR / Problem::ctrMerR;
						result[indId].nbnode = Problem::ctrAllNodes;
						result[indId].hit = Problem::pbMemo.hit;
						result[indId].miss = Problem::pbMemo.miss;
						result[indId].nbentryMem = Problem::pbMemo.nbentry;
						result[indId].nbBytesMem = Problem::pbMemo.bytesMem;
						result[indId].nbJobsMem = -1; // (Problem::pbMemo.bytesMem - Problem::pbMemo.nbentry*sizeof(SolvedPb)) / 2;
						result[indId].ram = get_ram_usage();
						//Strategy
						result[indId].rmFail = Stat::ctrMemRemoveFail;
						result[indId].minRmRatio = Stat::minRmRatio;
						result[indId].info = Problem::pbMemo.isMemLimitReached ? "MemLimitReached" : "";
						//DB2 : a second database for passive memorization (store partial sequences)
						result[indId].hit2 = Problem::sigmaMemo.hit;
						result[indId].miss2 = Problem::sigmaMemo.miss;
						result[indId].nbentryMem2 = Problem::sigmaMemo.nbentry;
						result[indId].nbBytesMem2 = Problem::sigmaMemo.bytesMem;

						if (Config::MALLOC_FAILED)result[indId].info += string("/MallocFailed");
						cout << "Solved in " << setw(6) << result[indId].t_cpu << "/" << result[indId].t_wall << 
							"s. \nTT=" << setw(7) << result[indId].tt <<
							"\t#MovL=" << result[indId].nbmovL <<
							"\t#MerL=" << result[indId].nbmerL <<
							"\tszMerL=" << result[indId].szmerL <<
							"\t#MovR=" << result[indId].nbmovR <<
							"\t#MerR=" << result[indId].nbmerR <<
							"\nszMerR=" << result[indId].szmerR <<
							"\t#All=" << result[indId].nbnode <<
							"\tHits=" << result[indId].hit <<
							"\tMiss=" << result[indId].miss <<
							"\tNbEntryMem=" << result[indId].nbentryMem <<
							"\t#NbJobsMem=" << result[indId].nbJobsMem <<
							"\n#NbBytesMem=" << result[indId].nbBytesMem <<
							"\tRAM=" << result[indId].ram <<
							"\tMinRatio=" << result[indId].minRmRatio <<
							"\tHits2=" << result[indId].hit2 <<
							"\tMiss2=" << result[indId].miss2 <<
							"\tNbEntryMem2=" << result[indId].nbentryMem2 <<
							"\t#NbBytesMem2=" << result[indId].nbBytesMem2 << endl<<endl;
						
						// stat ctr
						//cout << "\nLeft ========================\n";
						//FOR_E(l, Stat::rangeCtr1.size())
						//	cout << "(" << l << ":" << Stat::rangeCtr1[l].first << "/" << Stat::rangeCtr1[l].second << ")\n";
						//cout << "\nRight ========================\n";
						//FOR_E(l, Stat::rangeCtr2.size())
						//	cout << "(" << l << ":" << Stat::rangeCtr2[l].first << "/" << Stat::rangeCtr2[l].second << ")\n";

						// Save to file immediately to avoid program interruption
						ss.str("");
						ss << outFolder << "res" << sizeIns << fileSuffix;
						saveCurrentRes(ss.str(), result[indId]);

						//if (!resolve){
						//	resolve = true;
						//	goto resolve;
						//}

						//cerr << "after save res\n";
						if (Config::SOLVE_ONE) {
							//notify();
							FREE<Job>(&res, pb.N);
							//cerr << "after Free res\n";
							return 0;
						}

#ifdef PRINT_SOL
						if (Config::SOLVE_ONE) {
							cout << "Sol: \n"; copy(res, res + sizeIns, ostream_iterator<Job>(cout, "\n"));
						}
#endif                
						FREE<Job>(&res, pb.N);
					}
					else
					{
						cout << "Instance : Already solved. Make sure res files were cleaned.\n";
					}
				}// For each instance
			}// For each t
		}// For each r
	
        
		ss.str("");
		ss << outFolder<<"stat" << sizeIns << fileSuffix;
		makeStat(ss.str(), result);

        // For test purpose, we compare tt values to that we obtain in the history. 
        ss.str("");
		ss << outFolder << "sol" << sizeIns << fileSuffix;
        
        ifstream solFileIn(ss.str());
        if (!solFileIn){    // Sol file does not exist, create it for later use
            ofstream solFileOut(ss.str());
            copy(tts.begin(), tts.end(), ostream_iterator<int>(solFileOut, "\n"));
            solFileOut.close();
        }
#ifdef VRIFY_RES
        else{// Compare current sols to that in the file
            int ttFile=-1;
            bool flag = true;
            for (auto itt = 0; itt < tts.size(); itt++){
                solFileIn >> ttFile;
				if (itt == 0) { 
					// Skip to the corresponding tt
					int nbskip = Config::ID_START-1;
					if (Config::ONLY_HARDEST && nbskip < 20) nbskip = 20;
					FOR_E(i, nbskip) solFileIn >> ttFile; 
				}
                if (tts[itt] != ttFile){
                    // ! problem
                    cerr << "TT value is inconsistent : file = " << ss.str() << "; instance id: " << itt + 1 << "; " << tts[itt]<<"!="<<ttFile << endl;
                    cerr << "Program suspended..." << endl;
                    cin >> ttFile;
                    flag = false;
                    //break;
                }
            }
            if (flag) cout << "\nSol consistency test passed!\n"<<endl;
            solFileIn.close();
        }
#endif
    }
    return 0;
}

void writeResHeader(string filename){
	ofstream fres(filename, ios::app);
	if (!fres)cerr << "Cannot open res file.\n", getchar();
	fres <<
		"id\t" <<
		"tt\t" <<
		"t_cpu\t" <<
		/*"nbReplace\t" <<
		"nbRm\t" <<
		"szRm\t" <<
		"unused\t" <<
		"nbCut\t" <<
		"szCut\t" <<*/
		"nbnode\t" <<
		"hit\t" <<
		"miss\t" <<
		"nbentryMem\t" <<
		/*"nbJobsMem\t" <<
		"nbBytesMem\t" <<*/
		"t_wall\t" <<
		"ram\t" <<
		/*"minRmRatio\t" <<
		"rmFail\t" <<
		"info\t" <<
		"hit2\t" <<
		"miss2\t" <<
		"nbentryMem2\t" <<
		"nbBytesMem2\t" <<*/ endl;
	fres.close();
}

void saveCurrentRes(string filename, const Result & res){
	ofstream fres(filename, ios::app);
	if (!fres)cerr << "Cannot open res file.\n", getchar();
	fres << res.id <<
		"\t" << res.tt <<
		"\t" << res.t_cpu <<
		/*"\t" << res.nbmovL <<
		"\t" << res.nbmerL <<
		"\t" << res.szmerL <<
		"\t" << res.nbmovR <<
		"\t" << res.nbmerR <<
		"\t" << res.szmerR <<*/
		"\t" << res.nbnode <<
		"\t" << res.hit <<
		"\t" << res.miss <<
		"\t" << res.nbentryMem <<
		/*"\t" << res.nbJobsMem <<
		"\t" << res.nbBytesMem <<*/
		"\t" << res.t_wall <<
		"\t" << res.ram << setprecision(2) << fixed <<
		/*"\t" << res.minRmRatio <<
		"\t" << res.rmFail <<
		"\t" << res.info <<
		"\t" << res.hit2 <<
		"\t" << res.miss2 <<
		"\t" << res.nbentryMem2 <<
		"\t" << res.nbBytesMem2 <<*/ endl;
	fres.close();
}

void makeStat(string filename, Result* result, int nbIns ,int nbInsPerRT ){
	// line 1:min avg max avg_nodes_moved avg_nodes_merged avgsizemerged avg_total_nodes;  
	// line n: R T min avg max avgmove avgnmerge avgsize avgnbnodes
	Stat stat[21];                          
	double sumT, minT, maxT, sumTperRT, minTperRT, maxTperRT;
	sumT = maxT = 0;
	minT = INT_MAX;
	long long nbMove = 0, nbMerge = 0, szMer = 0, 
		nbNode = 0, nbHit = 0, nbMiss = 0, nbMem = 0, nbJobsMem=0, nbBytesMem=0;

	for (int i = 0; i < 20; i++){
		sumTperRT = maxTperRT = 0;
		minTperRT = INT_MAX;
		long long nbMovePerRT = 0, nbMergePerRT = 0, szMerPerRT = 0, 
			nbNodePerRT = 0, nbHitPerRT = 0, nbMissPerRT = 0, nbMemPerRT = 0, nbJobsMemPerRT = 0, nbBytesMemPerRT=0;
		for (int j = 0; j < nbInsPerRT; j++)// For each instance
		{
			int indId = i * 10 + j;
			// Save results and do statistics
			nbMovePerRT += result[indId].nbmovL + result[indId].nbmovR;
			nbMergePerRT += result[indId].nbmerL + result[indId].nbmerR;
			szMerPerRT += (result[indId].szmerL*result[indId].nbmerL + result[indId].szmerR*result[indId].nbmerR);// / (result[indId].nbmerL + result[indId].nbmerR);
			nbNodePerRT += result[indId].nbnode;
			nbHitPerRT += result[indId].hit;
			nbMissPerRT += result[indId].miss;
			nbMemPerRT += result[indId].nbentryMem;
			nbJobsMemPerRT += result[indId].nbJobsMem;
			nbBytesMemPerRT += result[indId].nbBytesMem;
			sumT += result[indId].t_cpu;
			sumTperRT += result[indId].t_cpu;
			if (result[indId].t_cpu > maxT) maxT = maxTperRT = result[indId].t_cpu;
			else if (result[indId].t_cpu > maxTperRT) maxTperRT = result[indId].t_cpu;
			if (result[indId].t_cpu < minT) minT = minTperRT = result[indId].t_cpu;
			else if (result[indId].t_cpu < minTperRT) minTperRT = result[indId].t_cpu;
		}
		// stat for this RT pair
		int ind = i + 1;
		stat[ind].R = (i / 4 + 1)*0.2;				//R
		stat[ind].T = (i % 4 + 1)*0.2;				//T
		stat[ind].tmin = minTperRT;					//min
		stat[ind].tavg = sumTperRT / nbInsPerRT;	//avg
		stat[ind].tmax = maxTperRT;					//max
		stat[ind].nbmovAvg = nbMovePerRT / nbInsPerRT;    //avg nb nodes moved
		stat[ind].nbmerAvg = nbMergePerRT / nbInsPerRT;   //avg nb nodes merged
		// Be careful on this stat. Avg size of merged node should not be divided by nbInsPerRT
		stat[ind].sizemerAvg = (nbMergePerRT==0?0:(szMerPerRT / nbMergePerRT));  //avg size nodes merged
		stat[ind].nbnodeAvg = nbNodePerRT / nbInsPerRT;   //avg total nb nodes
		stat[ind].hitAvg = nbHitPerRT / nbInsPerRT;
		stat[ind].missAvg = nbMissPerRT / nbInsPerRT;
		stat[ind].nbentryMemAvg = nbMemPerRT / nbInsPerRT;
		stat[ind].nbJobsMemAvg = nbJobsMemPerRT / nbInsPerRT;   
		stat[ind].nbBytesAvg = nbBytesMemPerRT / nbInsPerRT;    
		nbMove += nbMovePerRT;
		nbMerge += nbMergePerRT;
		nbNode += nbNodePerRT;
		szMer += szMerPerRT;
		nbHit += nbHitPerRT;
		nbMiss += nbMissPerRT;
		nbMem += nbMemPerRT;
		nbJobsMem += nbJobsMemPerRT;
		nbBytesMem += nbBytesMemPerRT;
	}
	// Global stat for all instances of the current size
	stat[0].R = 0;
	stat[0].T = 0;
	stat[0].tmin = minT;				//min
	stat[0].tavg = sumT / nbIns;		//avg
	stat[0].tmax = maxT;				//max
	stat[0].nbmovAvg = nbMove / nbIns;  //avg nb nodes moved
	stat[0].nbmerAvg = nbMerge / nbIns; //avg nb nodes merged
	stat[0].sizemerAvg = (nbMerge==0)?0:(szMer / nbMerge);   //avg merged node size
	stat[0].nbnodeAvg = nbNode / nbIns; //avg total nodes
	stat[0].hitAvg = nbHit / nbIns;
	stat[0].missAvg = nbMiss / nbIns;
	stat[0].nbentryMemAvg = nbMem / nbIns;
	stat[0].nbJobsMemAvg = nbJobsMem / nbIns;
	stat[0].nbBytesAvg = nbBytesMem / nbIns;
	
	// Write stat to file
	ofstream out(filename);
	//out << stat[0] << " " << stat[1] << " " << stat[2] << " " << stat[3] << " " << stat[4] << endl;
	for (int k = 0; k < 21; k++){
		out << stat[k].R << " "
			<< stat[k].T << " "
			<< stat[k].tmin << " "
			<< stat[k].tavg << " "
			<< stat[k].tmax << " "
			<< stat[k].nbmovAvg << " "
			<< stat[k].nbmerAvg << " "
			<< stat[k].sizemerAvg << " "
			<< stat[k].nbnodeAvg << " "
			<< stat[k].hitAvg << " "
			<< stat[k].missAvg << " "
			<< stat[k].nbentryMemAvg<<" "
			<< stat[k].nbJobsMemAvg << " "
			<< stat[k].nbBytesAvg << endl;
	}
}

void resetStatic(){
	Problem::resetStatic();
	//memset(Problem::ctrNodes, 0, sizeof(long long)*(MAX_N_JOBS + 1));
	//memset(Problem::ctrSplit, 0, sizeof(long long)*(MAX_N_JOBS + 1));
	Stat::resetStatic();
}

//////////////////////////////////////////////////////////////////////////
// Measure REAL CPU Time 
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq)){
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)){
		//  Handle error
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0){
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		double kernelT = (double)(c.dwLowDateTime |
			((unsigned long long)c.dwHighDateTime << 32)) * 0.0000001;
		double userT =	(double)(d.dwLowDateTime |
			((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
		//cout << "Kernel = " << kernelT << ";\t User = " << userT<<endl;
		return kernelT + userT;
	}
	else{
		//  Handle error
		cerr << "GetProcessTimes returns 0. Current WallTime = " << get_wall_time() <<"\n";
		return 0;
	}
}
ULONG64 get_cpu_cycle(){
	ULONG64 ctr;
	//To compile an application that uses this function, define _WIN32_WINNT as 0x0600 or later.
	if (QueryProcessCycleTime(GetCurrentProcess(), &ctr)){
		return ctr;
	}
	else{
		//  Handle error
		return 0;
	}
}

long long get_ram_usage(){
	static  HANDLE currProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	if (K32GetProcessMemoryInfo(currProc, &pmc, sizeof(pmc)))
	{
		//printf( "\tWorkingSetSize: %llu\n", pmc.WorkingSetSize );
		//printf( "\tPagefileUsage (potential usage, not actual usage. Not 0 even with pagefile disabled...): %llu\n", pmc.PagefileUsage );		//!! When working set increases, pagefileusage also increases, which corresponds to "committed" column in task manager.
		//printf("\tPageFaultCount: %llu\n", pmc.PageFaultCount);
		//printf("\tPeakWorkingSetSize: %llu\n",
		//	pmc.PeakWorkingSetSize);
		//printf( "\tQuotaPeakPagedPoolUsage: %llu\n", 
		//	pmc.QuotaPeakPagedPoolUsage );
		//printf( "\tQuotaPagedPoolUsage: %llu\n", 
		//	pmc.QuotaPagedPoolUsage );
		//printf( "\tQuotaPeakNonPagedPoolUsage: %llu\n", 
		//	pmc.QuotaPeakNonPagedPoolUsage );
		//printf( "\tQuotaNonPagedPoolUsage: %llu\n", 
		//	pmc.QuotaNonPagedPoolUsage );
		//printf( "\tPeakPagefileUsage: %llu\n", pmc.PeakPagefileUsage );
		return pmc.WorkingSetSize;
	}
	return -1;
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time, NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
#endif

/////////////////////////////////////////////////////////
// Some experiments 
/////////////////////////////////////////////////////////
void testFree(){
	// Observation: when free, small blocks are not released to os immediately: try unit_size=400, below i%10000;
	//						   big blocks yes: try unit_size = 40960, below i%1000
	int n = 1000000;
	int unit_size = 400;	
	long long ** ptrs = new long long *[n];
	FOR_E(i, n){
		//ptrs[i] = (long long *)malloc(unit_size);
		ptrs[i] = new long long [unit_size/sizeof(long long)];
		memset(ptrs[i], 1, unit_size);
	}
	cout << "Init finished. Press enter to continue.\n";
	getchar();
	cout << "Start freeing...\n";
	int size_del = 0;
	for (int i = 0; i < n; i++){
		//free(ptrs[i]);
		delete (ptrs[i]);
		size_del += unit_size;
		if (i % 10000 == 0){
			cout << size_del << " Bytes Freed\n";
			size_del = 0;
			get_ram_usage();
			getchar();
		}
	}
	cout << "All freed";
	getchar();
}

void eatMemoryBy1G(){
	while (true){
		FOR_E(i, 1000){
			FAIL_ON(NULL == calloc(_1M, 1));
			cout << i + 1 << "M" << endl;
			get_ram_usage();
			getchar();
		}
		cout << "1G eaten. Enter s to stop.";
		char cmd = getchar();
		cout << cmd<<endl;
		if (cmd == 's')break;
	}
}

void getHardInfo(){
	SYSTEM_INFO siSysInfo;
	// Copy the hardware information to the SYSTEM_INFO structure. 
	GetSystemInfo(&siSysInfo);
	// Display the contents of the SYSTEM_INFO structure. 
	printf("  Page size: %u\n", siSysInfo.dwPageSize); // return 4096! 4K
}

void setMaxMemLimit(long long limit){
	auto job = CreateJobObject(NULL, NULL);
	FAIL_ON(job == NULL);
	FAIL_ON(0==AssignProcessToJobObject(job, GetCurrentProcess()));
	JOBOBJECT_BASIC_LIMIT_INFORMATION basicInfo;
	basicInfo.LimitFlags = JOB_OBJECT_LIMIT_PROCESS_MEMORY;
	JOBOBJECT_EXTENDED_LIMIT_INFORMATION extInfo;
	extInfo.BasicLimitInformation = basicInfo;
	extInfo.ProcessMemoryLimit = limit;
	FAIL_ON(0==SetInformationJobObject(job, JobObjectExtendedLimitInformation, &extInfo, sizeof(extInfo)));
}
//	BOOL b = SetProcessWorkingSetSize(GetCurrentProcess(),100000000,200000000); This is not a hard limit, system only attempts to respect it.