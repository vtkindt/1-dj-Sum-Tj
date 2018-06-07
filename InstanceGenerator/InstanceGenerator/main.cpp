#include<iostream>
#include<fstream>
#include<random>
#include<iterator>
#include<functional>
#include<algorithm>
#include<numeric>
#include<string>
#include<sstream>
#include <cstdlib>
#include <ctime>
#define EPSILON 0.00000001
using namespace std;



// Generate random instances and write them to folder "[instance size]/"
// The total number of instances is nbInsPerRT*20 since we have 20 combinations of RT (see Potts 1982)
void Generate(int nbJobStart=100, int nbInsPerRT=10, int nbJobStep=100, int nbJobEnd=600, double r=-1, double t=-1)
{
	string folder = "data";
	//system(("MD " + folder).c_str());

    default_random_engine randgen(time(NULL));
    uniform_int_distribution<int> dist(1,100);
	//srand(time(NULL));
	auto genPRand = [&](){return dist(randgen); };
	stringstream ss;
    for(int n=nbJobStart; n<=nbJobEnd; n+=nbJobStep){
		system(("MD " + folder+"\\"+to_string(n)).c_str());
		       
		vector<int> vecP(n), vecD(n);
        for(double R=0.2; R<=1.0; R+=0.2){
			if (r > 0 && abs(R-r)>EPSILON)continue;
            for(double T=0.2; T<=0.8; T+=0.2){
				if (t > 0 && abs(T - t)>EPSILON)continue;
				for(int i=0; i<nbInsPerRT; i++){
					ss.str(string());
					ss << folder << "/"<<n << "/SDT_"<<n<<"_"<<R<<"_"<<T<<"_"<<(i+1)<<".txt";
					ofstream file(ss.str());


                    // nbInsPerRT instances per combination RT
					generate(vecP.begin(), vecP.end(), genPRand);
					//generate(vecP.begin(), vecP.end(), genPRand);
					int sumP = accumulate(vecP.begin(), vecP.end(), 0);
					int minD = (sumP*(1 - T - R / 2)), maxD = (sumP*(1 - T + R / 2));
					uniform_int_distribution<int> rdist(minD, maxD);
					//auto genDRand = [minD, maxD]() {int res = rand() % (maxD - minD + 1) + minD; return res > 0 ? res : 0; };
					auto genD = [&]() { return max(0, rdist(randgen)); };
					//generate(vecD.begin(), vecD.end(), genDRand);
					generate(vecD.begin(), vecD.end(), genD);
                    // To be harder: sort lpt + edd
                    //sort(vecP.begin(), vecP.end(), greater<int>());
                    //sort(vecD.begin(), vecD.end());
                    // Write instance
                    //file<<i+1<<endl;
					cout << R<<","<<T<<", Instance " << i + 1 << endl;
                    for(int j=0; j<n; j++){
                        file<<vecP[j]<<" "<<vecD[j]<<endl;
                    }
					file.close();
                }
            }
        }
    }
}



int main(){
	//int nbJobStart=100, int nbInsPerRT=10, int nbJobStep=100, int nbJobEnd=600, double r=-1, double t=-1
    //Generate(700,10,100,1000);
	//Generate(200,200,1,200, 0.2, 0.6); //generate for ml
	Generate(100, 10, 100, 1300);
    return 0;
}