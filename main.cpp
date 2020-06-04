#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "dna.h"
using namespace std;


int main() {

	string header = ">Header";
	string sequence = "CAGGTCGTTACCCC";

	DNA DNA1 = DNA(header, sequence);

	vector<string> result = DNA1.translate();

	for(int i = 0; i < (int)result.size(); i++ ){
	  cout<< result[i] <<endl;
	}

	return 0;
}
