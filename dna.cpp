// Name: Fehmina Hasan
// Date: 25th October 2017


#include <string>
#include <vector>
#include <cstdlib>
#include<fstream>
#include <iostream>
#include "dna.h"
using namespace std;


//default Constructor
DNA::DNA(){
  string hdr = "";  // hdr is assigned
  string seq = "";  // seq is assigned as an empty string
}

//Constructor with separate arguments for header and sequence
DNA::DNA(string header,string sequence){
  if (header[0] == '>') {
    hdr = header;
  }
  else throw runtime_error("Invalid");  

  //for loop to go through each sequence character
  for (int i=0;i < (int) sequence.size();i++){  
      if (sequence[i] == 'A' || sequence[i] == 'T' || sequence[i] == 'G' || sequence[i] == 'C'){
        seq += sequence[i];
      }
      else throw runtime_error("Invalid sequence");  
  }
}

//Constructor that expects an input file stream to a FASTA-formatted file
DNA::DNA(ifstream &infile){
  string header;
  getline(infile, header);
  string sequence;
  string line;
  while (getline(infile, line)) { //while loop
    sequence += line;
  }
	if(header[0] =='>') {
  }
  else throw runtime_error("Invalid output");  //exception handling
  for(int i = 0; i < (int) sequence.size(); i++){
		if(sequence[i] == 'A' || sequence[i] == 'G' || sequence[i] == 'T' || sequence[i] == 'C'){
    }
      else throw runtime_error("Invalid");   //exception handling
		}
  hdr = header;
  seq = sequence;
  infile.close();
}

//Getter methods that returns seq
string DNA::getSequence(){ 
  return seq;
}

//getter methods that returns hdr
string DNA::getHeader(){ 
  return hdr;
}

//Function that gives result in the fasta format by taking 3 characters at a time
string DNA::toFasta(int columns){
  string header = hdr;   
  string sequence = "";  
	int i = 0;
	while (i < (int)seq.length()) {  
		if (i + columns < (int)seq.length()) {
			sequence = sequence + seq.substr(i, columns) + "\n";
			i = i + columns;
		}
    else sequence = sequence + seq.substr(i) + "\n";
	}
	return header + "\n" + sequence + "\n";  
}

//this function gives the reverse compliment of the original seq
DNA DNA::revcomp(){
  string outcome = ""; 
  // for loop that goes from the last value to the first
	for (int a = (int) seq.length()-1; a >= 0; a--) {  
    if (seq[a] == 'A') outcome += 'T';
    if (seq[a] == 'T') outcome += 'A';
    if (seq[a] == 'C') outcome += 'G';
    if (seq[a] == 'G') outcome += 'C';
    if (seq[a] == 'N') outcome += 'N';
	}
	return DNA(hdr, outcome); 
}

//function Searches the sequence for `query`, a string, AND its reverse complement.
size_t DNA::find(string query, size_t start){

  //calls on the revcomp and getSequence function for the string
	string compiling = DNA(">", query).revcomp().getSequence(); 

   //it finds the postion of the 'query' from the start
	size_t orig = seq.find(query, start); 

  //it finds the postion of the 'compiling' from the start
	size_t backwards = seq.find(compiling, start); 
  if (backwards < orig) return backwards;

  if (orig < backwards) return orig;
  return string::npos;   //return statement
}


//function that converts the codons into amino acids and called in the translate function
string HelperFunction2(string words){
   //if statements used to convert to single letters
  if(words == "ATG") return "M"; 
  if(words == "TGT") return "W";
  if(words == "TTT" || words == "TTC") return "F";
  if(words == "TAT" || words == "TAC") return "Y";
  if(words == "CAT" || words == "CAC") return "H";
  if(words == "CAA" || words == "CAG") return "Q";
  if(words == "AAT" || words == "AAC") return "N";
  if(words == "AAA" || words == "AAG") return "K";
  if(words == "GAT" || words == "GAC") return "D";
  if(words == "GAA" || words == "GAG") return "E";
  if(words == "TGT" || words == "TGC") return "C";
  if(words == "AGT" || words == "AGC") return "S";
  if(words == "ATT" || words == "ATC" || words == "ATA") return "I";
  if(words == "GTT" || words == "GTC" || words == "GTA" || words == "GTG") return "V";
  if(words == "TCT" || words == "TCC" || words == "TCA" || words == "TCG") return "S";
  if(words == "CCT" || words == "CCC" || words == "CCA" || words == "CCG") return "P";
  if(words == "ACT" || words == "ACC" || words == "ACA" || words == "ACG") return "T";
  if(words == "GCT" || words == "GCC" || words == "GCA" || words == "GCG") return "A";
  if(words == "GGT" || words == "GGC" || words == "GGA" || words == "GGG") return "G";
  if(words == "CGT" || words == "CGC" || words == "CGA" || words == "CGG" || words == "AGA" || words == "AGG") return "R";
  if(words == "TTA" || words == "TTG" || words == "CTT" || words == "CTC" || words == "CTA" || words == "CTG") return "L";
  else return "STOP";
}


//function takes in the 3 character codons and converts them into amino acids and returns it in a new string
vector<string> DNA::translate(){
  string newseq = "";
  vector<string> list_of_new_seq_to_amino_acid;
  for(int i = 0; i < 3; i++){ //for loop
    string codons = "";   // new string
    for(int j = i; j < (int) seq.length(); j+=3){  
      newseq = seq.substr(j,3);
      string newseqcodon = HelperFunction2(newseq);
      if (newseqcodon == "STOP"){
        break;
      }
      else{
        codons += newseqcodon;
      }
    }
    if (codons != "STOP") list_of_new_seq_to_amino_acid.push_back(codons);
  }
  string revnewseq = revcomp().getSequence(); 
  for(int i = 0; i < 3; i++){ //for loop
    string codons = "";  // new string
    for(int j = i; j < (int) revnewseq.length();j+=3){ 
      string NRnewseq = revnewseq.substr(j,3);
      string newrevcodon = HelperFunction2(NRnewseq);
      if (newrevcodon == "STOP"){
        break;
      }
      else{
        codons += newrevcodon;
      }
    }
    if(codons != "STOP") list_of_new_seq_to_amino_acid.push_back(codons);
  }
  return list_of_new_seq_to_amino_acid;
}


//boolean function that gives the return statement if the d1 seq is the same as d2 seq
bool operator==(DNA d1, DNA d2){
  return d1.seq == d2.seq;  
}
