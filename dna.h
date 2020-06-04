#ifndef _dna_h
#define _dna_h

#include <iostream>
#include <string>
#include <vector>

/*
A class to represent DNA sequences, related to the FASTA format.
A DNA object has two member fields: a header, which is a sequence of characters beginning with the `>` character,
and a sequence, which is a sequence of characters containing only A, C, T, G, and N.
*/

class DNA {

public:

	/*
	Default constructor
	*/
	DNA();

	/*
	Constructor with separate arguments for header and sequence
	If `sequence` contains invalid characters, throws a std::runtime_error
	*/
	DNA(std::string header, std::string sequence);

	/*
	Constructor that expects an input file stream to a FASTA-formatted file
	If the first line is not a valid header (does not begin with a `>`)
	or if the sequence part of the file contains invalid characters
	(other than newlines, which are removed),
	throws a std::runtime_error
	*/
	DNA(std::ifstream &infile);

	/*
	Getter methods
	*/
	std::string getSequence();
	std::string getHeader();

	/*
	Produces FASTA-formatted output as a string,
	with the sequence wrapped to `columns` number of characters (default 80)
	*/
	std::string toFasta(int columns=80);

	/*
	Returns a new DNA object whose `header` is the same as this object's,
	and whose underlying `sequence` is the reverse complement as this object's.
	*/
	DNA revcomp();

	/*
	Searches the sequence for `query`, a string, AND its reverse complement.
	Returns the index of the first instance of `query` or its reverse complement,
	or `string::npos` if not found.
	*/
	size_t find(std::string query, size_t start);

	/*
	Generates all six possible amino acid sequences from this sequence.
	Returns them as a vector of strings, in unspecified order.
	*/
	std::vector<std::string> translate();

private:

	/*
	the equality operator does not test full equality.
	It returns true if and only if d1 and d2 have the same sequence
	(even if their headers differ), and false otherwise.
	*/
	friend bool operator==(DNA d1, DNA d2);
	

	/* instance variables */
	std::string seq;
	std::string hdr;

};


#endif
