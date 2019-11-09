#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
	GenomeImpl(const string& nm, const string& sequence);
	static bool load(istream& genomeSource, vector<Genome>& genomes);
	int length() const;
	string name() const;
	bool extract(int position, int length, string& fragment) const;
private:
	string gname;
	string seq;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
	this->gname = nm;
	this->seq = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
	if (!genomes.empty())
	{
		genomes.erase(genomes.begin(), genomes.end());
		genomes.pop_back();
	}

	string istreamString = "";	// the string that will hold each getline
	string newName = "";	// a string for the name of the genes
	string newGene = "first";	// a string for the genome code of the genes
	while (getline(genomeSource, istreamString))
	{
		if (istreamString.empty())	// if an empty string is found, return false
		{
			return false;
		}

		if (istreamString[0] == '>')	// if the first char is >
		{
			if (istreamString.size() == 1 || newGene == "")	// if that's the only char in the line, or there isn't anything in the gene
			{
				return false;
			}

			if (newGene != "first")	// if this isn't the first occurence of '>', then add the previous name and gene to the list and push it back
			{
				Genome* g = new Genome(newName, newGene);
				genomes.push_back(*g);
			}

			newName = istreamString.substr(1);	// set new name to whatever is in this line, minus >
			newGene = ""; // clear up new gene
		}
		else
		{
			if (newGene == "first")	// if the first line ends up being one that isn't a name line
			{
				return false;
			}

			for (int i = 0; i < istreamString.size(); i++)	// this method to check every character of each line
			{
				char c = toupper(istreamString[i]);
				switch (c)
				{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					newGene += c;	// add the character to newGene
					break;
				default:
					return false;
				}
			}
		}
	}

	Genome* g = new Genome(newName, newGene);
	genomes.push_back(*g);
	return true;
}

int GenomeImpl::length() const
{
	return seq.size();
}

string GenomeImpl::name() const
{
	return gname;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	int geneSize = this->length();
	if ((position + length) > geneSize)
	{
		return false;
	}
	fragment = seq.substr(position, length);
	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
	m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
	delete m_impl;
}

Genome::Genome(const Genome& other)
{
	m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
	GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
	delete m_impl;
	m_impl = newImpl;
	return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
	return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
	return m_impl->length();
}

string Genome::name() const
{
	return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
	return m_impl->extract(position, length, fragment);
}
