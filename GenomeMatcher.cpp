#include "provided.h"
#include "Trie.h"	// new
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

struct geneIndex : public DNAMatch
{
	int index;
};

class GenomeMatcherImpl
{
public:
	GenomeMatcherImpl(int minSearchLength);
	void addGenome(const Genome& genome);
	int minimumSearchLength() const;
	bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
	bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	int smallestSearch;
	vector<Genome> geneLibrary;
	Trie<geneIndex> genomeTrie;
};

bool geneMatchSort(GenomeMatch g1, GenomeMatch g2)
{
	if (g1.percentMatch != g2.percentMatch)
	{
		return g1.percentMatch > g2.percentMatch;
	}
	else
	{
		return g1.genomeName < g2.genomeName;
	}
}

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	smallestSearch = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	int seqLength = genome.length();
	string genomeExtract;
	for (int i = 0; i + smallestSearch <= seqLength; i++)
	{
		geneIndex d;
		d.genomeName = genome.name();
		genome.extract(i, smallestSearch, genomeExtract);
		d.position = i;
		d.length = smallestSearch;
		d.index = geneLibrary.size();
		genomeTrie.insert(genomeExtract, d);
	}	
	geneLibrary.push_back(genome);
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return smallestSearch;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if ((fragment.size() < minimumLength) || minimumLength < smallestSearch)
	{
		return false;
	}
	string smallestString = fragment.substr(0, minimumLength);	//the smallest string that must be correct of fragment

		vector<geneIndex> curMatch = genomeTrie.find(smallestString.substr(0, smallestSearch), exactMatchOnly);	// all geneIndexes who have have smallestSearch chars correct
		
		if (curMatch.empty())	// no small matches found, return false
		{
			return false;
		}

		for (int i = 0; i < curMatch.size(); i++)	// the loop to check the maximum length of each match
		{
			geneIndex* cur = &curMatch[i];	// point to a GeneIndex from curMatch
			const Genome* g = &geneLibrary[cur->index];	// extract its genome ptr
			for (int j = fragment.size(); j >= minimumLength ; j--)	// loop to check its maximum length, with j causing length to reduce until a true match is found
			{	
				// int i = cur->position;
				string uselessString;	// a string to add the resultant string of extract to
				if (g->extract(cur->position, j, uselessString))	// if this isn't true, then move on, otherwise, store the resultant string in useless string // this causes the crash
				{
					string smallFrag = fragment.substr(0, j);	// a substring of the same size as uselessString taken from the fragment
					if (uselessString == smallFrag)	// if they are equal, set length of the current geneIndex as j(highest correct length) and move on to next GeneIndex
					{
						cur->length = j;
						break;
					}
					else if (!exactMatchOnly)	// else, if a SNiP is acceptable
					{
						if (uselessString[0] != smallFrag[0])	// if the starting chars aren't the same, it can't be a SNiP, move on to next GeneIndex
						{
							break;
						}
						int mistakes = 0;	// to count the number of errors
						for (int s = 0; s < uselessString.size(); s++)	// a loop to check every corresponding char
						{
							if (uselessString[s] != smallFrag[s])	// if two chars in the same postion aren't equal, increase mistakes
							{
								mistakes++;
							}
						}
						if (mistakes <= 1)	// if mistakes = 0,1, set length of the current geneIndex as j(highest correct length) and move on to next GeneIndex
						{
							cur->length = j;
							break;
						}
					}
				}
			}
		}
		
		vector<DNAMatch> finalMatch;	// a vector with all DNA matches that have the correct length, and excluding duplicates
		for (int i = 0; i < curMatch.size(); i++)	// checking every geneIndex of curMatch
		{
			DNAMatch d = curMatch[i];	// extracting the DNAMatch part of the geneIndices in curMatch, since we don't need the Genome ptr
			DNAMatch* dup = nullptr;	// to check for duplicates
			for (int j = 0; j < finalMatch.size(); j++)
			{
				if (d.genomeName == finalMatch[j].genomeName)	// if two genomes have the same names, they are duplicates
				{
					dup = &(finalMatch[j]);
					break;
				}
			}
			if (dup == nullptr)	// if no duplicates were found, and the length of the DNAMatch is bigger than the minimumLength specified, add it to finalMatch
			{
				if (d.length >= minimumLength)
				{
					finalMatch.push_back(d);
				}
			}
			else	// duplicate found
			{
				if (d.length > dup->length)	// if the length of the current DNAIndex is larger than its duplicate's, replace the duplicates length and position with that of the current DNAMatch
				{
					dup->length = d.length;
					dup->position = d.position;
				}
				else if (d.length == dup->length)	// if they both have the same length, then choose the one that appeared first
				{
					if (d.position < dup->position)
					{
						dup->position = d.position;
					}
				}
			}
		}

		if (finalMatch.empty())	// if this final vector is empty, return false
		{
			return false;
		}

		matches.swap(finalMatch);	// replace the contents of matches with those of finalMatch

	return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < this->minimumSearchLength())
	{
		return false;
	}

	int noOfSequences = (query.length()) / fragmentMatchLength;

	vector<DNAMatch> matchVector;
	
	for (int i = 0; i < noOfSequences; i++)
	{
		vector<DNAMatch> minimatchVector;
		string seq;
		query.extract(i*fragmentMatchLength, fragmentMatchLength, seq);
		this->findGenomesWithThisDNA(seq, fragmentMatchLength, exactMatchOnly, minimatchVector);
		matchVector.insert(matchVector.end(), minimatchVector.begin(), minimatchVector.end());
	}

	if (matchVector.empty())
	{
		return false;
	}

	vector<GenomeMatch> percentVector;

	for (int i = 0; i < geneLibrary.size(); i++)
	{
		Genome curGene = geneLibrary[i];
		GenomeMatch g;		
		g.genomeName = curGene.name();
		double noOfMatches = 0;
		for (int j = 0; j < matchVector.size(); j++)
		{
			if (g.genomeName == matchVector[j].genomeName)
			{
				noOfMatches++;
			}
		}
		double percent = (noOfMatches * 100) / noOfSequences;
		g.percentMatch = percent;
		if (percent >= matchPercentThreshold)
		{
			percentVector.push_back(g);
		}
	}
	sort(percentVector.begin(), percentVector.end(), geneMatchSort);
	results.swap(percentVector);
	return true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
	m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
	delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
	m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
	return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
