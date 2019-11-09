#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <algorithm>

using namespace std;

template<typename ValueType>
struct Node
{
	Node()
	{
		label = '\0';	// default for root node
	}

	Node(char s)
	{
		label = s;	// label for child node
	}

	vector<ValueType> values;
	char label;
	vector<Node*> children;
};

template<typename ValueType>
class Trie
{
public:
	Trie()
	{
		Node<ValueType> *n = new Node<ValueType>();	// create a new empty node
		root = n;	// make this the root node		
	}
	~Trie()
	{
		destruct(root);	// call the destruct function to free all memory
	}
	void reset()
	{
		destruct(root);	// delete all memory of the Trie
		Node<ValueType> *n = new Node<ValueType>();	// give the Trie a new root
		root = n;
	}
	void insert(const string& key, const ValueType& value)
	{
		oldAdd(key, value, root);	// this should work		
	}

	vector<ValueType> find(const string& key, bool exactMatchOnly) const
	{
		if (exactMatchOnly)	// we only want exact matches
		{
			return rigidFind(key, root);
		}
		else	// snips are fine
		{
			return lenientFind(key, root);
		}
	}

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;

private:
	Node<ValueType>* root;
	void destruct(Node<ValueType>* n)	// post-order traversal
	{
		if (n == nullptr)
		{
			return;
		}


		int noOfElements = n->children.size();
		for (int i = 0; i < noOfElements; i++)
		{
			destruct(n->children[i]);
		}
		delete n;
	}

	void oldAdd(string s, ValueType val, Node<ValueType>* parent)	// if this is an existing string; always start with this
	{
		if (s.size() == 0)	// mistake has happened, return
		{
			return;
		}

		for (int i = 0; i < parent->children.size(); i++)	// checking if it exists
		{
			if (s[0] == parent->children[i]->label)	// does the first letter of the string equal a label of one of the children?
			{
				if (s.size() == 1)	// if only one letter is left
				{
					parent->children[i]->values.push_back(val);	// push it into the values vector of the child with the right label
					return;
				}

				return oldAdd(s.substr(1), val, parent->children[i]);	// otherwise just oldAdd, with one less letter, and using the node of the child
			}
		}

		// if program comes till here, it means nothing was found

		newAdd(s, val, parent);	// so we'll use a new add		
	}

	void newAdd(string s, ValueType val, Node<ValueType>* parent)	// if this is a new string 
	{
		if (s.size() == 0)	// mistake, return
		{
			return;
		}

		Node<ValueType> *n = new Node<ValueType>(s[0]);	// create a new node with label as the first letter of the string

		parent->children.push_back(n);	// push this new node into children of the current node

		if (s.size() == 1)	// if only one letter is left
		{
			n->values.push_back(val);	// push the required values into the new code
			return;
		}

		newAdd(s.substr(1), val, n);	// if more is left, reduce one letter, and use the newly made node as a parent

	}

	vector<ValueType> rigidFind(string key, Node<ValueType>* ptr) const	// only complete matches
	{
		vector<ValueType> result;	// empty vector

		if (key.size() == 0)	// mistake, return
		{
			return result;
		}

		for (int i = 0; i < ptr->children.size(); i++)	// checking if it exists
		{
			if (key[0] == ptr->children[i]->label)	// does the first letter of the string equal a label of one of the children?
			{
				if (key.size() == 1)	// if only one letter is left
				{
					result = ptr->children[i]->values;	// push it into the values vector of the child with the right label
					return result;
				}

				return rigidFind(key.substr(1), ptr->children[i]);	// otherwise just rigidFind, with one less letter, and using the node of the child
			}
		}
		return result;	// if we didn't find anything valid, return an empty vector
	}

	vector<ValueType> lenientFind(string key, Node<ValueType>* ptr, bool firstValid = false) const	// for SNIPS
	{
		vector<ValueType> result;

		if (key.size() == 0)	// mistake, return
		{
			return result;
		}

		for (int i = 0; i < ptr->children.size(); i++)	// checking if it exists
		{
			if (key[0] == ptr->children[i]->label)	// does the first letter of the string equal a label of one of the children?
			{
				if (key.size() == 1)	// if only one letter is left
				{
					vector<ValueType> toAdd = ptr->children[i]->values;
					result.insert(result.end(), toAdd.begin(), toAdd.end()); // add the values vector of this correct one to the current vector

				}
				else
				{
					vector<ValueType> extra = lenientFind(key.substr(1), ptr->children[i], true);	// this will perform a lenientFind, now with firstChar valid, taking the same string minus the first char, and the correct child node
					result.insert(result.end(), extra.begin(), extra.end());	// add this vector to our result vector
				}
			}

			else	// if the char doesn't exist
			{
				if (firstValid)	// if the first char was valid
				{
					if (key.size() == 1)	// if only one letter is left, just add anyways
					{
						vector<ValueType> toAdd = ptr->children[i]->values;
						result.insert(result.end(), toAdd.begin(), toAdd.end()); // add the current vector to result
					}
					else	// change to rigidFind, no more room for mistakes
					{
						vector<ValueType> extra = rigidFind(key.substr(1), ptr->children[i]);
						result.insert(result.end(), extra.begin(), extra.end());
					}
				}
			}
		}
		return result;
	}
};



#endif // TRIE_INCLUDED

