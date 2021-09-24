#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "Trie.h"

using namespace std;

vector<int> splitToInt(string line, string delims){
   string::size_type bi, ei;
   vector<int> words;
   bi = line.find_first_not_of(delims);
   while(bi != string::npos) {
        ei = line.find_first_of(delims, bi);
        if(ei == string::npos)
                ei = line.length();
                string aux = line.substr(bi, ei-bi);
                words.push_back(atoi(aux.c_str()));
                bi = line.find_first_not_of(delims, ei);
        }
  return words;
}

void readFile(ifstream *infile, map<uint32_t, vector<uint32_t> > &adjlist){

   string line, word;

   while(getline(*infile,line)){
        //cout<<"line "<<line<<endl;
        if(line.size() == 0)
                continue;
        vector<int> words = splitToInt(line, " ");
	vector<uint32_t> outlink;
	uint32_t currentNode = (uint32_t)words[0];
	for(unsigned j=1; j<words.size(); j++){
		outlink.push_back((uint32_t)words[j]);
	}
	sort(outlink.begin(), outlink.end());
	adjlist[currentNode] = outlink;
   }
}


int main(int argc, char *argv[]){

  if(argc < 2){
	cout<<"usage: ./testtrie file \n";
	exit(1);
  }

  map<uint32_t, vector<uint32_t> > adjlist;
  map<uint32_t, vector<uint32_t> >::iterator mit;

  ifstream *infile = new ifstream(argv[1]);
  if(!infile){
       cout<<" file "<<argv[1]<<" does not exists"<<endl;
       exit(1);
  }


  readFile(infile, adjlist);

  Trie t;
  for(mit=adjlist.begin(); mit != adjlist.end(); mit++){
      vector<uint32_t> outlink = mit->second;
      uint32_t v = mit->first;
      t.insert(v,outlink);
  }
  t.printTrie();
  t.makeEmpty();

}
