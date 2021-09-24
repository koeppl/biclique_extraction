#include <iostream>
#include <stdlib.h>
#include "AdjacencyMatrix.h"

int main(int argc, char **argv){
   int v=0;

   if(argc < 2 ){
	std::cout<<" usage: ./fromEu2Stream adjacencyMatrix.txt [0/1] (0 Claude Format, 1 adjacency format ). "<<std::endl;
	exit(1);
   }
   if(argc == 3){
	v = atoi(argv[2]);
   }

   AdjacencyMatrix adjMatrix;

   if(v == 0){
	adjMatrix.loadFileFranFormat(argv[1]); 
   } else if(v == 1){
	adjMatrix.loadFileAdjFormat(argv[1]); 
   }

   std::cout<<"number rows "<<adjMatrix.numberROWS<<std::endl;

   adjMatrix.dumpFile();

   return 0;

}
