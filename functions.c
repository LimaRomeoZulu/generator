#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <inttypes.h>
#include <getopt.h>

#include "RAxML/globalVariables.h"


int countTaxaInTopology(char fileName[1024])
{
  FILE 
    *f = myfopen(fileName, "rb");   

  int
    c,   
    taxaCount = 0;

  while((c = fgetc(f)) != EOF)
    {
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    ungetc(c, f);
	  else
	    {	      	      	  	      
	      do
		{		
		  c = fgetc(f);
		}
	      while(c != ':' && c != ')' && c != ',');	    

	      taxaCount++;	     	     
	    
	      ungetc(c, f);
	    }
	}
    }
 
  printf("Found a total of %d taxa in tree file %s\n", taxaCount, fileName);

  fclose(f);

  return taxaCount;
}

int* readHistogram(char fileName[1024])
{
  FILE 
    *f = myfopen(fileName, "rb");
  
  int 
	item = 0,
	count = 0,
  	arraySize = 1024;
	
	int *result;
	
	result = (int *) malloc(arraySize * sizeof(int *));
	
	while(fscanf(f,"%d",&item) == 1)  
	{  
		if(count<arraySize){
			result[count] = item;
		}
		else{
			arraySize *= 2;
			result = (int *) realloc(result, arraySize * sizeof(int *));
		}
	} 
	return result;
}

void copyArgs(char **targetVariable, char *sourceArgv)
{
	
	size_t length = strlen(sourceArgv)+1;
	*targetVariable = (char *)rax_malloc(length * sizeof(char *));
    memcpy(*targetVariable, sourceArgv, length);
}






