
#include <unistd.h>
#include "RAxML/axml.h"
#include "functions.c"
#include <stdio.h>
#include <math.h>


int main(int argc, char* argv[]) {
	
	rawdata			*rdta;
	cruncheddata *cdta;
	tree				 *tr;
	analdef			*adef;
	
	char *SpeciesTreePath = NULL,
	*GeneTreesPath = NULL,
	*output = NULL;
	//std::ofstream evaluationTrees;
	
	copyArgs(&SpeciesTreePath, argv[1]);
	copyArgs(&GeneTreesPath, argv[2]);
	copyArgs(&output, argv[3]);

	adef = (analdef *)rax_malloc(sizeof(analdef));
	rdta = (rawdata *)rax_malloc(sizeof(rawdata));
	cdta = (cruncheddata *)rax_malloc(sizeof(cruncheddata));
	tr	 = (tree *)rax_malloc(sizeof(tree));
	initAdef(adef);
	
	adef->restart = TRUE;
	//strcpy(bootStrapFile, optarg);
	adef->model = M_GTRCAT;
	adef->useInvariant = FALSE;
	adef->readTaxaOnly = TRUE;
	
	printf("%s\n","Reference Tree");
	extractTaxaFromTopology(tr, rdta, cdta, SpeciesTreePath); 
	getinput(adef, rdta, cdta, tr);
	checkOutgroups(tr, adef);	
	
	FILE 
		*fp = myfopen(output, "w"),
		*input = myfopen(GeneTreesPath, "r");
	
	//treeReadLen(fp, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);
	
	for(int i = 0; i < sizeof(tr->nameList); i++){
		int count = 0;
	    int read_counter;
		while(!feof(input))
	    {
	        fscanf(input, "%s",tr->nameList);
			if(read_counter >= 0) count++;
	    }
	    count--;
	    rewind(input);
		fprintf(fp,"%s, %d", tr->nameList[i], count);
	}
	myfclose(fp);
	//free allocated memory
	rax_free(SpeciesTreePath);
	rax_free(GeneTreesPath);
	rax_free(tr);
	rax_free(adef);
	rax_free(rdta);
	rax_free(cdta);
	
	return 0;
}
