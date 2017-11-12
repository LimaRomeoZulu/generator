
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
	*output = NULL,
	*word = (char*)malloc(sizeof(char) * (8));
	char ch;	
	
	FILE 
		*fp = myfopen(output, "w"),
		*input = myfopen(GeneTreesPath, "r"),
		*treeFile;
	
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
	

	/* now see how many small trees we have */
	treeFile = getNumberOfTrees(tr, GeneTreesPath, adef);
	checkTreeNumber(tr->numberOfTrees, GeneTreesPath);
	
	//treeReadLen(fp, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);
	
	for(int i = 1; i <= tr->mxtips; i++){
		int count = 0;
	    	int read_counter;
		while((ch = treeGetCh(input)) == '('){}
		ungetc(ch,input);
		while(!feof(input))
	   	{
	       	read_counter = fscanf(input, "%[^,:();%f]%*[%f:.,();]",word);
			if(read_counter >= 0) 
			{
				if(strcmp(word, tr->nameList[i])==0)	count++;
			}
	    }
	    rewind(input);
		fprintf(fp,"%s, %f \n", tr->nameList[i], (float)count/tr->numberOfTrees);
	}
	fclose(fp);
	fclose(input);
	//free allocated memory
	rax_free(tr);
	rax_free(adef);
	rax_free(rdta);
	rax_free(cdta);
	return 0;
}
