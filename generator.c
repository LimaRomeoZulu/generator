
#include <unistd.h>
//extern "C" {
#include "RAxML/axml.h"
//}
#include "functions.c"
#include <stdio.h>
#include <math.h>


void generateGeneTree(tree *referenceTree, tree *geneTree, int taxaNewSpeciesTree, int taxaGeneTree, analdef *adef){
	// Create result tree.
	nodeptr	
		p,
	tmp;
	
	int 
		ch, 
	leftChild,
	parent,
	lcount = 0; 
	boolean
	readBranches 		= FALSE,
	storeBranchLabels	= FALSE;
		
	setupGeneTree(geneTree, taxaNewSpeciesTree);
	
	geneTree->start			= geneTree->nodep[taxaNewSpeciesTree];
	geneTree->ntips			= 0;
	geneTree->nextnode		= geneTree->mxtips + 1;	
	
	//Copy the prefixSum Array to alter it for leafs that were already added
	geneTree->taxaOccurencePrefixSum = longDup(referenceTree->taxaOccurencePrefixSum, referenceTree->mxtips);
	determineLeafs(geneTree, taxaGeneTree, referenceTree->numberOfTrees);
	
	p = geneTree->nodep[(geneTree->nextnode)]; 
	
	leftChild = addElement(referenceTree->start->back, geneTree, p, readBranches, &lcount, adef, storeBranchLabels);
	float draw = (rand() / ((double)RAND_MAX + 1.0));
        //add taxa to tree
        if(draw <= prob)
        {
		parent = referenceTree->start->number;        	
		(geneTree->ntips)++;
               
        }		
	p = geneTree->nodep[leftChild];
	if (leftChild > 0 && parent > 0 && geneTree->ntips > 2)
	{
		tmp = geneTree->nodep[parent];
		geneTree->start = tmp;
		hookupDefault(p, tmp, geneTree->numBranches);
	}
	else if(leftChild == 0 && parent > 0)
	{
		printf("ERROR: Only one leaf\n");
		return;
	}
	else if(leftChild > 0 && parent == 0)
	{
		hookupDefault(p, p->next->next->back, geneTree->numBranches);
		geneTree->rooted = TRUE;
		geneTree->wasRooted		 = TRUE;
		p->next->next->back = (nodeptr) NULL;
	}
	
	return;
}

int main(int argc, char* argv[]) {
	
	rawdata			*rdta;
	cruncheddata 	*cdta;
	tree			*tr;
	analdef			*adef;
	//char modelChar = 'R';
	
	char *referenceSpeciesTreePath = NULL,
	*newSpeciesTreePath = NULL,
	*geneTreePath = NULL;
	//*histogramPath = NULL,
	//*rfMetrixPath = NULL;
	//std::ofstream evaluationTrees;

	int innerBranches = 0,
	innerNodeNumber = 0;
	
	copyArgs(&referenceSpeciesTreePath, argv[1]);
	copyArgs(&newSpeciesTreePath, argv[2]);
	copyArgs(&geneTreePath, argv[3]);
	//copyArgs(&histogramPath, argv[3]);
	//copyArgs(&rfMetrixPath, argv[4]);
	
	//int *histogram = NULL; 
	//double *rfMetrix = NULL; 
	
	//histogram = readHistogram(histogramPath);
	//rfMetrix = readRFDistance(rfMetrixPath);

	
	//TODO: check where it is freed
	adef = (analdef *)rax_malloc(sizeof(analdef));
	rdta = (rawdata *)rax_malloc(sizeof(rawdata));
	cdta = (cruncheddata *)rax_malloc(sizeof(cruncheddata));
	tr	 = (tree *)rax_malloc(sizeof(tree));
	initAdef(adef);
	
	adef->restart = TRUE;
	adef->model = M_GTRCAT;
	adef->useInvariant = FALSE;
	adef->readTaxaOnly = TRUE;
	
	printf("%s\n","Reference Tree");
	extractTaxaFromTopology(tr, rdta, cdta, newSpeciesTreePath); 
	getinput(adef, rdta, cdta, tr);
	checkOutgroups(tr, adef);	

	int taxaNewSpeciesTree = rdta->numsp;
	int taxaReferenceSpeciesTree = countTaxaInTopology(referenceSpeciesTreePath);
	
	FILE 
		*fp = myfopen(newSpeciesTreePath, "r"),
		*output = myfopen("../geneTrees.tre", "w");
	
	treeReadLen(fp, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);
	rewind(fp);
	//srand(time(NULL));
	int
		*taxonHasDeg = (int *)rax_calloc((2*tr->mxtips - 2),sizeof(int)),
		*taxonToLabel  = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),
		*eulerIndexToLabel = (int *)rax_malloc((4*tr->mxtips - 5) * sizeof(int)),
		*labelToTaxon = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),
		*taxonToEulerIndex  = (int *)rax_malloc((tr->mxtips) * sizeof(int)),
			*taxonToReduction = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
	prepareRefernceTree(tr, taxonToReduction, taxonHasDeg, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);
	getGeneTreeStatistics(tr, geneTreePath, adef, taxonToReduction, taxonHasDeg, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);

	int i = 0;
	double rf;
	Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
	printf("%s", tr->tree_string);
	for(i; i < tr->numberOfTrees; i++){
		tree *geneTree = (tree *)rax_malloc(sizeof(tree));
		int taxaGeneTree = tr->geneLeafDistributions[i];
		innerNodeNumber =0;
		innerBranches = 0;
		
		//init HashTable with all TaxaNames because we don't know at the moment which taxa won't be in the gene tree
		geneTree->nameHash 	= tr->nameHash;		
		geneTree->rdta 		= tr->rdta;
		geneTree->nameList	= tr->nameList;		
		geneTree->rooted	= FALSE;

		generateGeneTree(tr, geneTree, taxaNewSpeciesTree, taxaGeneTree, adef);
		Tree2String(geneTree->tree_string, geneTree, geneTree->start->back, FALSE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
		printf("%s", geneTree->tree_string);
		
		innerNodeNumber = geneTree->mxtips + 1;

  		relabelInnerNodes(geneTree->start->back, geneTree, &innerNodeNumber, &innerBranches);
		rf = calculateRFDistance(tr, geneTree, innerBranches, adef, taxonToReduction, taxonHasDeg, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);//rf = 0.0;
		printf("target rf distance: %f, current rf distance: %f \n", tr->geneRFDistances[i], rf);
		fprintf(output, "%s", geneTree->tree_string);
		printf("%s: %d\n","Gene Tree", i);
		rewind(fp);	
		
		//freeMultifurcations(geneTree);
		freeTree(geneTree);
		rax_free(geneTree);
	}
	
	//free allocated memory
	rax_free(referenceSpeciesTreePath);
	rax_free(newSpeciesTreePath);
	rax_free(geneTreePath);
	rax_free(taxonHasDeg);
	rax_free(taxonToReduction);
	rax_free(taxonToEulerIndex);
	rax_free(taxonToLabel);
	rax_free(eulerIndexToLabel);
	rax_free(labelToTaxon);
	freeReferenceTree(tr);
	rax_free(tr);
	rax_free(adef);
	rax_free(rdta);
	rax_free(cdta);
	fclose(fp);
	fclose(output);	

	return 0;
}
