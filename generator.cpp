#include <unistd.h>
extern "C" {
#include "RAxML/axml.h"
}
#include "generateGeneTree.cpp"
#include "treePreparation.cpp"
#include <stdio.h>
#include <math.h>
#include <random>

//TODO set disclaimer to RAxML

int main(int argc, char* argv[]) {
	
	rawdata			*rdta;
	cruncheddata 	*cdta;
	tree			*tr;
	analdef			*adef;
	
	char 
		*referenceSpeciesTreePath = NULL,
		*geneTreePath = NULL;

	int 
		numLeaf1 = 0,
		numLeaf2 = 0,
		i = 0,
		loopIterations,
		newTaxaCount = 0,
		oldTaxaCount = 0;
		
	double 
		rf,
		rf2,
		ratio;
		
	FILE 
		*fp = myfopen(referenceSpeciesTreePath, "r"),
		*output = myfopen("../geneTrees.tre", "w");
		
	std::default_random_engine generator;
		
	copyArgs(&referenceSpeciesTreePath, argv[1]);
	copyArgs(&geneTreePath, argv[2]);
	newTaxaCount = argv[3];

	adef = (analdef *)rax_malloc(sizeof(analdef));
	rdta = (rawdata *)rax_malloc(sizeof(rawdata));
	cdta = (cruncheddata *)rax_malloc(sizeof(cruncheddata));
	tr	 = (tree *)rax_malloc(sizeof(tree));
	initAdef(adef);
	
	adef->restart = TRUE;
	adef->model = M_GTRCAT;
	adef->useInvariant = FALSE;
	adef->readTaxaOnly = TRUE;
	
	//read in reference tree
	//extractTaxaFromTopology(tr, rdta, cdta, referenceSpeciesTreePath); 
	oldTaxaCount = extractTaxaFromTopologyAndAddNew(tr, rdta, cdta, referenceSpeciesTreePath, newTaxaCount); 
	getinput(adef, rdta, cdta, tr);
	checkOutgroups(tr, adef);	
	treeReadLen(fp, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);
	rewind(fp);
	
	if(oldTaxaCount < newTaxaCount)	ratio = (double)newTaxaCount / (double)oldTaxaCount;
	else ratio = 1.0;
	
	//enlarge reference Tree
	enlargeTree(tr, oldTaxaCount &generator);

	//Allocate memory for the information of rf calculation of the reference tree
	int
		*taxonToLabel  = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),
		*eulerIndexToLabel = (int *)rax_malloc((4*tr->mxtips - 5) * sizeof(int)),
		*labelToTaxon = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),
		*taxonToEulerIndex  = (int *)rax_malloc((tr->mxtips) * sizeof(int)),
		*taxonToReduction = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
	
	//Get the necessary information of the reference tree	
	prepareReferenceTree(tr, taxonToReduction, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);
	//Calulcate taxa distribution and RF distances for all gene trees
	getGeneTreeStatistics(tr, geneTreePath, adef, taxonToReduction, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);
	
	/*
		Build new geneTrees based on the gathered statistics
	*/
	for(; i < tr->numberOfTrees; i++){
		tree *geneTree = (tree *)rax_malloc(sizeof(tree));
		int taxaGeneTree = (int)((double)tr->geneLeafDistributions[i] * ratio);
		//maps the taxa number of the gene tree to the taxa number of the reference tree
		int* treeTaxa = (int *)rax_malloc((tr->ntips) * sizeof(int));
		
		std::uniform_int_distribution<int> distribution(0,geneTree->ntips-1);

		//init HashTable with all TaxaNames because we don't know at the moment which taxa won't be in the gene tree
		geneTree->nameHash 	= tr->nameHash;		
		geneTree->rdta 		= tr->rdta;
		geneTree->nameList	= tr->nameList;		
		geneTree->rooted	= FALSE;
		loopIterations = 10;
		
		//Generate one geneTree
		generateGeneTree(tr, geneTree, taxaGeneTree, adef, ratio, treeTaxa, &generator);
		
		Tree2String(geneTree->tree_string, geneTree, geneTree->start->back, FALSE, TRUE,     FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
		printf("%s", geneTree->tree_string);

		rf = calculateRFDistance(tr, geneTree, geneTree->numberOfBranches, taxonToReduction, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);//rf = 0.0;
		printf("target rf distance: %f, current rf distance: %f \n", tr->geneRFDistances[i], rf);

		//switch leafs until the RF distance is satisfying
		while(!((rf >= tr->geneRFDistances[i]*0.9) && (rf <= tr->geneRFDistances[i]*1.1)) && loopIterations != 0)
		{
			numLeaf1 = distribution(generator);
			numLeaf2 = distribution(generator);
			//if the same number occured try again
			while(numLeaf1 == numLeaf2) numLeaf2 = distribution(generator);
			
			switchLeafs(geneTree, numLeaf1, numLeaf2);
			rf2 = calculateRFDistance(tr, geneTree, geneTree->numberOfBranches, taxonToReduction, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);//rf = 0.0;

			//if the switch is an iprovement accept the changes
			if(((rf2 > rf) && (rf2 <= tr->geneRFDistances[i]) && (rf < tr->geneRFDistances[i])) || ((rf2 < rf) && (rf2 >= tr->geneRFDistances[i]) && (rf > tr->geneRFDistances[i])))
			{
				rf = rf2;
			}
			//if the switch is no improvement accept the changes with a probability	
			else
			{
				//TODO: Accept step back with certain propability
				//switch leafs back to the step before
				switchLeafs(geneTree, numLeaf1, numLeaf2, treeTaxa, taxonToReduction);
				loopIterations--;
			}
		}
		printf("target rf distance: %f, current rf distance: %f \n", tr->geneRFDistances[i], rf);

		fprintf(output, "%s", geneTree->tree_string);
		printf("%s: %d\n","Gene Tree", i);
		//rewind(fp);	
		
		//freeMultifurcations(geneTree);
		rax_free(treeTaxa)
		freeTree(geneTree);
		rax_free(geneTree);
	}
	
	//free allocated memory
	rax_free(referenceSpeciesTreePath);
	rax_free(geneTreePath);
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
