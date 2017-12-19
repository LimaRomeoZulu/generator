#include "treePreparation.h"

int sortIntegers(const void *a, const void *b)
{
	int 
		ia = *(int *)(a),
	ib = *(int *)(b);

	if(ia == ib)
		return 0;

	if(ib > ia)
		return -1;
	else
		return 1;
}

void getNumberOfTrees(tree *tr, FILE *f)
{
	int 
		trees = 0,
	ch;

	while((ch = fgetc(f)) != EOF)
		if(ch == ';')
			trees++;

	assert(trees > 0);

	tr->numberOfTrees = trees;

	rewind(f);
}

void copyArgs(char **targetVariable, char *sourceArgv)
{
	
	size_t length = strlen(sourceArgv)+1;
	*targetVariable = (char *)rax_malloc(length * sizeof(char *));
	memcpy(*targetVariable, sourceArgv, length);
}

int extractTaxaFromTopologyAndAddNew(tree *tr, rawdata *rdta, cruncheddata *cdta, char fileName[1024], int newTaxaCount)
{
	FILE 
		*f = myfopen(fileName, "rb");

	char 
		**nameList,
	buffer[nmlngth + 2]; 

	int
		i = 0,
	c,
	taxaSize = 1024,
	taxaCount = 0,
	oldTaxaCount = 0;

	nameList = (char**)rax_malloc(sizeof(char*) * taxaSize);  

	while((c = fgetc(f)) != ';')
	{

		if(c == '(' || c == ',')
		{
			c = fgetc(f);
			if(c ==  '(' || c == ',')
			{
				ungetc(c, f);
			}
			else
			{	
				i = 0;	
	
				do
				{
					buffer[i++] = c;
					c = fgetc(f);
				}
				while(c != ':' && c != ')' && c != ',');
				buffer[i] = '\0';	    

				if(taxaCount == taxaSize)
				{		  
					taxaSize *= 2;
					nameList = (char **)rax_realloc(nameList, sizeof(char*) * taxaSize, FALSE);		 
				}

				nameList[taxaCount] = (char*)rax_malloc(sizeof(char) * (strlen(buffer) + 1));
				strcpy(nameList[taxaCount], buffer);
	
				taxaCount++;
			
				ungetc(c, f);
			}
		}
	}

	oldTaxaCount = taxaCount;
	
	while(newTaxaCount > taxaCount)
	{
		sprintf(buffer, "Taxon_%d", taxaCount);
		
		if(taxaCount == taxaSize)
		{		  
			taxaSize *= 2;
			nameList = (char **)rax_realloc(nameList, sizeof(char*) * taxaSize, FALSE);		 
		}

		nameList[taxaCount] = (char*)rax_malloc(sizeof(char) * (strlen(buffer) + 1));
		strcpy(nameList[taxaCount], buffer);

		taxaCount++;
	}

	/* BEGIN ensuring no taxon occurs twice */
	{
		char 
			**taxList = (char **)rax_malloc(sizeof(char *) * (size_t)taxaCount); 

		for(i = 0; i < taxaCount; ++i)
			taxList[i] = nameList[i]; 

		qsort(taxList, taxaCount, sizeof(char**), sortIntegers); 

		for(i = 1; i < taxaCount; ++i)
		{	
			if(strcmp(taxList[i], taxList[i-1]) == 0)
			{
				printf("\n\nA taxon labelled by %s appears twice in the first tree of tree collection, exiting ...\n\n", taxList[i]);
				exit(-1);
			}
		}
		rax_free(taxList);
	}
	/* END */


	printf("Found a total of %d taxa in first tree of tree collection \n", taxaCount);
	printf("Expecting all remaining trees in collection to have the same taxon set\n");

	rdta->numsp = taxaCount;

	tr->nameList = (char **)rax_malloc(sizeof(char *) * (taxaCount + 1));  
	for(i = 1; i <= taxaCount; i++)
		tr->nameList[i] = nameList[i - 1];
  
	rax_free(nameList);

	tr->rdta       = rdta;
	tr->cdta       = cdta;

	if (rdta->numsp < 4)
	{    
		printf("TOO FEW SPECIES, tree contains only %d species\n", rdta->numsp);
		assert(0);
	}

	tr->nameHash = initStringHashTable(10 * taxaCount);
	for(i = 1; i <= taxaCount; i++)   
	{
		printf("add [%s]\n", tr->nameList[i]);
		addword(tr->nameList[i], tr->nameHash, i);
	}

	fclose(f);
  
	return oldTaxaCount;
}

void prepareReferenceTree(tree *tr, int *taxonToReduction, int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon,  int *eulerIndexToLabel)
{
	int 
		newcount = 0; //counter used for correct traversals

	/* Preorder-Traversal of the large tree */
	preOrderTraversal(tr->start->back,tr->mxtips, tr->start->number, taxonToLabel, labelToTaxon, &newcount);

	newcount = 0; //counter set to 0 to be now used for Eulertraversal

	/* Init taxonToEulerIndex and taxonToReduction */
	int 
		ix;

	for(ix = 0; ix < tr->mxtips; ++ix)    
		taxonToEulerIndex[ix] = -1;    
  
	for(ix = 0; ix < (2*tr->mxtips - 2); ++ix)    
		taxonToReduction[ix] = -1;    


	/* Eulertraversal of the large tree*/
	unrootedEulerTour(tr->start->back,tr->mxtips, eulerIndexToLabel, taxonToLabel, &newcount, taxonToEulerIndex);

	/* Creating RMQ Datastructure for efficient retrieval of LCAs, using Johannes Fischers Library rewritten in C
	Following Files: rmq.h,rmqs.c,rmqs.h are included in Makefile.RMQ.gcc */
	RMQ_succinct(eulerIndexToLabel,4*tr->mxtips - 5);
	
}
void setupGeneTree(tree *geneTree, int taxaGeneTree){
	nodeptr	p0, p, q;
	int
		i,
	j,	
	tips,
	inter;
		
	tips = taxaGeneTree;
	inter = tips-1;
		
	geneTree->numberOfTrees = -1;
	geneTree->rooted = FALSE;
	//gives the maximum number of branches for one leaf
	geneTree->numBranches = 0;
	//gives the total number of inner branches
	geneTree->numberOfBranches = 0;	

	geneTree->treeStringLength = 
		2 * (size_t)tips + //parentheses
			2 * (size_t)tips * 64 + //branche lengths with : and . and branch labels and
				(size_t)tips + //commas
					1 + //closing semicolon 
						(size_t)tips * nmlngth; //taxon names

	geneTree->tree_string	= (char*)rax_calloc(geneTree->treeStringLength, sizeof(char)); 
		
	if (!(p0 = (nodeptr) rax_malloc((tips + 3*inter) * sizeof(node))))
	{
		printf("ERROR: Unable to obtain sufficient tree memory\n");
		return ;
	}

	if (!(geneTree->nodep = (nodeptr *) rax_malloc((2*tips) * sizeof(nodeptr))))
	{
		printf("ERROR: Unable to obtain sufficient tree memory, too\n");
		return ;
	}
	
	geneTree->nodep[0] = (node *) NULL;		/* Use as 1-based array */

	for (i = 1; i <= tips; i++)
	{
		p = p0++;

		p->hash	 =	KISS32(); /* hast table stuff */
		p->x			=	0;
		p->number =	i;
		p->next	 =	p;
		p->back	 = (node *)NULL;
		p->bInf	 = (branchInfo *)NULL;
		p->support = -1;

		geneTree->nodep[i] = p;
	}

	for (i = tips + 1; i <= tips + inter; i++)
	{
		q = (node *) NULL;
		for (j = 1; j <= 3; j++)
		{	 
			p = p0++;
			if(j == 1)
				p->x = 1;
			else
				p->x =	0;
				
			p->number = i;
			p->next	 = q;
			p->bInf	 = (branchInfo *)NULL;
			p->back	 = (node *) NULL;
			p->hash	 = 0;
			p->support = -2;

			q = p;
		}
		p->next->next->next 	= p;
		geneTree->nodep[i] 		= p;
	}

	geneTree->likelihood	= unlikely;
	geneTree->start			= (node *) NULL;
	geneTree->ntips			= 0;
	geneTree->nextnode		= 0;
}
void getTaxaDistribution(tree *tr, FILE  *input)
{
	char
		*word = NULL;
	word = (char*)rax_malloc(25*sizeof(char));
	char ch;


	unsigned long count;	
	for(int i = 1; i <= tr->mxtips; i++){
		count = 0;
		int read_counter;
		while((ch = treeGetCh(input)) == '('){}
		ungetc(ch,input);
		while(!feof(input))
		{
			read_counter = fscanf(input, "%[^,:();%f]%*[%f:.,();]", word);
			if(read_counter >= 0) 
			{
				if(strcmp(word, tr->nameList[i])==0)	count++;
			}
		}
		rewind(input);
		if(i == 1)
		{
			tr->taxaOccurencePrefixSum[i] = count;
		}
		else
		{
			tr->taxaOccurencePrefixSum[i] = tr->taxaOccurencePrefixSum[i-1] + count;
		}
		tr->nodep[i]->rec_distr = (float)count/tr->numberOfTrees;
	}	rax_free(word);
}


float calculateRFDistance(tree *tr, tree *geneTree, int numberOfSplits, int *taxonToReduction, int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon, int *eulerIndexToLabel, int* geneTreeTaxa)
{
	int *taxonHasDeg = (int *)rax_calloc((2*tr->mxtips - 2),sizeof(int));
	if(geneTree->ntips < 3)
	{
		return -1.0;
	} 
	int
		firstTaxon;           

	float
		rec_rf,
	maxRF;

	/* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	/* note that here we need to pay attention, since the RF distance is not normalized 
	by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	will potentially contain less bipartitions. 
	Hence the normalization factor is obtained as n-3 + numberOfSplits, where n-3 is the number 
	of bipartitions of the pruned down large reference tree for which we know that it is 
	bifurcating/strictly binary */
	  
	maxRF = (double)(2 * numberOfSplits);
	  
	/* now get the index of the first taxon of the small tree.
	we will use this to unambiguously store the bipartitions 
	*/
	  
	firstTaxon = geneTree->start->number;
	  
	/***********************************************************************************/
	/* Reconstruction Step */
	  
	/* Init hashtable to store Bipartitions of the induced subtree */
	/* 
	using geneTree->ntips instead of geneTree->mxtips yields faster code 
	e.g. 120 versus 128 seconds for 20,000 small trees on my laptop 
	*/
	hashtable
		*s_hash = initHashTable(geneTree->ntips * 4);
	  
	/* geneTreeTaxa[geneTree->ntips]; 
	Stores all taxa numbers from geneTree into an array called geneTreeTaxa: (Index) -> (Taxonnumber)  */
	  
	int 
		newcount = 0,
		newcount2 = 0;
	  
	/* seq2[2*geneTree->ntips - 2]; 
	stores PreorderSequence of the reference geneTree: (Preorderindex) -> (Taxonnumber) */
	int* 
		seq2 = (int *)rax_malloc((2*geneTree->ntips - 2) * sizeof(int));
	/* used to store the vectorLength of the bitvector */
	unsigned int 
		vectorLength;
	  
	/* extract all taxa of the geneTree and store it into an array, 
	also store all counts of taxa and nontaxa in taxonToReduction */
	rec_extractTaxa(geneTreeTaxa, taxonToReduction, geneTree->start, geneTree->mxtips, &newcount, &newcount2);
	  
	rec_extractTaxa(geneTreeTaxa, taxonToReduction, geneTree->start->back, geneTree->mxtips, &newcount, &newcount2);
	  
	/* counter is set to 0 to correctly preorder traverse the small tree */
	newcount = 0;
	  
	/* Preordertraversal of the small tree and save its sequence into seq2 for later extracting the bipartitions, it
	also stores information about the degree of every node */
	  
	rec_preOrderTraversalMulti(geneTree->start->back,geneTree->mxtips, geneTree->start->number, seq2, taxonHasDeg, &newcount);
	  
	/* calculate the bitvector length */
	if(geneTree->ntips % MASK_LENGTH == 0)
		vectorLength = geneTree->ntips / MASK_LENGTH;
	else
		vectorLength = 1 + (geneTree->ntips / MASK_LENGTH); 
	  
	unsigned int 
		**bitVectors = rec_initBitVector(geneTree, vectorLength);
	  
	/* store all non trivial bitvectors using an subtree approach for the induced subtree and 
	store it into a hashtable, this method was changed for multifurcation */
	rec_extractBipartitionsMulti(bitVectors, seq2, newcount,tr->mxtips, vectorLength, geneTree->ntips, 
	firstTaxon, s_hash, taxonToReduction, taxonHasDeg, numberOfSplits);
	  
	/* counter is set to 0 to be used for correctly storing all EulerIndices */
	newcount = 0; 
	  
	/* geneTreeTaxonToEulerIndex[geneTree->ntips]; 
	Saves all first Euler indices for all Taxons appearing in small Tree: 
	(Index) -> (Index of the Eulertour where the taxonnumber of the small tree first appears) */
	int* 
		geneTreeTaxonToEulerIndex = (int *)rax_malloc((geneTree->ntips) * sizeof(int));
	  
	/* seq[(geneTree->ntips*2) - 1] 
	Stores the Preordersequence of the induced small tree */
	int* 
		seq = (int *)rax_malloc((2*geneTree->ntips - 1) * sizeof(int));
	  
	int ix; 
	/* iterate through all small tree taxa */
	for(ix = 0; ix < geneTree->ntips; ix++) 
	{        
		int 
			taxanumber = geneTreeTaxa[ix];
	      
		/* To create geneTreeTaxonToEulerIndex we filter taxonToEulerIndex for taxa in the small tree*/
		geneTreeTaxonToEulerIndex[newcount] = taxonToEulerIndex[taxanumber-1]; 
	      
		/* Saves all Preorderlabel of the geneTree taxa in seq*/
		seq[newcount] = taxonToLabel[taxanumber-1];
	      
		newcount++;
	}
	  
	/* sort the euler indices to correctly calculate LCA */
	//quicksort(geneTreeTaxonToEulerIndex,0,newcount - 1);             
	  
	qsort(geneTreeTaxonToEulerIndex, newcount, sizeof(int), sortIntegers);
	  
	//printf("newcount2 %i \n", newcount2);      
	/* Iterate through all small tree taxa */
	for(ix = 1; ix < newcount; ix++)
	{  
		/* query LCAs using RMQ Datastructure */
		seq[newcount - 1 + ix] =  eulerIndexToLabel[query(geneTreeTaxonToEulerIndex[ix - 1],geneTreeTaxonToEulerIndex[ix])]; 	 
	      
		/* Used for dynamic programming. We save an index for every inner node:
		For example the reference tree has 3 inner nodes which we saves them as 0,1,2.
		Now we calculate for example 5 LCA to construct the induced subtree, which are also inner nodes. 
		Therefore we mark them as 3,4,5,6,7  */
	      
		taxonToReduction[labelToTaxon[seq[newcount - 1 + ix]] - 1] = newcount2;
	      
		newcount2 += 1;
	}
	  
	/* sort to construct the Preordersequence of the induced subtree */
	//quicksort(seq,0,(2*geneTree->ntips - 2));
	  
	qsort(seq, (2 * geneTree->ntips - 2) + 1, sizeof(int), sortIntegers);
	  
	/* calculates all bipartitions of the reference small tree and count how many bipartition it shares with the induced small tree */
	int 
		rec_bips = rec_findBipartitions(bitVectors, seq,(2*geneTree->ntips - 1), labelToTaxon, tr->mxtips, vectorLength, geneTree->ntips, firstTaxon, s_hash, taxonToReduction);
	  
	rec_rf = (double)(2 * (numberOfSplits - rec_bips)) / maxRF;
	  
	assert(numberOfSplits >= rec_bips);
	  
	rec_freeBitVector(geneTree, bitVectors);
	rax_free(bitVectors);
	  
	freeHashTable(s_hash);
	rax_free(s_hash);
	rax_free(taxonHasDeg);
	rax_free(seq);
	rax_free(seq2);
	rax_free(geneTreeTaxonToEulerIndex);

	return rec_rf;
}

void getGeneTreeStatistics(tree *tr, char *geneTreeFileName, analdef *adef, int *taxonToReduction, int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon, int *eulerIndexToLabel)
{
	FILE 
		*input = myfopen(geneTreeFileName, "r");
    
	tree 
		*smallTree = (tree *)rax_malloc(sizeof(tree));
	
	float
		rf_dist;
	
	/* now see how many small trees we have */
	getNumberOfTrees(tr, input);
	
	tr->geneLeafDistributions = (int *)rax_calloc(tr->numberOfTrees + 1, sizeof(int)); 
	tr->geneRFDistances = (float *)rax_calloc(tr->numberOfTrees + 1, sizeof(float)); 
	tr->taxaOccurencePrefixSum = (unsigned long *)rax_calloc(tr->mxtips + 1, sizeof(unsigned long)); 

	
	getTaxaDistribution(tr, input);
	allocateMultifurcations(tr, smallTree);
	
	
	for(int i = 0; i < tr->numberOfTrees;  i++)
	{
		int numberOfSplits = readMultifurcatingTree(input, smallTree, adef, FALSE); //LR set to false
		if(numberOfSplits > 0)
		{
			 int* smallTreeTaxa = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
			rf_dist = calculateRFDistance(tr, smallTree, numberOfSplits, taxonToReduction, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel, smallTreeTaxa);
			rax_free(smallTreeTaxa);
		}	
		tr->geneLeafDistributions[i] = smallTree->ntips;
		tr->geneRFDistances[i] = rf_dist;
	}
	freeMultifurcations(smallTree);
	rax_free(smallTree);
	fclose(input);
}

void enlargeTree(tree *tr, int oldTaxaCount, std::default_random_engine *generator)
{
	int
		i,
		leaf,
		lastInner;
	
	std::uniform_int_distribution<int> distribution(1,oldTaxaCount);
	
	//Get the next inner node
	lastInner = (tr->nextnode)++;
	
	if(tr->nodep[lastInner]->back != NULL) assert(0);
	
	
	for(i = oldTaxaCount + 1; i <= tr->ntips; i++)
	{
		leaf = distribution(*generator);
		
		hookupDefault(tr->nodep[leaf]->back, tr->nodep[lastInner], tr->numBranches);
		hookupDefault(tr->nodep[leaf], tr->nodep[lastInner]->next, tr->numBranches);
		hookupDefault(tr->nodep[i], tr->nodep[lastInner]->next->next, tr->numBranches);
		
		lastInner = (tr->nextnode)++;
		//TODO adjust distribution
	}

}

void freeTree(tree *tr)
{
	rax_free(tr->tree_string);
	rax_free(tr->nodep[1]);
	rax_free(tr->nodep);
	rax_free(tr->taxaOccurencePrefixSum);
}

void freeReferenceTree(tree *tr)
{
	//freeHashTable((hashtable*)tr->nameHash);
	//rax_free(tr->nameHash->table);
	rax_free(tr->nameHash);
	for(int i =1; i <= tr->ntips; i++)
	{
		rax_free(tr->nameList[i]);
	}
	rax_free(tr->nameList);
	rax_free(tr->geneLeafDistributions);
	rax_free(tr->geneRFDistances);
	freeTree(tr);
	
}
