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
#include <random>

extern "C" {
#include "RAxML/globalVariables.h"
#include "RAxML/rmq.h"
}

static int sortIntegers(const void *a, const void *b)
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

void getNumberOfTrees(tree *tr, FILE *f, analdef *adef)
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
	geneTree->numBranches = 0;

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
	geneTree->mxtips		= taxaGeneTree;
	geneTree->ntips			= 0;
	geneTree->nextnode		= 0;
}

/*
	@param n 	returns -1 if an error occured, 
				returns 0 if the taxa was not added to the tree due to statistical representation, 
				returns the number of the taxa that can be added to then inner branch
*/
static int addElementLen (FILE *fp, tree *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount, analdef *adef, boolean storeBranchLabels, float prob)
{	 
	nodeptr	
		q,
		tmp;
	int			
		n = -1,
		inner, 
		ch, 
		fres, 
		leftChild,
		rightChild;
		
	if ((ch = treeGetCh(fp)) == '(') 
	{ 
		inner = (tr->nextnode)++;
		if (inner > 2*(tr->mxtips) - 2) 
		{
			if (tr->rooted || inner > 2*(tr->mxtips) - 1) 
			{
				printf("ERROR: Too many internal nodes.	Is tree rooted?\n");
				printf("			 Deepest splitting should be a trifurcation.\n");
				return -1;
			}
			else 
			{
				if(readNodeLabels)
				{
					printf("The program will exit with an error in the next source code line\n");
					printf("You are probably trying to read in rooted trees with a RAxML option \n");
					printf("that for some reason expects unrooted binary trees\n\n");
				}

				assert(!readNodeLabels);
				tr->rooted = TRUE;
			}
		}
			
		q = tr->nodep[inner];
		
		//add left child
		leftChild = addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount, adef, storeBranchLabels, prob);
		if (! treeNeedCh(fp, ',', "in"))				return -1;
		
		//add right child
		rightChild = addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount, adef, storeBranchLabels, prob);
		if (leftChild < 0 || rightChild < 0)
		{
			return -1;
		}
		else if(leftChild > 0 && rightChild == 0)
		{
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			return leftChild;
		}
		else if(leftChild == 0 && rightChild > 0)
		{
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			return rightChild;
		}
		else if(leftChild == 0 && rightChild == 0)
		{
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			return 0;
		}

		else
		{
			tmp = tr->nodep[leftChild];
			if (tr->start->number > leftChild)	tr->start = tmp;
			hookupDefault(q->next, tmp, tr->numBranches);
			
			tmp = tr->nodep[rightChild];
			if (tr->start->number > rightChild)	tr->start = tmp;
			hookupDefault(q->next->next, tmp, tr->numBranches);
			
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			
			return inner;
		}
			
		if(readNodeLabels)
		{
			char label[64];
			int support;

			if(treeGetLabel (fp, label, 10, FALSE))
			{	
				int val = sscanf(label, "%d", &support);
			
				assert(val == 1);

				/*printf("LABEL %s Number %d\n", label, support);*/
				p->support = q->support = support;
				/*printf("%d %d %d %d\n", p->support, q->support, p->number, q->number);*/
				assert(p->number > tr->mxtips && q->number > tr->mxtips);
				*lcount = *lcount + 1;
			}
		}
		else	
			(void) treeFlushLabel(fp);
	}
	else 
	{	 
		//float draw = 0.5;
		float draw = (rand() / ((double)RAND_MAX + 1.0));
		//add taxa to tree

		ungetc(ch, fp);		
		if ((n = treeFindTipName(fp, tr, TRUE)) <= 0)					return -1;	
		fres = treeFlushLen(fp, tr);
		if(!fres) return -1;
		if(draw <= prob)
		{
			(tr->ntips)++;
			return n;
		}
		//don't add taxa. This results in one inter node less. 
		else
		{
			return 0;
		}
	}
					
} 


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
			if(c ==	'(' || c == ',')
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
			count++;
		}
		else{
			arraySize *= 2;
			result = (int *) realloc(result, arraySize * sizeof(int *));
		}
	} 
	return result;
}
double* readRFDistance(char fileName[1024])
{
	FILE 
		*f = myfopen(fileName, "rb");
	
	double 
		item = 0;
	int
	count = 0,
	arraySize = 1024;
	
	double *result;
	
	result = (double *) malloc(arraySize * sizeof(double *));
	
	while(fscanf(f,"%lf",&item) == 1)	
	{	
		if(count<arraySize){
			result[count] = item;
			count++;
		}
		else{
			arraySize *= 2;
			result = (double *) realloc(result, arraySize * sizeof(double *));
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
void prepareRefernceTree(tree *tr, int *taxonToReduction, int *taxonHasDeg    , int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon,  int *eulerIndexToLabel)
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
double calculateRFDistance(tree *tr, tree *geneTree, int numberOfSplits, analdef *adef, int *taxonToReduction, int *taxonHasDeg, int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon, int *eulerIndexToLabel)
{
	if(geneTree->ntips < 3)
	{
		return -1.0;
	} 

	  
      int
	  innerNodeNumber = geneTree->mxtips + 1;

	  int
	    firstTaxon;           

	  double
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
	  int* 
	    geneTreeTaxa = (int *)rax_malloc((geneTree->ntips) * sizeof(int));
	  
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
	  
	  rax_free(geneTreeTaxa);
	  rax_free(seq);
	  rax_free(seq2);
	  rax_free(geneTreeTaxonToEulerIndex);

	return rec_rf;
}

void getTaxaDistribution(tree *tr, FILE  *input, analdef *adef)
{
	char
		*word = (char*)malloc(25*sizeof(char));
	char ch;

	/* now see how many small trees we have */
	getNumberOfTrees(tr, input, adef);
	unsigned long count;	
	for(int i = 1; i <= tr->mxtips; i++){
		count = 0;
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
		tr->taxaOccurencePrefixSum[i] = tr->taxaOccurencePrefixSum[i-1] + count;
		tr->nodep[i]->rec_distr = (float)count/tr->numberOfTrees;
	}
	free(word);
}

void getGeneTreeStatistics(tree *tr, char *geneTreeFileName, analdef *adef, int *taxonToReduction, int *taxonHasDeg    , int *taxonToEulerIndex, int *taxonToLabel, int *labelToTaxon, int *eulerIndexToLabel)
{
    FILE 
	*input = myfopen(geneTreeFileName, "r");
    
	tree 
      *smallTree = (tree *)rax_malloc(sizeof(tree));
	
	float
		rf_dist;
	
	getTaxaDistribution(tr, input, adef);
	tr->geneLeafDistributions = (int *)rax_malloc(sizeof(int) * (tr->numberOfTrees + 1)); 
	tr->geneRFDistances = (float *)rax_malloc(sizeof(float) * (tr->numberOfTrees + 1)); 

	allocateMultifurcations(tr, smallTree);
	
	
	for(int i = 0; i < tr->numberOfTrees;  i++)
	{
        int numberOfSplits = readMultifurcatingTree(input, smallTree, adef, FALSE); //LR set to false
	if(numberOfSplits > 0)
	{
		rf_dist = calculateRFDistance(tr, smallTree, numberOfSplits, adef, taxonToReduction, taxonHasDeg, taxonToEulerIndex, taxonToLabel, labelToTaxon, eulerIndexToLabel);
	}	
	tr->geneLeafDistributions[i] = smallTree->ntips;
		tr->geneRFDistances[i] = rf_dist;
	}
  	freeMultifurcations(smallTree);
  	rax_free(smallTree);
	fclose(input);
}


void freeTree(tree *tr)
{
	rax_free(tr->tree_string);
	rax_free(tr->nodep[1]);
	rax_free(tr->nodep);
}

void freeReferenceTree(tree *tr)
{
	freeHashTable((hashtable*)tr->nameHash);
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


