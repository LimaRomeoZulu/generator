#include "treePreparation.h"
#include "geneTreeGeneration.h"

/*
@param n 	returns -1 if an error occured, 
returns 0 if the taxa was not added to the tree due to statistical representation, 
returns the number of the taxa that can be added to then inner branch
*/
static int addElement (nodeptr ref, tree *tr, boolean readBranchLengths, int *lcount, analdef *adef, boolean storeBranchLabels, boolean* inTreeMapping)
{	 
	nodeptr	
		q,
		tmp;
	int			
		inner, 
		leftChild,
		rightChild;
		
	if (isTip(ref->number, tr->mxtips))
	{
		//add taxa to tree
		if(inTreeMapping[(ref->number)-1])	return ref->number;
		//don't add taxa. This results in one inter node less. 
		else return 0;
	} 
	else
	{ 	
		//add left child
		leftChild = addElement(ref->next->back, tr, readBranchLengths, lcount, adef, storeBranchLabels, inTreeMapping);
		
		//add right child
		rightChild = addElement(ref->next->next->back, tr, readBranchLengths, lcount, adef, storeBranchLabels, inTreeMapping);
		if (leftChild < 0 || rightChild < 0)
		{
			return -1;
		}
		else if(leftChild > 0 && rightChild == 0)
		{
			return leftChild;
		}
		else if(leftChild == 0 && rightChild > 0)
		{
			return rightChild;
		}
		else if(leftChild == 0 && rightChild == 0)
		{
			return 0;
		}
		else
		{
			inner = (tr->nextnode)++;
			q = tr->nodep[inner];
			tmp = tr->nodep[leftChild];				

			if(isTip(leftChild, tr->mxtips))
			{
				if (tr->start->number > leftChild)	tr->start = tmp;
			}
			hookupDefault(q->next, tmp, tr->numBranches);
			
			tmp = tr->nodep[rightChild];
			if(isTip(rightChild, tr->mxtips))
			{
				if (tr->start->number > rightChild)	tr->start = tmp;
			}
			hookupDefault(q->next->next, tmp, tr->numBranches);

			tr->numberOfBranches++;
	
			return inner;
		}
	}
} 

void determineLeafs(tree *geneTree, int *taxaGeneTree, boolean *inTreeMapping, int const numberOfTrees, std::default_random_engine *generator)
{
	int j = geneTree->mxtips; 
	unsigned long draw = 0;
	float drawFloat = 0.0;
	std::uniform_int_distribution<unsigned long> distribution(geneTree->taxaOccurencePrefixSum[1], geneTree->taxaOccurencePrefixSum[geneTree->mxtips]);

	if(*taxaGeneTree == 0)
	{
		return;
	}
	else
	{
		draw = distribution(*generator);
		while((geneTree->taxaOccurencePrefixSum[j] > draw) && (j > 1))
		{
			j--;
		}
		//if the leaf was already choosen go either up or down an take the next fitting neighbor. 
		if((inTreeMapping[j-1] == TRUE))
		{
			//get a number between 0 and 1
			std::uniform_real_distribution<float> distributionFloat(0.0,1.0);
			drawFloat = distributionFloat(*generator);
			if(drawFloat < 0.5)
			{
				while((inTreeMapping[j-1]) && (j > 1))	j--;
			}
			else
			{
				while((inTreeMapping[j-1]) && (j < geneTree->mxtips))	j++;
			}
		}
		//At this point you have either found a suitable leaf or you try again. 
		if((inTreeMapping[j-1] == FALSE))
		{
			inTreeMapping[j-1] = TRUE;
			if(j != geneTree->mxtips)
			{
				geneTree->taxaOccurencePrefixSum[j] = geneTree->taxaOccurencePrefixSum[j+1];
			}
			*taxaGeneTree = *taxaGeneTree - 1;
		}
		else
		{
			determineLeafs(geneTree, taxaGeneTree, inTreeMapping, numberOfTrees, generator);
			return;
		}
	}
	determineLeafs(geneTree, taxaGeneTree, inTreeMapping, numberOfTrees, generator);
	return;
}

void generateGeneTree(tree *referenceTree, tree *geneTree, int taxaGeneTree, analdef *adef, double ratio, boolean* inTreeMapping, std::default_random_engine *generator){

	// Create result tree.
	nodeptr	
		p,
		tmp;
	
	int 
		leftChild = 0,
		start = 0,
		lcount = 1,
		i; 
	boolean
		readBranches 		= FALSE,
		storeBranchLabels	= FALSE;
		
	setupGeneTree(geneTree, referenceTree->ntips);
	
	geneTree->start			= geneTree->nodep[referenceTree->ntips];
	geneTree->ntips			= taxaGeneTree;
	geneTree->mxtips 		= referenceTree->mxtips;
	geneTree->nextnode		= geneTree->mxtips + 1;	
	
	//Copy the prefixSum Array to alter it for leafs that were already added
	geneTree->taxaOccurencePrefixSum = longDup(referenceTree->taxaOccurencePrefixSum, referenceTree->mxtips);
	
	//determine the leafs that will be in the tree
	determineLeafs(geneTree, &taxaGeneTree, inTreeMapping, referenceTree->numberOfTrees, generator);
	
	//loop recursively through the reference Tree and add Taxa that were set before
	leftChild = addElement(referenceTree->start->back, geneTree, readBranches, &lcount, adef, storeBranchLabels, inTreeMapping);
	
	//check if reference->start would be in the tree 
	if(inTreeMapping[referenceTree->start->number - 1])
	{
		start = referenceTree->start->number;        	
	}		
	p = geneTree->nodep[leftChild];

	if(leftChild > 0 && start > 0 )
	{
		tmp = geneTree->nodep[start];
		geneTree->start = tmp;
		hookupDefault(p, tmp, geneTree->numBranches);
	}
	//if start is not included the next inner node must be pulled up to avoid a rooted tree
	else if(leftChild > 0 && start == 0)
	{
		if(!isTip(p->next->back->number, geneTree->mxtips))
		{
			tmp = p->next->back;
			hookupDefault(p->next,tmp->next->next->back, geneTree->numBranches);
		}
		else
		{
			tmp = p->next->next->back;
			hookupDefault(p->next->next,tmp->next->next->back, geneTree->numBranches);
		}
			hookupDefault(p,tmp->next->back, geneTree->numBranches);
			tmp->next->next->back = (nodeptr) NULL;
			tmp->next->back = (nodeptr) NULL;
			tmp->back = (nodeptr) NULL;
			geneTree->numberOfBranches--;	
	}
	//counting of inner nodes -1 to get the inner branches
	geneTree->numberOfBranches--;	
	
	return;
}

void switchLeafs(tree *tr, int numLeaf1, int numLeaf2, int* geneTreeTaxa)
{

	int 
		newcount = 0,
		newcount2 = 0;

	nodeptr
		leaf1 = (nodeptr) NULL,
		leaf2 = (nodeptr) NULL,
		tmp = (nodeptr) NULL;

	leaf1 = tr->nodep[geneTreeTaxa[numLeaf1-1]];
	leaf2 = tr->nodep[geneTreeTaxa[numLeaf2-1]];

	tmp = leaf2->back;
	hookupDefault(leaf1->back, leaf2, tr->numBranches);
	hookupDefault(tmp, leaf1, tr->numBranches);

}

unsigned long * longDup(unsigned long const *src, size_t len)
{
   unsigned long *p = (unsigned long *) rax_malloc((len+1) * sizeof(unsigned long));
   memcpy(p, src, (len+1) * sizeof(unsigned long));
   return p;
}



