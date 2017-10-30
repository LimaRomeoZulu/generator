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
			(tr->nextnode)--;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			return leftChild;
		}
		else if(leftChild == 0 && rightChild > 0)
		{
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			(tr->nextnode)--;
			return rightChild;
		}
		else if(leftChild == 0 && rightChild == 0)
		{
			if (! treeNeedCh(fp, ')', "in"))				return -1;
			(void) treeFlushLabel(fp);
			treeFlushLen(fp, tr);
			(tr->nextnode)--;
			return 0;
		}

		else
		{
			tmp = tr->nodep[leftChild];
			if (tr->start->number > leftChild)	tr->start = tmp;
			(tr->ntips)++;
			hookupDefault(q->next, tmp, tr->numBranches);
			
			tmp = tr->nodep[rightChild];
			if (tr->start->number > rightChild)	tr->start = tmp;
			(tr->ntips)++;
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

void copyArgs(char **targetVariable, char *sourceArgv)
{
	
	size_t length = strlen(sourceArgv)+1;
	*targetVariable = (char *)rax_malloc(length * sizeof(char *));
	memcpy(*targetVariable, sourceArgv, length);
}






