#ifndef __GENE_TREE_GENERATION_INCLUDED__
#define __GENE_TREE_GENERATION_INCLUDED__   

#include <random>

static int addElement (nodeptr ref, tree *tr, boolean readBranchLengths, int *lcount, analdef *adef, boolean storeBranchLabels, int *treeTaxa, boolean* inTreeMapping);
void determineLeafs(tree *geneTree, int *taxaGeneTree, int const numberOfTrees, std::default_random_engine *generator);
void generateGeneTree(tree *referenceTree, tree *geneTree, int taxaGeneTree, analdef *adef, double ratio, int *treeTaxa, boolean* inTreeMapping, std::default_random_engine *generator);
void switchLeafs(tree *tr, int numLeaf1, int numLeaf2);
unsigned long * longDup(unsigned long const *src, size_t len);

#endif //__GENE_TREE_GENERATION_INCLUDED__
