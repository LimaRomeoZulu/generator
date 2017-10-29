
#include <unistd.h>
#include "RAxML/axml.h"
#include "functions.c"
#include <stdio.h>
#include <math.h>


void generateGeneTree(FILE *fp, tree *geneTree, int taxaNewSpeciesTree, int taxaGeneTree, analdef *adef){
	// Create result tree.
	nodeptr	
		p;
	
	int 
		i, 
		ch, 
		addedTaxa,
		lcount = 0; 
	boolean
		readBranches 		= FALSE,
		readNodeLabels 		= TRUE,
		topologyOnly		= TRUE,
		completeTree		= TRUE,
		storeBranchLabels	= FALSE;
		
	setupGeneTree(geneTree, taxaGeneTree);
	
	geneTree->start			= geneTree->nodep[taxaGeneTree];
	geneTree->ntips			= 0;
	geneTree->nextnode		= geneTree->mxtips + 1;	
	
	p = geneTree->nodep[(geneTree->nextnode)++]; 
	
	float prob = float)taxaGeneTree/(float)taxaNewSpeciesTree;
	
	//loop to the first occurence of '(' which is the first subtree/inner node
	while((ch = treeGetCh(fp)) != '(')
	{
		if(ch == EOF)
		{
			printf("RAxML could not find a single \"(\" in what is supposed to be your tree file\n");
			errorExit(-1);
		}						
	};
	//add the left child
	addedTaxa = addElementLen(fp, geneTree, p, readBranches, readNodeLabels, &lcount, adef, storeBranchLabels, prob);
	if (n < 0) 							assert(0);
	if (! treeNeedCh(fp, ',', "in")) 	assert(0);
		
	//add the right child
	addedTaxa = addElementLen(fp, geneTree, p->next, readBranches, readNodeLabels, &lcount, adef, storeBranchLabels, prob);
	if (n < 0) 							assert(0);
	if (! geneTree->rooted) 
	{
		if ((ch = treeGetCh(fp)) == ',') 
		{ 
			//add the parent
			addedTaxa = addElementLen(fp, geneTree, p->next->next, readBranches, readNodeLabels, &lcount, adef, storeBranchLabels, prob);
			if (n < 0 )					assert(0);			
		}
		else 
		{			
			/*	A rooted format */
			geneTree->rooted = TRUE;
			geneTree->wasRooted		 = TRUE;
		
			if (ch != EOF)	(void) ungetc(ch, fp);
		}	
	}
	else 
	{						
		p->next->next->back = (nodeptr) NULL;
		geneTree->wasRooted		 = TRUE;		
	}
	
	
	Tree2String(geneTree->tree_string, geneTree, geneTree->start->back, FALSE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
	printf("%s", geneTree->tree_string);

	// How deep is the current token nested in the tree?
	//int depth = 0;
	//std::vector<int> depthCounter;
	/*auto childrenNotVisited	 = std::vector<size_t>( referenceSpeciesTree.node_count(), 0 );
	auto childrenInserted	 = std::vector<size_t>( referenceSpeciesTree.node_count(), 0 );
	float p = ((float)taxaGeneTree * ((float)taxaNewSpeciesTree/(float)taxareferenceSpeciesTree))/(float)taxaNewSpeciesTree;
	
	std::cout << "p: " << p << std::endl; 
	
	for (auto it : eulertour(referenceSpeciesTree)) {
	if (it.node().is_leaf()) {
	depth++;
			
	//float draw = 0.5;
	float draw = (rand() / (RAND_MAX + 1.0));
	//leaf is added to the tree
	//std::cout << draw << std::endl;
	if(draw <= p){
	node.name = it.node().data_cast<DefaultNodeData>()->name;
	int idx = it.node().index();
	//node.name = std::to_string(idx);
	node.depth = depth;

	broker.push_top( node );
	node = NewickBrokerElement();
	--childrenNotVisited[it.link().outer().node().index()];
	++childrenInserted[it.link().outer().node().index()];
	}
	//leaf is not added but was considered
	else{
	--childrenNotVisited[it.link().outer().node().index()];
	}
	}
	else{
	if(std::find(depthCounter.begin(), depthCounter.end(), it.node().index()) != depthCounter.end()) {
	depth--;
	if((childrenNotVisited[it.node().index()] == 0)&&(childrenInserted[it.node().index()]>1)){
	node.name = "a";
	node.depth = depth;
	broker.push_top( node );
	node = NewickBrokerElement();
	--childrenNotVisited[it.link().outer().node().index()];
	++childrenInserted[it.link().outer().node().index()];
	}
	if((childrenNotVisited[it.node().index()] == 0)&&(childrenInserted[it.node().index()]<1)){
	node = broker.top();
	broker.pop_top();
	node.depth = depth;
	broker.push_top( node );
	node = NewickBrokerElement();
	++childrenInserted[it.link().outer().node().index()];
	}

	} else {
	depthCounter.push_back(it.node().index());
	depth++;
	childrenNotVisited[it.node().index()] = it.node().rank();
	childrenInserted[it.node().index()] = 0;
	}
	}
	}

	node.name = "root";
	node.depth = 0;
	broker.push_top( node );
	node = NewickBrokerElement();
	
	return broker;*/
	return;
}

int main(int argc, char* argv[]) {
	
	rawdata			*rdta;
	cruncheddata *cdta;
	tree				 *tr;
	analdef			*adef;
	char modelChar = 'R';
	
	char *referenceSpeciesTreePath = NULL,
	*newSpeciesTreePath = NULL,
	*histogramPath = NULL,
	*RFMetrixPath = NULL;
	//std::ofstream evaluationTrees;
	
	copyArgs(&referenceSpeciesTreePath, argv[1]);
	copyArgs(&newSpeciesTreePath, argv[2]);
	copyArgs(&histogramPath, argv[3]);
	copyArgs(&RFMetrixPath, argv[4]);
	
	int *histogram = NULL; 
	float *rfMetrix = NULL; 
	
	histogram = readHistogram(histogramPath);
	//std::ifstream file_rfMetrix(RFMetrixPath);
	
	//TODO: check where it is freed
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
	extractTaxaFromTopology(tr, rdta, cdta, newSpeciesTreePath); 
	getinput(adef, rdta, cdta, tr);
	checkOutgroups(tr, adef);	

	int taxaNewSpeciesTree = rdta->numsp;
	int taxaReferenceSpeciesTree = countTaxaInTopology(referenceSpeciesTreePath);
	
	FILE 
		*fp = myfopen(newSpeciesTreePath, "r");
	
	//treeReadLen(fp, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);
	
	for(int i = 0 ; i < 1; i++){
		tree *geneTree = (tree *)rax_malloc(sizeof(tree));;
		//int taxaGeneTree = (int)roundf((float)histogram[i] * (float)taxaReferenceSpeciesTree/(float)taxaNewSpeciesTree) ;
		int taxaGeneTree = taxaNewSpeciesTree;
		//init HashTable with all TaxaNames because we don't know at the moment which taxa won't be in the gene tree
		geneTree->nameHash	= initStringHashTable(10*tr->mxtips);
		geneTree->nameHash 	= tr->nameHash;		
		geneTree->rdta 		= tr->rdta;
		geneTree->nameList	= tr->nameList;		

		generateGeneTree(fp, geneTree, taxaNewSpeciesTree, taxaGeneTree, adef);

		
		printf("%s: %d\n","Gene Tree", i);
		//std::cout << PrinterCompact().print( geneTree ) << std::endl;
		
		//auto writer = DefaultTreeNewickWriter();
		//writer.to_file(geneTree, "gene_Tree_" + std::to_string(i) + ".tre");
	
		rax_free(geneTree);
	}
	
	/*utils::InputStream instream(utils::make_unique<utils::FileInputSource>(referenceSpeciesTreePath));
	auto itTree = NewickInputIterator(instream, DefaultTreeNewickReader());

	
	while (itTree) { // iterate over the set of evaluation trees
	Tree const& tree = *itTree;
	for (auto it : eulertour(tree)) {
	if (it.node().is_leaf()) {
	++leafs;
	}
	else{
	++vertices;
	}
	}
	counting << tree_count << "," << vertices <<"," << leafs << std::endl;
	vertices = 0;
	leafs = 0;
	++tree_count;
	++itTree;
	}*/
	
	//free allocated memory
	rax_free(referenceSpeciesTreePath);
	rax_free(newSpeciesTreePath);
	rax_free(histogramPath);
	rax_free(RFMetrixPath);
	rax_free(tr);
	rax_free(adef);
	rax_free(rdta);
	rax_free(cdta);
	
	return 0;
}
