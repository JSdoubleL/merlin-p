#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER

// Thing for running the python script (copied from Constrained-INC github)
#define ERROR -1
#define SUCCESS1 0 
static const char SCRIPT[] = "IKC.py ";
static const char LOG[] = "Script failed to run...\n";
#define PRINT_AND_RETURN(p, r)\
    do{\
      printf("%s\n", p);\
      return r;\
    } while(0) 
#define SYSCAL(err_val, log, format,...)\
    do{\
      sprintf(sys_buffer, format, __VA_ARGS__);\
      if(system(sys_buffer) != SUCCESS1) PRINT_AND_RETURN(log, err_val);\
    } while(0)

#include <map>
#include <string>
#include <list>
#include <vector>
#include <queue>
using namespace std;
class VariableManager;
class HierarchicalClusterNode;
class HierarchicalCluster
{
	public:
	HierarchicalCluster();
	~HierarchicalCluster();
	int setOutputDir(const char*);
	int cluster(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double);
	int cluster_ikc(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double);
	HierarchicalClusterNode* getRoot();
	int setVariableManager(VariableManager* p);
	private:
	int estimatePairwiseDist(map<int,HierarchicalClusterNode*>&);
	int mergePairs(map<int,HierarchicalClusterNode*>&);
	int mergePairs_LazyDelete(map<int,HierarchicalClusterNode*>&);
	int clearNeighborsFromList(HierarchicalClusterNode* node, map<string,int>& deleteMe);
	int generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>&);
	int generateModules_ikc(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>&);
	int calculatePercentVarianceExplained(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int calculateSilhouetteIndex(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int populateMembers(map<string,int>* members,HierarchicalClusterNode* node);
	int read_ikc_modules(string currMoudleFName,map<int,map<string,int>*>& modules);
	int findRoot(HierarchicalClusterNode*n);
	double threshold;
	map<string,HierarchicalClusterNode*> backup;
	char outputDir[1024];
	HierarchicalClusterNode* root;
	int globalNodeID;
	int globalPairID;
	double** distvalues;
	int** neighbors;
	int* visited;
	int treenodecnt;

	char currModuleFName[1024];
	char sys_buffer[10000];

	struct Pair
	{
		int node1;
		int node2;
		double value;
	};
	class ComparePair
	{
		public:
		bool operator() (Pair* p1, Pair* p2)
		{
			if(p1->value > p2->value)
			{
				return true;
			}
			return false;
		}
	};
	priority_queue<HierarchicalCluster::Pair*,vector<HierarchicalCluster::Pair*>, HierarchicalCluster::ComparePair> myqueue;
	VariableManager* vMgr;
	
};
#endif
