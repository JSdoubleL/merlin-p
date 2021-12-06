#ifndef _META_LEARNER_
#define _META_LEARNER_

// Thing for running the python script (copied from Constrained-INC github)
#define ERROR -1
#define SUCCESS1 0 
static const char SCRIPT[] = "testscript.py ";
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
#include "CommonTypes.H"
using namespace std;


//This is gamma map for each variable within each datapoint
typedef map<int,INTDBLMAP*> DATA_GAMMA_MAP;
typedef map<int,INTDBLMAP*>::iterator DATA_GAMMA_MAP_ITER;

class FactorManager;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class Graph;
class MetaLearner
{
	public:
	MetaLearner();
	~MetaLearner();
	int setInputFName(const char*);
	int setMaxFactorSize(int);
	int setMaxFactorSize_Approx(int);
	int setPenalty(double);
	int setBeta1(double aval);
	int setBeta_Motif(double bval);
	int setLambda(double bval);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setPreRandomizeSplit();
	int setTargetList(const char*);
	int setNoPrune(bool);
	int init();
	int initPartitions(int);
	int initEdgePriorMeta();
	int initEdgePriorMeta(map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePrior);
	int doCrossValidation(int);
	int start(int);
	int start_gradualMBIncrease(int);
	int start_gradualMBIncrease_RankRegulators(int);
	int getTopRegs(map<int,int>&);
	int setGlobalEvidenceManager(EvidenceManager*);
	int setHoldOutEvManager(EvidenceManager*);
	int setVariableManager(VariableManager*);
	int setOutputDirName(const char*);
	int setClusteringThreshold(double);
	int setSpecificFold(int);
	int showModelParameters(int);
	int getPredictionError_Holdout(int);
	int getPredictionError_CrossValid(int);
	int setRandom(bool);
	int readModuleMembership(const char*);
	int setDefaultModuleMembership();

	//int setMotifGraph(const char*);
	//int setChIPGraph(const char*);
	//int setBeta_ChIP(double bval);
	int setPriorGraph_All(const char*);
	private:
	
	//map<int,INTDBLMAP*> edgePriors_Motif;
	//map<int,INTDBLMAP*> edgePriors_ChIP;
	//map<string,map<string,double>*> priorGraph_Motif;
	//map<string,map<string,double>*> priorGraph_ChIP;
	//double beta_chip;
	//int initEdgePriorMeta_Motif();
	//int initEdgePriorMeta_ChIP();
	int initEdgePriorMeta_All();
	map <string,double> betamap;
	map <string,map<int,INTDBLMAP*>*> edgepriormap;
	map <string,map<string,map<string,double>*>*> priorgraphmap;
		
	double getInitPLLScore();
	double getPLLScore();
	double getInitPrior();
	double getPriorChange();
	double getEdgePrior(int tfID, int targetID);
	double getPLLScore_Datapoint(int datapoint);
	int clearFoldSpecData();
	int initEdgeSet(bool);
	int initVariableComponents();
	int initCondsetMap_Nopool();
	INTINTMAP* getConditionSet(int cind);

	int collectMoves(int);
	int collectMoves(int,int);
	int sortMoves();
	int makeMoves();
	int attemptMove(MetaMove*,map<int,INTINTMAP*>&);
	int getNewPLLScore(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement,Potential**);
	double getNewPLLScore_Condition(int cId, int csetId, int vId, Potential* newPot,int);
	double getNewPLLScore_Condition(int csetId, int vId, Potential* newPot);
	double getNewPLLScore_Condition(int csetId, int vId, int uID, Potential* newPot);
	double getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uID, Potential* newPot);
	double getNewPLLScore_Condition_Noprior(int csetId, int vId, Potential* newPot);
	double getMBPrior(int targetID);

	int genCondSetKey(INTINTMAP& condSet, string& aKey);

	int cleanUp();
	int dumpAllGraphs(int,int,int);
	int getEdgeVars(string&, char*, char*);
	bool checkMBSize(INTINTMAP*,int,int,int);
	bool checkMBSize(int,int,int,int);
	int populateGraphsFromFile();
	
	int normalizeWeight(INTDBLMAP&);
	//Update the potentials to reflect the new covariance and mean estimates and also recompute the likelihood.
	int updatePotentials();

	int setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph);
	int redefineModules();
	int redefineModules_Global();
	double getModuleContrib(string& tgtName, string& tfName);
	double getModuleContribLogistic(string& tgtName, string& tfName);
	int initPhysicalDegree();
	int getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet);
	double getModuleWideScoreImprovement(int cid,INTINTMAP& conditionSet,Variable* u, Variable* v,map<string,int>* moduleMembers);
	
	map<int,FactorManager*> fgMgrSet;
	map<int,PotentialManager*> potMgrSet;
	map<int,EvidenceManager*> evMgrSet;
	map<int,FactorGraph*> fgGraphSet;

	EvidenceManager* globalEvMgr;
	EvidenceManager* holdoutEvMgr;

	map<int,string> outLocMap;

	VariableManager* varManager;
	int maxFactorSize;
	int maxFactorSizeApprox;
	int maxMBCntGreedy;
	char inputFName[1024];
	double penalty;
	double convThreshold;
	char restrictedFName[1024];
	map<string,INTINTMAP*> edgeConditionMap;
	FactorGraph* fg;
	vector<MetaMove*> moveSet;
	map<string,double> edgeValue;
	map<int,INTINTMAP*> condsetMap;
	map<string,int> condsetKeyIDMap;
	map<int,string> condsetIDKeyMap;
	map<string,int> restrictedVarList;
	map<int,INTDBLMAP*> currPLLMap;
	char trueGraphFName[1024];
	map<int,double> finalScores;
	map<int,INTINTMAP*> affectedVariables;
	bool preRandomizeSplit;
	gsl_rng* rnd;
	char outputDirName[1024];
	map<int,INTINTMAP*> datasetInd;
	map<string,double> edgePresenceProb;
	map<int,double> varNeighborhoodPrior;
	//These are the edges that got added
	map<string,int> edgeUpdates;
	double beta1;
	double beta_motif;
	double lambda;
	map<string,int> reqdTargetList;
	bool random;
	bool randomData;
	bool noprune;
	map<int,map<string,int>*> moduleGeneSet;
	map<string,int> geneModuleID;
	map<int,map<string,int>*> moduleIndegree;
	map<string,int> regulatorModuleOutdegree;
	map<int,int> idVidMap;
	int globalID;
	double clusterThreshold;
	map<string,int> variableStatus;
	int currFold;
	char currModuleFName[1024];
	char foldoutDirName[1024];
	char inGraphFName[1024];
	int specificFold;
	char sys_buffer[10000];
};
#endif
