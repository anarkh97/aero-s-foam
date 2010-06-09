// ----------------------------------------------------------------
// HB - 08/10/03
// HB - Modified 02/28/04
// ----------------------------------------------------------------
#ifndef _NODALMORTARSHAPEFCT_H_
#define _NODALMORTARSHAPEFCT_H_

// STL
#include <vector>
#include <utility>
#include <map>

class Connectivity;
class FaceElement;
class FaceElemSet;
class FFIPolygon;
class LMPCons;

class NodalMortarShapeFct {
  private:
        typedef std::pair<int, double> node_pair_t;
        std::vector<node_pair_t > NodalData;
	//int nFaceElems;  
	//int* FaceElemsId; 
	//FaceElement** PtrFaceElems;

	std::vector<int> LinkedSlaveNodes;
	std::vector<int> LinkedMasterNodes;

	std::vector<double> SlaveMPCCoeffs;
	std::vector<double> MasterMPCCoeffs;
        double MPCRhs;
	
  public:
	// Constructors
	// ~~~~~~~~~~~~
	NodalMortarShapeFct();
	NodalMortarShapeFct(int ref_node, double ref_coeff);
	
	// Destructors
	// ~~~~~~~~~~~
	~NodalMortarShapeFct();

        // Initialization & clean/clear methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void ClearData();
        void Reset();

	// Set methods
	// ~~~~~~~~~~~
        void SetRefData(int ref_node, double ref_coeff);

	// Get methods
	// ~~~~~~~~~~~
        int GetnNodes()  { return NodalData.size(); }
        int GetRefNode() { return NodalData[0].first; }
        double GetRefNodalCoeff() { return NodalData[0].second; }

        int GetNodeId(int i)        { return NodalData[i].first; }
        double GetNodalCoeff(int i) { return NodalData[i].second; }

 	// Mortar LMPC methods
	// ~~~~~~~~~~~~~~~~~~~
	void MakeSlaveLink(Connectivity*, FaceElemSet*, bool Dual=false);
	void MakeMasterLink(Connectivity*, FaceElemSet*, Connectivity*, FFIPolygon*);
	void BuildMortarLMPC(Connectivity*, FaceElemSet*, Connectivity*, FFIPolygon*);
        void BuildMortarCtcLMPC(Connectivity*, FaceElemSet*, Connectivity*, FFIPolygon*);

        LMPCons* CreateMortarLMPCons(int, int, double, int* SlaveLlToGlNodeMap=0, 
                                                       int* MasterLlToGlNodeMap=0);

        LMPCons* CreateMortarCtcLMPCons(int lmpcnum, double rhs, int* SlaveLlToGlNodeMap=0,
                                                                 int* MasterLlToGlNodeMap=0);

	// Print, display methods
	// ~~~~~~~~~~~~~~~~~~~~~~
	void print(int* SlaveLlToGlNodeMap=0, int* MasterLlToGlNodeMap=0);

        // Test wet FSI methods (NEED TO BE REDESIGNED)
        // ~~~~~~~~~~~~~~~~~~~~
	void BuildWetFSICoupling(Connectivity*, FaceElemSet*, Connectivity*, FFIPolygon*);
        LMPCons* CreateWetFSICons(int* SlaveLlToGlNodeMap=0, int* MasterLlToGlNodeMap=0);
};

#endif

