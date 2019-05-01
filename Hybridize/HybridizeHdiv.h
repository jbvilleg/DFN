#ifndef HybridizeHDiv_h
#define HybridizeHDiv_h

#include <stdio.h>
#include <map>
#include "pzcompel.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;

struct HybridizeHDiv {
    int fHDivWrapMatid = -10;
    int fLagrangeInterface = -9;
    int fInterfaceMatid = -8;
    int fNState = 1;
    
    HybridizeHDiv() = default;
    
    HybridizeHDiv(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    void ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid);
    
    /// compute material ids for the periferal material objects
    void ComputePeriferalMaterialIds(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// return true if a material id is a peripheral material
    bool IsPeriferalMaterialId(int matid)
    {
        return matid == fHDivWrapMatid || matid == fLagrangeInterface || matid == fInterfaceMatid;
    }
    /// split the connects between flux elements and create a dim-1 pressure element
    void HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// Create interface elements with material id InterfaceMatid
    void CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// create a multiphysics mesh for the hybrid formulation using the materials of another mesh and the given atomic meshes
    TPZCompMesh * CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// group and condense the elements
    static void GroupElements(TPZCompMesh *cmesh_Hybrid);
    
    /// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes
    void InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh
    void InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier = 1.);
    
    /// clones the atomic meshes in meshvec_HDiv and creates a multiphysics hybrid mesh
    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > Hybridize(TPZCompMesh *cmesh_Multiphysics, TPZVec<TPZCompMesh *> &meshvec_HDiv, bool group_elements=true, double Lagrange_term_multiplier = 1.);
    
    /// verify the consistency of the solution of the flux mesh
    static void VerifySolutionConsistency(TPZCompMesh *fluxmesh, std::ostream &out);
private:
    
    std::tuple<int64_t,int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
public:
    
    static TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side);
    
};

#endif
