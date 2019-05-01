
#include "HybridizeHdiv.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzintel.h"
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzgeoelside.h"
#include "pztrnsform.h"
#include "tpzintpoints.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZNullMaterial.h"

    
HybridizeHDiv::HybridizeHDiv(TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
}
    
void HybridizeHDiv::ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid){
    
}
    
    /// compute material ids for the periferal material objects
void HybridizeHDiv::ComputePeriferalMaterialIds(TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
}
    

    /// split the connects between flux elements and create a dim-1 pressure element
void HybridizeHDiv::HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
}
    
    /// Create interface elements with material id InterfaceMatid
void HybridizeHDiv::CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
}
    
    /// create a multiphysics mesh for the hybrid formulation using the materials of another mesh and the given atomic meshes
TPZCompMesh * HybridizeHDiv::CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_Hybrid, double Lagrange_term_multiplier ){
    
}
    
    /// group and condense the elements
 void HybridizeHDiv::GroupElements(TPZCompMesh *cmesh_Hybrid){
    
}
    
    /// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes
void HybridizeHDiv::InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
}
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh
void HybridizeHDiv::InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier  ){
    
}
    
    /// clones the atomic meshes in meshvec_HDiv and creates a multiphysics hybrid mesh
std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > HybridizeHDiv::Hybridize(TPZCompMesh *cmesh_Multiphysics, TPZVec<TPZCompMesh *> &meshvec_HDiv, bool group_elements, double Lagrange_term_multiplier ){
    
    //copia a malha de pressao e fluxo
    TPZManVector<TPZCompMesh *> meshvec_Hybrid(meshvec_HDiv.size(),0);
    for (int i = 0; i<meshvec_HDiv.size(); i++) {
        meshvec_Hybrid[i] = meshvec_HDiv[i]->Clone();
    }
    
    //Nuevos Materiales
    fHDivWrapMatid = 100;
    fLagrangeInterface = 101;
    fInterfaceMatid = 102;
    
    //Insertando Nuevos Materiales
    TPZCompMesh *pressure_mesh = meshvec_Hybrid[1];
    TPZFNMatrix<1, STATE> xk(fNState, fNState, 0.), xb(fNState, fNState, 0.), xc(fNState, fNState, 0.), xf(fNState, 1, 0.);

    int dim_p= pressure_mesh->Dimension();
    if (dim_p ==2) {
        TPZMat1dLin *mat= new TPZMat1dLin(fLagrangeInterface);
        mat->SetMaterial(xk, xc, xb, xf);
        pressure_mesh->InsertMaterialObject(mat);
    }
    if (dim_p ==3) {
        TPZMat2dLin *mat= new TPZMat2dLin(fLagrangeInterface);
        mat->SetMaterial(xk, xc, xf);
        pressure_mesh->InsertMaterialObject(mat);
    }
    
    TPZCompMesh *flux_mesh = meshvec_Hybrid[0];
    
    int dim_f= flux_mesh->Dimension();
    if (dim_f ==2) {
        TPZMat1dLin *mat= new TPZMat1dLin(fHDivWrapMatid);
        mat->SetMaterial(xk, xc, xb, xf);
        flux_mesh->InsertMaterialObject(mat);
    }
    if (dim_f ==3) {
        TPZMat2dLin *mat= new TPZMat2dLin(fHDivWrapMatid);
        mat->SetMaterial(xk, xc, xf);
        flux_mesh->InsertMaterialObject(mat);
    }
    
    //HybridizeInternalSides
    TPZGeoMesh *gmesh = flux_mesh->Reference();
    gmesh->ResetReference();
    flux_mesh->LoadReferences();
    
    int64_t nel = flux_mesh->NElements();
    std::list<std::tuple<int64_t, int> > pressures;
    
    for (int iel=0; iel<nel; iel++) {
        TPZCompEl *  cel = flux_mesh->Element(iel);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel || intel->Reference()->Dimension() != dim_f) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side<gel->NSides(); side++) {
            if (gel->SideDimension(side) != dim_f - 1) {
                continue;
            }
           
            TPZCompElSide celside(intel,side);
            TPZCompElSide neighcomp = RightElement(intel, side);
            if (neighcomp) {
                //
            }
        }
    }
}
    
    /// verify the consistency of the solution of the flux mesh
void HybridizeHDiv::VerifySolutionConsistency(TPZCompMesh *fluxmesh, std::ostream &out){
    
}
    
std::tuple<int64_t,int> HybridizeHDiv::SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid){
    
        if (fHDivWrapMatid == 0 || fLagrangeInterface == 0) {
            std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
            DebugStop();
        }
        TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
        
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    
    TPZInterpolatedElement * intleft = dynamic_cast<TPZInterpolatedElement *>(left.Element());
    TPZInterpolatedElement * intright = dynamic_cast<TPZInterpolatedElement *>(right.Element());
    
    intleft->SetSideOrient(left.Side(), 1);
    intright->SetSideOrient(right.Side(), 1);
    
    
    
    
}
    

TPZCompElSide HybridizeHDiv::RightElement(TPZInterpolatedElement *intel, int side){
    
    bool isrestrained = false;
    TPZConnect &c = intel->SideConnect(0, side);
    if (c.HasDependency()) {
        isrestrained = true;
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel,side);
    
    if (isrestrained == true) {
         TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
        if (!celside) DebugStop();
        TPZGeoEl *neigh = celside.Element()->Reference();
        /// we assume that a larger element should not be a boundary element
        if (neigh->Dimension() != gel->Dimension()) {
            DebugStop();
        }
        return celside;
    }
    else{
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 1);
        if (celstack.size()==1) {
            TPZGeoEl *neigh = celstack[0].Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                return celstack[0];
            }
        }

    }
    return TPZCompElSide();
}
    
