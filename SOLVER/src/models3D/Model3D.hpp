//
//  Model3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/23/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of all 3D models

#ifndef Model3D_hpp
#define Model3D_hpp

#include <vector>
#include "eigen_sem.hpp"

class ExodusMesh;
class LocalMesh;
class Quad;

class Model3D {
public:
    // constructor
    Model3D(const std::string &modelName): mModelName(modelName) {
        // nothing
    }
    
    // destructor
    virtual ~Model3D() = default;
    
    // apply to Quad
    virtual void applyTo(Quad &quad) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // get model name
    const std::string &getModelName() const {
        return mModelName;
    }
    
protected:
    // mpi-dependent: rebuild after domain re-partitioning
    virtual bool isMPI_Dependent() const {
        return false;
    }
    
    // model name
    const std::string mModelName;
    
    
    ////////////////////////////// static //////////////////////////////
    // compute spz in reference
    static const eigen::DMatX3 &computeRefSPZ(Quad &quad);
    
public:
    // build from inparam
    static void
    buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                 std::vector<std::shared_ptr<const Model3D>> &models3D,
                 bool rebuilding);
};

#endif /* Model3D_hpp */
