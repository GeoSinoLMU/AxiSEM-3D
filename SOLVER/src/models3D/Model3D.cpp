//
//  Model3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of all 3D models

#include "Model3D.hpp"
#include "Quad.hpp"

// compute spz in reference
const eigen::DMatX3 &Model3D::computeRefSPZ(Quad &quad) {
    // use static to avoid repeated computation
    static int quadTag = -1;
    static eigen::DMatX3 spzRef;
    if (quadTag != quad.getGlobalTag()) {
        // quad nr
        const eigen::IRowN &pointNr = quad.getPointNr();
        // GLL coordinates
        const eigen::DMat2N &pointSZ = quad.getPointSZ();
        // allocate
        spzRef = eigen::DMatX3(pointNr.sum(), 3);
        // structured to flattened
        int row = 0;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int nr = pointNr(ipnt);
            // identical (s, z)
            spzRef.block(row, 0, nr, 1).fill(pointSZ(0, ipnt));
            spzRef.block(row, 2, nr, 1).fill(pointSZ(1, ipnt));
            // linearly varying phi
            spzRef.block(row, 1, nr, 1) =
            eigen::DColX::LinSpaced(nr, 0, 2. * numerical::dPi / nr * (nr - 1));
            row += nr;
        }
        quadTag = quad.getGlobalTag();
    }
    return spzRef;
}


#include "Volumetric3D.hpp"
#include "Geometric3D.hpp"
#include "OceanLoad3D.hpp"
#include "inparam.hpp"
#include "timer.hpp"

// build from inparam
void Model3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             std::vector<std::shared_ptr<const Model3D>> &models3D,
             bool rebuilding) {
    // check empty
    if (inparam::gInparamModel.
        get<std::string>("list_of_3D_models") == "NONE") {
        return;
    }
    
    // count
    int modelCount = inparam::gInparamModel.get<int>("list_of_3D_models:[?]");
    
    // loop over model list
    int mindexActive = -1;
    for (int mindex = 0; mindex < modelCount; mindex++) {
        // model name and key in inparam
        std::string keyInparam = "list_of_3D_models";
        const std::string &modelName = inparam::gInparamModel.
        get<std::string>(keyInparam + ":{" + bstring::toString(mindex) + "}");
        keyInparam += ":[" + bstring::toString(mindex) + "]:" + modelName;
        
        // start building
        timer::gPreloopTimer.begin("Building 3D model: " + modelName);
        
        // activated or not
        if (!inparam::gInparamModel.
            getWithDefault(keyInparam + ":activated", true)) {
            timer::gPreloopTimer.message("Model is deactivated.");
            timer::gPreloopTimer.ended("Building 3D model: " + modelName);
            continue;
        }
        mindexActive++;
        
        // do not rebuild if model is mpi-independent
        if (rebuilding) {
            if (!models3D[mindexActive]->isMPI_Dependent()) {
                timer::gPreloopTimer.message("Model is MPI-independent, "
                                             "no need to rebuild.");
                timer::gPreloopTimer.ended("Building 3D model: " + modelName);
                continue;
            }
        }
        
        // try volumetric first
        std::shared_ptr<const Model3D>
        model = Volumetric3D::buildInparam(exodusMesh, localMesh,
                                           modelName, keyInparam);
        // try geometric
        if (model == nullptr) {
            model = Geometric3D::buildInparam(exodusMesh, localMesh,
                                              modelName, keyInparam);
        }
        
        // try ocean load
        if (model == nullptr) {
            model = OceanLoad3D::buildInparam(exodusMesh, localMesh,
                                              modelName, keyInparam);
        }
        
        // unknown model class
        if (model == nullptr) {
            throw std::runtime_error("Model3D::buildInparam || "
                                     "Unknown 3D model: " + modelName);
        }
        
        // replace or push back
        if (rebuilding) {
            models3D[mindexActive].reset(model.get());
        } else {
            models3D.push_back(model);
        }
        timer::gPreloopTimer.ended("Building 3D model: " + modelName);
    }
}
