//
//  OceanLoad3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models

#include "OceanLoad3D.hpp"
#include "Quad.hpp"
#include "vicinity.hpp"

// apply to Quad
void OceanLoad3D::applyTo(Quad &quad) const {
    // only surface
    int surfEdge = quad.getSurfaceEdge();
    if (surfEdge == -1) {
        return;
    }
    
    // nr and edge point
    const eigen::IRowN &pointNr = quad.getPointNr();
    const std::vector<int> &ipnts = vicinity::constants::gEdgeIPnt[surfEdge];
    
    // cardinal coordinates in reference geometry
    // use static to avoid repeated computation
    static int quadTag = -1;
    static eigen::DMatX3 spzRef;
    if (quadTag != quad.getGlobalTag()) {
        // GLL coordinates
        const eigen::DMat2N &pointSZ = quad.getPointSZ();
        // allocate
        spzRef = eigen::DMatX3(pointNr(ipnts).sum(), 3);
        // structured to flattened
        int row = 0;
        for (int ipnt: ipnts) {
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
    
    // get data
    eigen::DColX sumRhoDepth;
    if (getSumRowDepth(spzRef, quad.getNodalSZ(), sumRhoDepth)) {
        // flattened to structured
        eigen::arP_DColX sumRD;
        int row = 0;
        for (int ip = 0; ip < spectral::nPED; ip++) {
            int nr = pointNr(ipnts[ip]);
            sumRD[ip] = sumRhoDepth.block(row, 0, nr, 1);
            row += nr;
        }
        // set to Quad
        quad.getOceanLoadPtr()->addSumRhoDepth(sumRD);
    }
}


#include "StructuredGridO3D.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const OceanLoad3D> OceanLoad3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             const std::string &modelName, const std::string &keyInparam) {
    // short alias
    const InparamYAML &gm = inparam::gInparamModel;
    const std::string &root = keyInparam;
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "StructuredGridO3D") {
        // file name
        const std::string &fname = gm.get<std::string>(root + ":nc_data_file");
        
        ////////////// coords //////////////
        const std::string &rootc = root + ":coordinates";
        // horizontal
        bool sourceCentered = false, ellipticity = false;
        std::array<std::string, 2> crdVarNames;
        std::array<int, 2> shuffleData;
        sg_tools::inparamHorizontal<2>(gm, rootc, modelName, className,
                                       sourceCentered, ellipticity,
                                       crdVarNames, shuffleData);
        // units
        double lengthUnit = 1., angleUnit = 1.;
        sg_tools::inparamUnits(gm, rootc, lengthUnit, angleUnit);
        
        ////////////// data //////////////
        const std::string &rootd = root + ":data_sum_rho_depth";
        const std::string &dataVarName = gm.get<std::string>(rootd + ":nc_var");
        double factor = gm.getWithDefault(rootd + ":factor", 1.);
        
        // construct
        return std::make_shared
        <const StructuredGridO3D>(modelName, fname, crdVarNames, shuffleData,
                                  sourceCentered, ellipticity,
                                  lengthUnit, angleUnit, dataVarName, factor);
    } else {
        // other models
    }
    
    // unknown class
    return nullptr;
}
