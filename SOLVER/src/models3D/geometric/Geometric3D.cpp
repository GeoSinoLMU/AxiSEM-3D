//
//  Geometric3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D geometric models

#include "Geometric3D.hpp"
#include "Quad.hpp"

// apply to Quad
void Geometric3D::applyTo(Quad &quad) const {
    // cardinal coordinates in reference geometry
    const eigen::DMatX3 &spzRef = computeRefSPZ(quad);
    
    // get undulation
    eigen::DColX unds;
    if (getUndulation(spzRef, quad.getNodalSZ(), unds)) {
        // flattened to structured
        eigen::arN_DColX und;
        int row = 0;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int nr = quad.getPointNr()(ipnt);
            und[ipnt] = unds.block(row, 0, nr, 1);
            row += nr;
        }
        // set to Quad
        quad.getUndulationPtr()->addUndulation(und);
    }
}


#include "StructuredGridG3D.hpp"
#include "Ellipticity.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const Geometric3D> Geometric3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             const std::string &modelName, const std::string &keyInparam) {
    // short alias
    const InparamYAML &gm = inparam::gInparamModel;
    const std::string &root = keyInparam;
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "StructuredGridG3D") {
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
        // vertical
        bool useDepth = false, depthSolid = false;
        sg_tools::inparamVertical(gm, rootc, modelName, className,
                                  useDepth, depthSolid);
        const std::string &rd = useDepth ? "depth" : "radius";
        double interface = gm.get<double>(rootc + ":" + rd + ":interface");
        double min = gm.get<double>(rootc + ":" + rd + ":min_max:[0]");
        double max = gm.get<double>(rootc + ":" + rd + ":min_max:[1]");
        if (interface >= max || interface <= min) {
            throw std::runtime_error
            ("Geometric3D::buildInparam || "
             + rd + ":interface must lie within " + rd + ":min_max."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
        // units
        double lengthUnit = 1., angleUnit = 1.;
        sg_tools::inparamUnits(gm, rootc, lengthUnit, angleUnit);
        
        ////////////// data //////////////
        const std::string &dataVarName =
        gm.get<std::string>(root + ":undulation_data:nc_var");
        double factor =
        gm.getWithDefault(root + ":undulation_data:factor", 1.);
        
        // construct
        return std::make_shared
        <const StructuredGridG3D>(modelName, fname, crdVarNames, shuffleData,
                                  sourceCentered, ellipticity,
                                  useDepth, depthSolid,
                                  interface, min, max, lengthUnit, angleUnit,
                                  dataVarName, factor);
    } else if (className == "Ellipticity") {
        // ellipticity can be added only once
        static bool ellipticityAdded = false;
        if (ellipticityAdded) {
            throw
            std::runtime_error("Geometric3D::buildInparam || Ellipticity "
                               "model cannot be added more than once.");
        }
        if (io::gVerboseWarnings) {
            if (geodesy::isCartesian()) {
                io::cout <<
                bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                 "model will be ignored for a Cartesian mesh.");
                
            }
            if (geodesy::getOuterFlattening() < numerical::dEpsilon) {
                io::cout <<
                bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                 "model will be ignored with zero flattening.");
                
            }
        }
        ellipticityAdded = true;
        return std::make_shared<const Ellipticity>(modelName);
    }
    
    // unknown class
    return nullptr;
}
