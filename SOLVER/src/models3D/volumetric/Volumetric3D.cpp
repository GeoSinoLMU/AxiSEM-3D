//
//  Volumetric3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D volumetric models
//  density, velocity, elasticity, attenuation

#include "Volumetric3D.hpp"
#include "Quad.hpp"

// apply to Quad
void Volumetric3D::applyTo(Quad &quad) const {
    // cardinal coordinates in reference geometry
    eigen::DMatX3 spzUse = computeRefSPZ(quad);
    if (usingUndulatedGeometry()) {
        // get undulation from quad
        eigen::arN_DColX und = quad.getUndulation();
        // structured to flattened
        int row = 0;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int nr =  quad.getPointNr()(ipnt);
            // extend to 3D
            if (und[ipnt].rows() == 1) {
                und[ipnt] = und[ipnt].replicate(nr, 1);
            }
            // change coords based on CS type
            if (geodesy::isCartesian()) {
                // add undulation to z in Cartesian
                spzUse.block(row, 2, nr, 1) += und[ipnt];
            } else {
                // add undulation to r in spherical
                double s = spzUse(row, 0);
                double z = spzUse(row, 2);
                double theta = acos(z / std::max(sqrt(s * s + z * z),
                                                 numerical::dEpsilon));
                spzUse.block(row, 0, nr, 1) += und[ipnt] * sin(theta);
                spzUse.block(row, 2, nr, 1) += und[ipnt] * cos(theta);
            }
            row += nr;
        }
    }
    
    // get properties
    std::vector<std::string> propKeys;
    std::vector<ReferenceKind> refKinds;
    eigen::IMatXX inScopes;
    eigen::DMatXX propValues;
    if (getProperties(spzUse, quad.getNodalSZ(), propKeys,
                      refKinds, inScopes, propValues)) {
        // property loop
        int nprop = (int)propKeys.size();
        for (int iprop = 0; iprop < nprop; iprop++) {
            // flattened to structured
            eigen::arN_IColX inScope;
            eigen::arN_DColX propValue;
            int row = 0;
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                int nr =  quad.getPointNr()(ipnt);
                inScope[ipnt] = inScopes.block(row, iprop, nr, 1);
                propValue[ipnt] = propValues.block(row, iprop, nr, 1);
                row += nr;
            }
            // set to material
            quad.getMaterialPtr()->
            addProperty3D(propKeys[iprop], refKinds[iprop], inScope, propValue);
        }
    }
}

// Bond transformation for rotating Cijkl
eigen::DMat66 Volumetric3D::bondTransformation(const eigen::DMat66 &inCijkl,
                                               double alpha, double beta,
                                               double gamma) {
    // R
    eigen::DMat33 R1, R2, R3, R;
    R1 << 1., 0., 0., 0., cos(alpha), sin(alpha), 0., -sin(alpha), cos(alpha);
    R2 << cos(beta), 0., sin(beta), 0., 1., 0., -sin(beta), 0, cos(beta);
    R3 << cos(gamma), sin(gamma), 0., -sin(gamma), cos(gamma), 0., 0., 0., 1.;
    R = R1 * R2 * R3;
    
    // K
    eigen::DMat33 K1, K2, K3, K4;
    K1.array() = R.array().pow(2.);
    K2 <<
    R(0, 1) * R(0, 2), R(0, 2) * R(0, 0), R(0, 0) * R(0, 1),
    R(1, 1) * R(1, 2), R(1, 2) * R(1, 0), R(1, 0) * R(1, 1),
    R(2, 1) * R(2, 2), R(2, 2) * R(2, 0), R(2, 0) * R(2, 1);
    K3 <<
    R(1, 0) * R(2, 0), R(1, 1) * R(2, 1), R(1, 2) * R(2, 2),
    R(2, 0) * R(0, 0), R(2, 1) * R(0, 1), R(2, 2) * R(0, 2),
    R(0, 0) * R(1, 0), R(0, 1) * R(1, 1), R(0, 2) * R(1, 2);
    K4 <<
    R(1, 1) * R(2, 2) + R(1, 2) * R(2, 1),
    R(1, 2) * R(2, 0) + R(1, 0) * R(2, 2),
    R(1, 0) * R(2, 1) + R(1, 1) * R(2, 0),
    R(2, 1) * R(0, 2) + R(2, 2) * R(0, 1),
    R(2, 2) * R(0, 0) + R(2, 0) * R(0, 2),
    R(2, 0) * R(0, 1) + R(2, 1) * R(0, 0),
    R(0, 1) * R(1, 2) + R(0, 2) * R(1, 1),
    R(0, 2) * R(1, 0) + R(0, 0) * R(1, 2),
    R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0);
    eigen::DMat66 K;
    K.block(0, 0, 3, 3) = K1;
    K.block(0, 3, 3, 3) = 2. * K2;
    K.block(3, 0, 3, 3) = K3;
    K.block(3, 3, 3, 3) = K4;
    
    // rotation
    return K * inCijkl * K.transpose();
}


#include "StructuredGridV3D.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const Volumetric3D> Volumetric3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             const std::string &modelName, const std::string &keyInparam) {
    // short alias
    const InparamYAML &gm = inparam::gInparamModel;
    const std::string &root = keyInparam;
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "StructuredGridV3D") {
        // file name
        const std::string &fname = gm.get<std::string>(root + ":nc_data_file");
        
        ////////////// coords //////////////
        const std::string &rootc = root + ":coordinates";
        // horizontal
        bool sourceCentered = false, ellipticity = false;
        std::array<std::string, 3> crdVarNames;
        std::array<int, 3> shuffleData;
        sg_tools::inparamHorizontal<3>(gm, rootc, modelName, className,
                                       sourceCentered, ellipticity,
                                       crdVarNames, shuffleData);
        // vertical
        bool useDepth = false, depthSolid = false;
        sg_tools::inparamVertical(gm, rootc, modelName, className,
                                  useDepth, depthSolid);
        const std::string &rd = useDepth ? "depth" : "radius";
        crdVarNames[2] = gm.get<std::string>(rootc + ":" + rd + ":nc_var");
        shuffleData[2] = gm.get<int>(rootc + ":" + rd + ":data_rank");
        bool undulated = false;
        if (!useDepth) {
            undulated =
            gm.getWithDefault(rootc + ":radius:undulated_geometry", false);
        }
        // units
        double lengthUnit = 1., angleUnit = 1.;
        sg_tools::inparamUnits(gm, rootc, lengthUnit, angleUnit);
        // center
        bool center =
        gm.getWithDefault(rootc + ":element_center_in_scope", false);
        
        ////////////// properties //////////////
        // size
        const std::string &rootp = root + ":properties";
        int nprop = gm.get<int>(rootp + ":[?]");
        // info
        std::vector<std::tuple<
        std::string, std::string, double, ReferenceKind>> propertyInfo;
        for (int iprop = 0; iprop < nprop; iprop++) {
            // get key
            const std::string &key = gm.get<std::string>
            (rootp + ":{" + bstring::toString(iprop) + "}");
            const std::string &rootpi =
            rootp + ":[" + bstring::toString(iprop) + "]:" + key;
            // var, factor, ref
            const std::string &vname = gm.get<std::string>(rootpi + ":nc_var");
            double factor = gm.getWithDefault(rootpi + ":factor", 1.);
            ReferenceKind ref = gm.getWithLimits<ReferenceKind>
            (rootpi + ":reference_kind", {
                {"ABS", ReferenceKind::ABS},
                {"REF1D", ReferenceKind::REF1D},
                {"REF3D", ReferenceKind::REF3D},
                {"REF_PERTURB", ReferenceKind::REF_PERTURB}});
            propertyInfo.push_back({key, vname, factor, ref});
        }
        
        // construct
        return std::make_shared
        <const StructuredGridV3D>(modelName, fname, crdVarNames, shuffleData,
                                  sourceCentered, ellipticity,
                                  useDepth, depthSolid, undulated,
                                  lengthUnit, angleUnit, center, propertyInfo);
    } else {
        // other models
    }
    
    // unknown class
    return nullptr;
}
