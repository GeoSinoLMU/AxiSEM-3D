//
//  StructuredGridV3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D volumetric models based on structured grid

#include "StructuredGridV3D.hpp"

// constructor
StructuredGridV3D::
StructuredGridV3D(const std::string &modelName, const std::string &fname,
                  const std::array<std::string, 3> &crdVarNames,
                  const std::array<int, 3> &shuffleData,
                  bool sourceCentered, bool ellipticity,
                  bool useDepth, bool depthSolid, bool undulated,
                  double lengthUnit, double angleUnit, bool center,
                  const std::vector<std::tuple<std::string, std::string,
                  double, ReferenceKind>> &propertyInfo):
Volumetric3D(modelName), mFileName(fname), mCrdVarNames(crdVarNames),
mSourceCentered(sourceCentered), mEllipticity(ellipticity),
mUseDepth(useDepth), mDepthSolid(depthSolid), mUndulatedGeometry(undulated),
mElementCenter(center) {
    // read info
    std::transform(propertyInfo.begin(), propertyInfo.end(),
                   std::back_inserter(mPropertyKeys),
                   [](const auto &tuple) {return std::get<0>(tuple);});
    std::transform(propertyInfo.begin(), propertyInfo.end(),
                   std::back_inserter(mPropertyVarNames),
                   [](const auto &tuple) {return std::get<1>(tuple);});
    std::transform(propertyInfo.begin(), propertyInfo.end(),
                   std::back_inserter(mPropertyFactors),
                   [](const auto &tuple) {return std::get<2>(tuple);});
    std::transform(propertyInfo.begin(), propertyInfo.end(),
                   std::back_inserter(mPropertyReferenceKinds),
                   [](const auto &tuple) {return std::get<3>(tuple);});
    
    // verify depth
    if (mUndulatedGeometry && mUseDepth) {
        throw std::runtime_error
        ("StructuredGridV3D::StructuredGridV3D || "
         "Cannot use depth as vertical coordinate in undulated geometry. || "
         "Model name = " + mModelName);
    }
    
    // verify full anisotropy
    for (const std::string &key: mPropertyKeys) {
        if (key.front() == 'C') {
            // if any CIJ appreas, all CIJ must appear
            const std::vector<std::string> &allCIJ = {
                "C11", "C12", "C13", "C14", "C15", "C16",
                "C22", "C23", "C24", "C25", "C26",
                "C33", "C34", "C35", "C36",
                "C44", "C45", "C46",
                "C55", "C56",
                "C66"};
            mIndexCIJ = std::make_unique<eigen::IMat66>();
            for (const std::string &CIJ: allCIJ) {
                auto it = std::find(mPropertyKeys.begin(),
                                    mPropertyKeys.end(), CIJ);
                if (it == mPropertyKeys.end()) {
                    throw std::runtime_error
                    ("StructuredGridV3D::StructuredGridV3D || "
                     "All 21 elements in Cijkl must be presented. || "
                     "Missing " + CIJ + ". || "
                     "Model name = " + mModelName);
                }
                int i = bstring::cast<int>
                (CIJ.substr(1, 1), "StructuredGridV3D::StructuredGridV3D") - 1;
                int j = bstring::cast<int>
                (CIJ.substr(2, 1), "StructuredGridV3D::StructuredGridV3D") - 1;
                (*mIndexCIJ)(i, j) = (int)(it - mPropertyKeys.begin());
                (*mIndexCIJ)(j, i) = (int)(it - mPropertyKeys.begin());
            }
            break;
        }
    }
    
    // init grid
    std::vector<std::pair<std::string, double>> data;
    for (int iprop = 0; iprop < propertyInfo.size(); iprop++) {
        data.push_back({mPropertyVarNames[iprop], mPropertyFactors[iprop]});
    }
    mGrid = std::make_unique<StructuredGrid<3, double>>
    (mFileName, mCrdVarNames, data, shuffleData);
    // coordinate units
    sg_tools::constructUnits(*mGrid, mSourceCentered, true,
                             lengthUnit, angleUnit);
    // longitude range
    mLon360 = sg_tools::constructLon360(*mGrid, mSourceCentered, mModelName);
}

// get properties
bool StructuredGridV3D::getProperties(const eigen::DMatX3 &spz,
                                      const eigen::DMat24 &nodalSZ,
                                      std::vector<std::string> &propKeys,
                                      std::vector<ReferenceKind> &refKinds,
                                      eigen::IMatXX &inScopes,
                                      eigen::DMatXX &propValues) const {
    //////////////////////// coords ////////////////////////
    // check center
    const auto &gridCrds = mGrid->getGridCoords();
    if (mElementCenter) {
        if (!inplaneScope<eigen::DCol2>
            (nodalSZ.rowwise().mean(),
             mSourceCentered, gridCrds[0].front(), gridCrds[0].back(),
             true, gridCrds[2].front(), gridCrds[2].back(),
             mUseDepth, mDepthSolid)) {
            return false;
        }
    }
    
    // check min/max
    if (!inplaneScope(nodalSZ,
                      mSourceCentered, gridCrds[0].front(), gridCrds[0].back(),
                      true, gridCrds[2].front(), gridCrds[2].back(),
                      mUseDepth, mDepthSolid)) {
        return false;
    }
    
    // compute grid coords
    eigen::DMatX3 crdGrid =
    coordsFromMeshToModel(spz, mSourceCentered, mEllipticity, mLon360,
                          mUseDepth, mDepthSolid, mModelName);
    
    //////////////////////// values ////////////////////////
    // header
    propKeys = mPropertyKeys;
    refKinds = mPropertyReferenceKinds;
    
    // allocate
    int nProperties = (int)mPropertyKeys.size();
    int nCardinals = (int)spz.rows();
    inScopes = eigen::IMatXX::Zero(nCardinals, nProperties);
    propValues = eigen::DMatXX::Zero(nCardinals, nProperties);
    
    // point loop
    const static double err = std::numeric_limits<double>::lowest();
    const eigen::DRowX &valOut = eigen::DRowX::Constant(nProperties, err);
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        const eigen::DRowX &val = mGrid->compute(crdGrid.row(ipnt), valOut);
        // check scope
        if (val(0) > err * .9) {
            inScopes.row(ipnt).fill(1);
            propValues.row(ipnt) = val;
        }
    }
    
    // rotate full anisotropy
    if (!mSourceCentered && mIndexCIJ) {
        // get depth back to radius
        if (mUseDepth) {
            double R = (mDepthSolid ? geodesy::getOuterSolidRadius() :
                        geodesy::getOuterRadius());
            crdGrid.col(2).array() = R - crdGrid.col(2).array();
        }
        // compute back azimuth
        const eigen::DColX &baz = geodesy::backAzimuth(crdGrid, mEllipticity);
        
        // point loop
        for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
            // copy input
            eigen::DMat66 inCijkl = eigen::DMat66::Zero();
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    inCijkl(i, j) = propValues(ipnt, (*mIndexCIJ)(i, j));
                }
            }
            // transform
            const eigen::DMat66 &outCijkl =
            bondTransformation(inCijkl, 0., 0., -baz(ipnt));
            // copy output
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    propValues(ipnt, (*mIndexCIJ)(i, j)) = outCijkl(i, j);
                }
            }
        }
    }
    
    // return true if any point/property is in scope
    return (bool)inScopes.maxCoeff();
}

// verbose
std::string StructuredGridV3D::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    // head
    ss << sg_tools::verboseHead(mModelName, "StructuredGridV3D", mFileName);
    
    // coords
    const auto &gcrds = mGrid->getGridCoords();
    ss << sg_tools::
    verboseCoords(mSourceCentered, true, mUseDepth,
                  {mCrdVarNames[0], mCrdVarNames[1], mCrdVarNames[2]},
                  {gcrds[0].front(), gcrds[1].front(), gcrds[2].front()},
                  {gcrds[0].back(), gcrds[1].back(), gcrds[2].back()});
    
    // options
    if (!mSourceCentered) {
        ss << boxEquals(4, 29, "ellipticity correction", mEllipticity);
    }
    if (mUseDepth) {
        ss << boxEquals(4, 29, "depth below solid surface", mDepthSolid);
    } else {
        ss << boxEquals(4, 29, "radial geometry", mUndulatedGeometry ?
                        "undulated" : "reference");
    }
    ss << boxEquals(4, 29, "force element center in scope", mElementCenter);
    
    // properties
    ss << boxSubTitle(2, "Properties");
    // width
    int widthK = std::max(vector_tools::maxLength(mPropertyKeys), 3);
    int widthV = std::max(vector_tools::maxLength(mPropertyVarNames), 6);
    // table
    ss << "      " << std::left;
    ss << std::setw(widthK) << "KEY" << " | ";
    ss << std::setw(widthV) << "NC-VAR" << " | ";
    ss << "REF KIND" << " | " << "RANGE" << "\n";
    const auto &minMax = mGrid->getDataRange();
    for (int iprop = 0; iprop < mPropertyKeys.size(); iprop++) {
        ss << "    * " << std::left;
        ss << std::setw(widthK) << mPropertyKeys[iprop] << " | ";
        ss << std::setw(widthV) << mPropertyVarNames[iprop] << " | ";
        ss << std::setw(8) <<
        sReferenceKindStr.at(mPropertyReferenceKinds[iprop]) << " | ";
        ss << range(minMax(iprop, 0), minMax(iprop, 1)) << "\n";
    }
    return ss.str();
}
