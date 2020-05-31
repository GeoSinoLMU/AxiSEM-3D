//
//  StructuredGridO3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models based on structured grid

#include "StructuredGridO3D.hpp"

// constructor
StructuredGridO3D::
StructuredGridO3D(const std::string &modelName, const std::string &fname,
                  const std::array<std::string, 2> &crdVarNames,
                  const std::array<int, 2> &shuffleData,
                  bool sourceCentered, bool ellipticity,
                  double lengthUnit, double angleUnit,
                  const std::string &dataVarName, double factor):
OceanLoad3D(modelName), mFileName(fname), mCrdVarNames(crdVarNames),
mSourceCentered(sourceCentered), mEllipticity(ellipticity),
mDataVarName(dataVarName), mFactor(factor) {
    // init grids
    std::vector<std::pair<std::string, double>> varInfo;
    varInfo.push_back({mDataVarName, mFactor});
    mGrid = std::make_unique<StructuredGrid<2, double>>
    (mFileName, mCrdVarNames, varInfo, shuffleData);
    // coordinate units
    sg_tools::constructGridUnits(*mGrid, mSourceCentered, false,
                                 lengthUnit, angleUnit);
    // longitude range
    mLon360 = sg_tools::constructLon360(*mGrid, mSourceCentered, mModelName);
}

// get sum(rho * depth)
bool StructuredGridO3D::getSumRowDepth(const eigen::DMatX3 &spz,
                                       const eigen::DMat24 &nodalSZ,
                                       eigen::DColX &sumRowDepth) const {
    //////////////////////// coords ////////////////////////
    // check min/max
    const auto &gridCrds = mGrid->getGridCoords();
    if (!sg_tools::coordsScope1D
        (nodalSZ, false, false,
         mSourceCentered, gridCrds[0].front(), gridCrds[0].back(),
         false, 0., 0.)) {
        return false;
    }
    
    // compute grid coords
    const eigen::DMatX3 &crdGrid = sg_tools::
    coordsToGrid(spz, mSourceCentered, mEllipticity, mLon360,
                 false, false, mModelName);
    
    //////////////////////// values ////////////////////////
    // allocate and fill with zero
    int nCardinals = (int)spz.rows();
    sumRowDepth = eigen::DColX::Zero(nCardinals);
    
    // point loop
    double err = std::numeric_limits<double>::lowest();
    bool oneInScope = false;
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        const eigen::DRow2 &horizontal = crdGrid.block(ipnt, 0, 1, 2);
        double val = mGrid->compute(horizontal, err);
        // check scope
        if (val > err * .9) {
            sumRowDepth(ipnt) = val;
            oneInScope = true;
        }
    }
    return oneInScope;
}

// verbose
std::string StructuredGridO3D::verbose() const {
    std::stringstream ss;
    // head
    ss << sg_tools::verboseHead(mModelName, "StructuredGridO3D", mFileName);
    
    // coords
    const auto &gcrds = mGrid->getGridCoords();
    ss << sg_tools::verboseCoords(mSourceCentered, false, false,
                                  {mCrdVarNames[0], mCrdVarNames[1]},
                                  {gcrds[0].front(), gcrds[1].front()},
                                  {gcrds[0].back(), gcrds[1].back()});
    
    // options
    int width = 15;
    if (!mSourceCentered) {
        width = 22;
        ss << bstring::boxEquals(4, width, "ellipticity correction",
                                 mEllipticity);
    }
    
    // data
    ss << bstring::boxSubTitle(2, "Data for sum(rho * depth)");
    ss << bstring::boxEquals(4, width, "NetCDF variable", mDataVarName);
    typedef Eigen::Tensor<double, 0, Eigen::RowMajor> Tensor0;
    double min = ((Tensor0)mGrid->getGridData().minimum())(0);
    double max = ((Tensor0)mGrid->getGridData().maximum())(0);
    ss << bstring::boxEquals(4, width, "data range", bstring::range(min, max));
    return ss.str();
}
