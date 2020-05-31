//
//  sg_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/6/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  tools for structured-grid based model

#ifndef sg_tools_hpp
#define sg_tools_hpp

#include "inparam.hpp"
#include "geodesy.hpp"
#include "StructuredGrid.hpp"
#include <iomanip> // required by gcc

namespace sg_tools {
    ////////////////////////// inparam //////////////////////////
    // horizontal coords
    template <int NDIM>
    void inparamHorizontal(const InparamYAML &gm, const std::string &rootc,
                           const std::string &modelName,
                           const std::string &className,
                           bool &sourceCentered, bool &ellipticity,
                           std::array<std::string, NDIM> &crdVarNames,
                           std::array<int, NDIM> &shuffleData) {
        if (gm.contains(rootc + ":distance")) {
            sourceCentered = true;
            crdVarNames[0] = gm.get<std::string>(rootc + ":distance:nc_var");
            crdVarNames[1] = gm.get<std::string>(rootc + ":azimuth:nc_var");
            shuffleData[0] = gm.get<int>(rootc + ":distance:data_rank");
            shuffleData[1] = gm.get<int>(rootc + ":azimuth:data_rank");
        } else if (gm.contains(rootc + ":latitude")) {
            sourceCentered = false;
            crdVarNames[0] = gm.get<std::string>(rootc + ":latitude:nc_var");
            crdVarNames[1] = gm.get<std::string>(rootc + ":longitude:nc_var");
            shuffleData[0] = gm.get<int>(rootc + ":latitude:data_rank");
            shuffleData[1] = gm.get<int>(rootc + ":longitude:data_rank");
            ellipticity =
            gm.getWithDefault(rootc + ":latitude:ellipticityticity", false);
        } else {
            throw std::runtime_error
            ("sg_tools::readHorizontal || "
             "Either coordinates:distance or coordinates:latitude must exist."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
    }
    
    // vertical coords
    inline void inparamVertical(const InparamYAML &gm, const std::string &rootc,
                                const std::string &modelName,
                                const std::string &className,
                                bool &useDepth, bool &depthSolid) {
        if (gm.contains(rootc + ":depth")) {
            useDepth = true;
            depthSolid =
            gm.getWithDefault(rootc + ":depth:below_solid_surface", true);
        } else if (gm.contains(rootc + ":radius")) {
            useDepth = false;
        } else {
            throw std::runtime_error
            ("sg_tools::readVertical || "
             "Either coordinates:radius or coordinates:depth must exist."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
    }
    
    // coord units
    inline void inparamUnits(const InparamYAML &gm, const std::string &rootc,
                             double &lengthUnit, double &angleUnit) {
        lengthUnit = gm.getWithOptions<double>(rootc + ":length_unit", {
            {"m", 1.}, {"km", 1e3}});
        angleUnit = gm.getWithLimits<double>(rootc + ":angle_unit", {
            {"degree", numerical::dDegree}, {"radian", 1.}});
    }
    
    
    ////////////////////////// construct //////////////////////////
    template <class Grid>
    void constructGridUnits(Grid &grid, bool sourceCentered, bool vertical,
                            double lengthUnit, double angleUnit) {
        // transform lambda expr
        auto vectorTimesFactor = [](std::vector<double> &vec, double factor) {
            std::transform(vec.begin(), vec.end(), vec.begin(),
                           [factor](double x) -> double {return factor * x;});};
        
        // horizontal
        if (sourceCentered) {
            // R
            if (geodesy::isCartesian()) {
                vectorTimesFactor(grid.getGridCoords()[0], lengthUnit);
            } else {
                vectorTimesFactor(grid.getGridCoords()[0], angleUnit);
            }
            // T
            vectorTimesFactor(grid.getGridCoords()[1], angleUnit);
        } else {
            double angle = angleUnit / numerical::dDegree;
            // lat
            vectorTimesFactor(grid.getGridCoords()[0], angle);
            // lon
            vectorTimesFactor(grid.getGridCoords()[1], angle);
        }
        
        // vertical
        if (vertical) {
            vectorTimesFactor(grid.getGridCoords()[2], lengthUnit);
        }
    }
    
    template <class Grid>
    bool constructLon360(const Grid &grid, bool sourceCentered,
                         const std::string &modelName) {
        if (sourceCentered) {
            // not used in source-centered
            return false;
        }
        const std::vector<double> &lon = grid.getGridCoords()[1];
        if (lon.front() < 0. - numerical::dEpsilon &&
            lon.back() > 180. + numerical::dEpsilon) {
            throw std::runtime_error
            ("sg_tools::checkLon360 || "
             "Longitude range must be either [-180, 180] or [0, 360]. || "
             "Model name = " + modelName);
        }
        return (lon.back() > 180. + numerical::dEpsilon);
    }
    
    
    ////////////////////////// coords //////////////////////////
    // check scope 1D
    template <class Mat2X>
    bool coordsScope1D(const Mat2X &sz, bool useDepth, bool depthSolid,
                       bool checkR, double minR, double maxR,
                       bool checkZ, double minZ, double maxZ) {
        // compute grid CS
        Mat2X gcrd;
        if (geodesy::isCartesian()) {
            gcrd = sz;
        } else {
            // sz -> tr
            gcrd = geodesy::sz2rtheta(sz, false, 0, 1, 1, 0);
        }
        // radius -> depth
        if (useDepth) {
            double router = (depthSolid ? geodesy::getOuterSolidRadius() :
                             geodesy::getOuterRadius());
            gcrd.row(1) = router - gcrd.row(1).array();
        }
        // check R
        if (checkR) {
            if (gcrd.row(0).maxCoeff() < minR ||
                gcrd.row(0).minCoeff() > maxR) {
                return false;
            }
        }
        // check Z
        if (checkZ) {
            if (gcrd.row(1).maxCoeff() < minZ ||
                gcrd.row(1).minCoeff() > maxZ) {
                return false;
            }
        }
        return true;
    }
    
    // compute grid coords
    inline eigen::DMatX3
    coordsToGrid(const eigen::DMatX3 &spz,
                 bool sourceCentered, bool ellipticity, bool lon360,
                 bool useDepth, bool depthSolid,
                 const std::string &modelName) {
        eigen::DMatX3 crdGrid;
        if (sourceCentered) {
            if (geodesy::isCartesian()) {
                crdGrid = spz;
            } else {
                // spz -> RTZ
                crdGrid = geodesy::sz2rtheta(spz, true, 0, 2, 2, 0);
                crdGrid.col(1) = spz.col(1);
            }
        } else {
            // correct spz for Cartesian
            // NOTE: Geographic and Cartesian contradict each other.
            //       We correct spz using s as arc-length and z as radius.
            //       Without such "bending", sqrt(s*s+z*z) on the surface
            //       will exceed the outer radius.
            if (geodesy::isCartesian()) {
                // z must be significantly larger than s
                if ((spz.array().col(2) < spz.array().col(0) * 10.).any()) {
                    throw std::runtime_error
                    ("sg_tools::computeGridCoords || "
                     "Invalid geographic location in Cartesian mesh. || "
                     "Set PLANET_RADIUS in Salvus mesher data file (*.bm). || "
                     "Model name = " + modelName);
                }
                eigen::DMatX3 spzBend(spz.rows(), spz.cols());
                const auto &t = spz.col(0).cwiseQuotient(spz.col(2));
                spzBend.col(0) = spz.col(2).array() * t.array().sin();
                spzBend.col(2) = spz.col(2).array() * t.array().cos();
                spzBend.col(1) = spz.col(1);
                // spz to lat, lon, r
                crdGrid = geodesy::spz2llr(spzBend, ellipticity, lon360);
            } else {
                // spz to lat, lon, r
                crdGrid = geodesy::spz2llr(spz, ellipticity, lon360);
            }
        }
        
        // depth (must be in reference geometry)
        if (useDepth) {
            double R = (depthSolid ? geodesy::getOuterSolidRadius() :
                        geodesy::getOuterRadius());
            crdGrid.col(2).array() = R - crdGrid.col(2).array();
        }
        return crdGrid;
    }
    
    
    ////////////////////////// verbose //////////////////////////
    // head
    inline std::string verboseHead(const std::string &modelName,
                                   const std::string &className,
                                   const std::string &fileName) {
        std::stringstream ss;
        ss << bstring::boxSubTitle(0, modelName + " ", '~');
        ss << bstring::boxEquals(2, 11, "class name", className);
        ss << bstring::boxEquals(2, 11, "NetCDF file", fileName);
        return ss.str();
    }
    
    // coords
    inline std::string
    verboseCoords(bool sourceCentered, bool vertical, bool useDepth,
                  const std::vector<std::string> &crdVarNames,
                  const std::vector<double> &crdMin,
                  const std::vector<double> &crdMax) {
        std::stringstream ss;
        ss << bstring::boxSubTitle(2, "Coordinates");
        std::vector<std::string> crdNames;
        if (sourceCentered) {
            crdNames.push_back("distance");
            crdNames.push_back("azimuth");
        } else {
            crdNames.push_back("latitude");
            crdNames.push_back("longitude");
        }
        if (vertical) {
            crdNames.push_back(useDepth ? "depth" : "radius");
        }
        // width
        int widthCrd = std::max(vector_tools::maxLength(crdNames), 5);
        int widthVar = std::max(vector_tools::maxLength(crdVarNames), 6);
        // table
        ss << "      " << std::left;
        ss << std::setw(widthCrd) << "COORD" << " | ";
        ss << std::setw(widthVar) << "NC-VAR" << " | SCOPE" << "\n";
        for (int idim = 0; idim < crdNames.size(); idim++) {
            ss << "    * " << std::left;
            ss << std::setw(widthCrd) << crdNames[idim] << " | ";
            ss << std::setw(widthVar) << crdVarNames[idim] << " | ";
            ss << bstring::range(crdMin[idim], crdMax[idim]) << "\n";
        }
        return ss.str();
    }
}

#endif /* sg_tools_hpp */
