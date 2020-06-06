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
    // coord units
    template <class Grid>
    void constructUnits(Grid &grid, bool sourceCentered, bool vertical,
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
    
    // longitude range
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
            ("sg_tools::constructLon360 || "
             "Longitude range must be either [-180, 180] or [0, 360]. || "
             "Model name = " + modelName);
        }
        return (lon.back() > 180. + numerical::dEpsilon);
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
