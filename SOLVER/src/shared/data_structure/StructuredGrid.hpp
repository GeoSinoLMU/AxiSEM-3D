//
//  StructuredGrid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  data on a structured grid
//  class parameters:
//  - RANK: number of dimensions
//  - DTYPE: datatype of data

#ifndef StructuredGrid_hpp
#define StructuredGrid_hpp

#include "NetCDF_Reader.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "eigen_tools.hpp"
#include "vector_tools.hpp"

template <int D, typename T>
class StructuredGrid {
public:
    // constructor
    StructuredGrid(const std::string &ncfile,
                   const std::array<std::string, D> &coordVarNames,
                   const std::vector<std::pair<std::string, double>> &varInfo,
                   const std::array<int, D> &shuffleData):
    mNV((int)varInfo.size()) {
        // verify shuffle data
        std::array<int, D> shuffleSorted = shuffleData;
        std::sort(shuffleSorted.begin(), shuffleSorted.end());
        for (int idim = 0; idim < D; idim++) {
            if (shuffleSorted[idim] != idim) {
                throw std::runtime_error("StructuredGrid::StructuredGrid || "
                                         "Invalid data ranks to shuffle.");
            }
        }
        
        // read data
        timer::gPreloopTimer.begin("Reading grid data");
        timer::gPreloopTimer.message("data file: " + io::popInputDir(ncfile));
        std::vector<std::vector<double>> gridCoords;
        if (mpi::root()) {
            // open
            NetCDF_Reader reader(io::popInputDir(ncfile));
            
            // read coords
            for (const std::string &cvName: coordVarNames) {
                std::vector<double> crd;
                reader.readVectorDouble(cvName, crd);
                gridCoords.push_back(crd);
            }
            
            // coord dimensions
            Eigen::array<Eigen::DenseIndex, D> dimsCrd;
            for (int idim = 0; idim < D; idim++) {
                dimsCrd[idim] = gridCoords[idim].size();
            }
            
            // allocate data
            Eigen::array<Eigen::DenseIndex, 1 + D> dimsData;
            dimsData[0] = mNV; // # variables comes first
            std::copy(dimsCrd.begin(), dimsCrd.end(), dimsData.begin() + 1);
            mGridData = Eigen::Tensor<T, 1 + D, Eigen::RowMajor>(dimsData);
            
            // shape
            Eigen::array<Eigen::DenseIndex, 1 + D> start, count;
            start.fill(0);
            count = mGridData.dimensions();
            count[0] = 1;
            
            // read data
            for (int ivar = 0; ivar < mNV; ivar++) {
                // read to double
                Eigen::Tensor<double, D, Eigen::RowMajor> dGridData;
                reader.readTensorDouble(varInfo[ivar].first, dGridData);
                dGridData = dGridData * varInfo[ivar].second;
                // round
                if (std::is_integral<T>::value) {
                    dGridData = dGridData.round();
                }
                // copy data
                start[0] = ivar;
                mGridData.slice(start, count).reshape(dimsCrd)
                = dGridData.shuffle(shuffleData).template cast<T>();
            }
            
            // check
            for (int idim = 0; idim < D; idim++) {
                // size
                if (gridCoords[idim].size() != mGridData.dimension(idim + 1)) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Incompatible grid sizes in coordinates and data. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncfile);
                }
                
                // size
                if (gridCoords[idim].size() < 2) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Too few grid points; at least two are required. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncfile);
                }
                
                // sorted
                if (!vector_tools::isSortedUnique(gridCoords[idim])) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Coordinates are not ascendingly sorted. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncfile);
                }
            }
        }
        timer::gPreloopTimer.ended("Reading grid data");
        
        // broadcast
        timer::gPreloopTimer.begin("Broadcasting grid data");
        mpi::bcast(gridCoords);
        // cast to array
        for (int idim = 0; idim < D; idim++) {
            mGridCoords[idim] = gridCoords[idim];
        }
        
        // below we broadcast the tensor manually
        // add mpi::bcastEigenTensor if needed in the future
        Eigen::array<Eigen::DenseIndex, 1 + D> dims;
        dims[0] = mNV;
        for (int idim = 0; idim < D; idim++) {
            dims[idim + 1] = mGridCoords[idim].size();
        }
        // allocate
        if (mpi::rank() != 0) {
            mGridData.resize(dims);
        }
        // broadcast tensor data
        mpi::bcast(mGridData.data(), (int)mGridData.size());
        
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfoTensor(mGridData, "structured_nu"));
        timer::gPreloopTimer.ended("Broadcasting grid data");
    }
    
    // in scope
    template <typename csArray>
    bool inScope(const csArray &targetCoords) const {
        // linear interp 0-1 points
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        for (int idim = 0; idim < D; idim++) {
            try {
                vector_tools::
                linearInterpSorted(mGridCoords[idim], targetCoords[idim],
                                   index0, index1, factor0, factor1);
            } catch (...) {
                // location out of range
                return false;
            }
        }
        return true;
    }
    
    // compute
    typedef Eigen::Matrix<T, 1, Eigen::Dynamic> TRowV;
    template <typename csArray>
    TRowV compute(const csArray &targetCoords,
                  const TRowV &valueForOutOfRange) const {
        // linear interp 0-1 points
        static std::array<std::array<int, D>, 2> index01;
        static std::array<std::array<double, D>, 2> factor01;
        for (int idim = 0; idim < D; idim++) {
            try {
                vector_tools::
                linearInterpSorted(mGridCoords[idim], targetCoords[idim],
                                   index01[0][idim], index01[1][idim],
                                   factor01[0][idim], factor01[1][idim]);
            } catch (...) {
                // location out of range
                return valueForOutOfRange;
            }
        }
        
        // shapes
        static Eigen::array<Eigen::DenseIndex, 1> shapeV = {mNV};
        static Eigen::array<Eigen::DenseIndex, 1 + D> start, count;
        start.fill(0);
        count.fill(1);
        count[0] = mNV;
        
        // loop over 2^D points
        typedef Eigen::Matrix<double, 1, Eigen::Dynamic> DRowV;
        DRowV result = DRowV::Zero(mNV);
        for (int deci = 0; deci < pow(2, D); deci++) {
            int left = deci;
            double factorCombine = 1.;
            for (int idim = 0; idim < D; idim++) {
                start[idim + 1] = index01[left % 2][idim];
                factorCombine *= factor01[left % 2][idim];
                left /= 2;
            }
            const Eigen::Tensor<double, 1, Eigen::RowMajor> &row =
            factorCombine * mGridData.slice(start, count).reshape(shapeV)
            .template cast<double>();
            result += Eigen::Map<const DRowV>(row.data(), mNV);
        }
        
        // round
        if (std::is_integral<T>::value) {
            result = result.array().round();
        }
        return result.template cast<T>();
    }
    
    // compute single
    template <typename csArray>
    T compute(const csArray &targetCoords, T valueForOutOfRange,
              int varIndex = 0) const {
        TRowV outVal = TRowV::Zero(mNV);
        outVal(varIndex) = valueForOutOfRange;
        return compute(targetCoords, outVal)(varIndex);
    }
    
    // get
    const std::array<std::vector<double>, D> &getGridCoords() const {
        return mGridCoords;
    }
    
    // get
    std::array<std::vector<double>, D> &getGridCoords() {
        return mGridCoords;
    }
    
    // get
    const Eigen::Tensor<T, 1 + D, Eigen::RowMajor> &getGridData() const {
        return mGridData;
    }
    
private:
    // grid coords
    std::array<std::vector<double>, D> mGridCoords;
    // grid data
    Eigen::Tensor<T, 1 + D, Eigen::RowMajor> mGridData;
    // number of variables
    const int mNV;
};

#endif /* StructuredGrid_hpp */
