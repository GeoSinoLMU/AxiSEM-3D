//
//  GLLPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  GLL point for preloop processing
//  generator of Point and boundary conditions in core

#include "GLLPoint.hpp"
#include "Domain.hpp"

// point
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "Mass1D.hpp"
#include "Mass3D.hpp"
#include "MassOceanLoad1D.hpp"
#include "MassOceanLoad3D.hpp"
#include "geodesy.hpp"

// boundary
#include "SolidFluidBoundary.hpp"
#include "SolidFluidCoupling1D.hpp"
#include "SolidFluidCoupling3D.hpp"
#include "ABC.hpp"
#include "AbsorbingBoundary.hpp"
#include "ClaytonSolid1D.hpp"
#include "ClaytonSolid3D.hpp"
#include "ClaytonFluid1D.hpp"
#include "ClaytonFluid3D.hpp"
#include "AxialBoundary.hpp"
#include "FluidSurfaceBoundary.hpp"

// release to domain
void GLLPoint::
release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain) {
    //////////////////////////// reduce ////////////////////////////
    op1D_3D::tryReduceTo1D(mMassFluid);
    op1D_3D::tryReduceTo1D(mMassSolid);
    op1D_3D::tryReduceTo1D(mNormalSFU);
    op1D_3D::tryReduceTo1D(mNormalSFA);
    for (auto itm = mClaytonABC.begin(); itm != mClaytonABC.end(); itm++) {
        for (auto itv = itm->second.begin(); itv != itm->second.end(); itv++) {
            op1D_3D::tryReduceTo1D(std::get<1>(*itv));
            op1D_3D::tryReduceTo1D(std::get<2>(*itv));
            op1D_3D::tryReduceTo1D(std::get<3>(*itv));
        }
    }
    op1D_3D::tryReduceTo1D(mNormalTop);
    op1D_3D::tryReduceTo1D(mSumRhoDepth);
    
    
    //////////////////////////// point ////////////////////////////
    // solid point
    if (mMassSolid.rows() > 0) {
        // mass
        std::unique_ptr<const Mass> mass;
        if (mSumRhoDepth.rows() == 0) {
            if (mMassSolid.rows() == 1) {
                // 1D mass
                mass = std::make_unique<const Mass1D>(mMassSolid(0));
            } else {
                // 3D mass
                mass = std::make_unique<const Mass3D>(mMassSolid);
            }
        } else {
            if (mMassSolid.rows() == 1 && mSumRhoDepth.rows() == 1 &&
                mNormalTop.rows() == 1) {
                // 1D mass with ocean load
                eigen::DCol2 nsz;
                nsz << mNormalTop(0, 0), mNormalTop(0, 2);
                const eigen::DCol2 &nrt = geodesy::sz2rtheta(nsz, false);
                mass = std::make_unique<const MassOceanLoad1D>
                (mMassSolid(0), mSumRhoDepth(0) * nrt(0), nrt(1));
            } else {
                // 3D mass with ocean load
                // mass of ocean
                eigen::DColX massOcean;
                const eigen::DColX &area = mNormalTop.rowwise().norm();
                op1D_3D::times(mSumRhoDepth, area, massOcean);
                // unit normal
                const eigen::DMatX3 &unitNormal =
                mNormalTop.array().colwise() / area.array();
                // mass
                mass = std::make_unique<const MassOceanLoad3D>
                (op1D_3D::to3D(mMassSolid, mNr), op1D_3D::to3D(massOcean, mNr),
                 op1D_3D::to3D(unitNormal, mNr));
            }
        }
        // point
        mSolidPoint = std::make_shared<SolidPoint>
        (mNr, mCoords.transpose(), mGlobalTag, mass, timeScheme);
        // release
        domain.addSolidPoint(mSolidPoint);
    }
    
    // fluid point
    if (mMassFluid.rows() > 0) {
        // mass
        std::unique_ptr<const Mass> mass;
        if (mMassFluid.rows() == 1) {
            // 1D mass
            mass = std::make_unique<const Mass1D>(mMassFluid(0));
        } else {
            // 3D mass
            mass = std::make_unique<const Mass3D>(mMassFluid);
        }
        // point
        mFluidPoint = std::make_shared<FluidPoint>
        (mNr, mCoords.transpose(), mGlobalTag, mass, timeScheme);
        // release
        domain.addFluidPoint(mFluidPoint);
    }
    
    // check empty
    if (!mSolidPoint && !mFluidPoint) {
        throw std::runtime_error("GLLPoint::release || "
                                 "Point is neither solid nor fluid.");
    }
    
    
    //////////////////////////// boundaries ////////////////////////////
    // solid-fluid
    if (mSolidPoint && mFluidPoint) {
        std::unique_ptr<const SolidFluidCoupling> sfc = nullptr;
        // if a rank contains pure fluid, mNormalSFU is unitialized
        // in such case, mNormalSFU = mNormalSFA
        if (mNormalSFU.rows() == 0) {
            mNormalSFU = mNormalSFA;
        }
        // 1D or 3D
        if (mNormalSFA.rows() == 1 && mMassFluid.rows() == 1) {
            // 1D coupling
            sfc = std::make_unique<SolidFluidCoupling1D>
            (mSolidPoint, mFluidPoint,
             mNormalSFU(0, 0), mNormalSFU(0, 2),
             mNormalSFA(0, 0), mNormalSFA(0, 2),
             mMassFluid(0));
        } else {
            // 3D coupling
            sfc = std::make_unique<SolidFluidCoupling3D>
            (mSolidPoint, mFluidPoint,
             op1D_3D::to3D(mNormalSFU, mNr), op1D_3D::to3D(mNormalSFA, mNr),
             op1D_3D::to3D(mMassFluid, mNr));
        }
        domain.getSolidFluidBoundary()->addSolidFluidCoupling(sfc);
    }
    
    // Clayton ABC
    for (auto itm = mClaytonABC.begin(); itm != mClaytonABC.end(); itm++) {
        for (auto itv = itm->second.begin(); itv != itm->second.end(); itv++) {
            bool fluid = std::get<0>(*itv);
            const eigen::DMatX3 &nABC = std::get<1>(*itv);
            const eigen::DColX &rhoVp = std::get<2>(*itv);
            const eigen::DColX &rhoVs = std::get<3>(*itv);
            // fluid
            if (fluid) {
                std::unique_ptr<const ClaytonFluid> abc = nullptr;
                if (nABC.rows() == 1 && rhoVp.rows() == 1) {
                    // 1D fluid
                    abc = std::make_unique<const ClaytonFluid1D>
                    (mFluidPoint, rhoVp(0), nABC.row(0).norm());
                } else {
                    // 3D fluid
                    abc = std::make_unique<const ClaytonFluid3D>
                    (mFluidPoint,
                     op1D_3D::to3D(nABC, mNr).rowwise().norm(),
                     op1D_3D::to3D(rhoVp, mNr));
                }
                domain.getAbsorbingBoundary()->addClaytonFluid(abc);
            } else {
                std::unique_ptr<const ClaytonSolid> abc = nullptr;
                if (nABC.rows() == 1 && rhoVp.rows() == 1 &&
                    rhoVs.rows() == 1) {
                    // 1D solid
                    eigen::DCol2 nsz;
                    nsz << nABC(0, 0), nABC(0, 2);
                    const eigen::DCol2 &nrt = geodesy::sz2rtheta(nsz, false);
                    abc = std::make_unique<const ClaytonSolid1D>
                    (mSolidPoint, rhoVp(0), rhoVs(0), nrt(0), nrt(1));
                } else {
                    // 3D solid
                    const eigen::DColX &area = nABC.rowwise().norm();
                    const eigen::DMatX3 &unitNormal =
                    nABC.array().colwise() / area.array();
                    abc = std::make_unique<const ClaytonSolid3D>
                    (mSolidPoint,
                     op1D_3D::to3D(rhoVp, mNr), op1D_3D::to3D(rhoVs, mNr),
                     op1D_3D::to3D(area, mNr), op1D_3D::to3D(unitNormal, mNr));
                }
                domain.getAbsorbingBoundary()->addClaytonSolid(abc);
            }
        }
    }
    
    // sponge ABC
    if (abc.sponge()) {
        // compute gamma based on the closest absorbing boundary
        double gammaMax = -1.;
        for (const std::string &key: abc.getBoundaryKeys()) {
            // data
            const auto &date = abc.getSpongeData(key);
            double outer = std::get<0>(date);
            double span = std::get<1>(date);
            double U0 = std::get<2>(date);
            // coord of me
            double coord = 0.;
            if (geodesy::isCartesian()) {
                coord = (key == "RIGHT" ? mCoords(0) : mCoords(1));
            } else {
                const eigen::DCol2 &rt = geodesy::sz2rtheta(mCoords, false);
                coord = (key == "RIGHT" ? rt(1) : rt(0));
            }
            // gamma
            double distToOuter = 1. / span * (outer - coord);
            if (distToOuter > 1.) {
                // point is inside the inner boundary, skip
                // there is no need to check distToOuter < 0.
                continue;
            }
            static const double piHalf = numerical::dPi / 2.;
            double gamma = U0 * pow(cos(piHalf * distToOuter), 2.);
            // use the max (closest)
            gammaMax = std::max(gamma, gammaMax);
        }
        // release
        if (gammaMax > numerical::dEpsilon) {
            if (mSolidPoint) {
                std::unique_ptr<const SpongeSolid> sponge =
                std::make_unique<const SpongeSolid>(mSolidPoint, gammaMax);
                domain.getAbsorbingBoundary()->addSpongeSolid(sponge);
            }
            if (mFluidPoint) {
                std::unique_ptr<const SpongeFluid> sponge =
                std::make_unique<const SpongeFluid>(mFluidPoint, gammaMax);
                domain.getAbsorbingBoundary()->addSpongeFluid(sponge);
            }
        }
    }
    
    // axial boundary
    if (mAxial) {
        if (mSolidPoint) {
            domain.getAxialBoundary()->addPoint(mSolidPoint);
        }
        if (mFluidPoint) {
            domain.getAxialBoundary()->addPoint(mFluidPoint);
        }
    }
    
    // fluid surface without ABC
    if (mFluidPoint && mSurface &&
        mClaytonABC.find("TOP") == mClaytonABC.end()) {
        domain.getFluidSurfaceBoundary()->addPoint(mFluidPoint);
    }
    
    // free dummy memory
    mMassFluid.resize(0);
    mMassSolid.resize(0);
    mNormalSFU.resize(0, 3);
    mNormalSFA.resize(0, 3);
    mClaytonABC.clear();
    mNormalTop.resize(0, 3);
    mSumRhoDepth.resize(0);
}
