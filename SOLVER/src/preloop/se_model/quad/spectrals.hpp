//
//  spectrals.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/18/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element constants

#ifndef spectrals_hpp
#define spectrals_hpp

#include "eigen_sem.hpp"
#include <vector>

namespace spectrals {
    ///////////////////////////// internal /////////////////////////////
    // copied from salvus
    namespace internal {
#if _NPOL == 1
        // nPol = 1
        const std::vector<double> pGLL = {
            -1.000000000000, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, +1.000000000000};
        const std::vector<double> wGLL = {
            +1.000000000000, +1.000000000000};
        const std::vector<double> wGLJ = {
            +0.666666666667, +1.333333333333};
        const std::vector<double> gGLL = {
            -0.500000000000, -0.500000000000, +0.500000000000, +0.500000000000};
        const std::vector<double> gGLJ = {
            -0.500000000000, -0.500000000000, +0.500000000000, +0.500000000000};
#elif _NPOL == 2
        // nPol = 2
        const std::vector<double> pGLL = {
            -1.000000000000, +0.000000000000, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, +0.200000000000, +1.000000000000};
        const std::vector<double> wGLL = {
            +0.333333333333, +1.333333333333, +0.333333333333};
        const std::vector<double> wGLJ = {
            +0.111111111111, +1.388888888889, +0.500000000000};
        const std::vector<double> gGLL = {
            -1.500000000000, -0.500000000000, +0.500000000000, +2.000000000000,
            +0.000000000000, -2.000000000000, -0.500000000000, +0.500000000000,
            +1.500000000000};
        const std::vector<double> gGLJ = {
            -1.333333333333, -0.333333333333, +0.333333333333, +2.083333333333,
            -0.416666666667, -2.083333333333, -0.750000000000, +0.750000000000,
            +1.750000000000};
#elif _NPOL == 3
        // nPol = 3
        const std::vector<double> pGLL = {
            -1.000000000000, -0.447213595500, +0.447213595500, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.261203874964, +0.546918160678, +1.000000000000};
        const std::vector<double> wGLL = {
            +0.166666666667, +0.833333333333, +0.833333333333, +0.166666666667};
        const std::vector<double> wGLJ = {
            +0.033333333333, +0.614297739604, +1.085702260396, +0.266666666667};
        const std::vector<double> gGLL = {
            -3.000000000000, -0.809016994375, +0.309016994375, -0.500000000000,
            +4.045084971875, +0.000000000000, -1.118033988750, +1.545084971875,
            -1.545084971875, +1.118033988750, +0.000000000000, -4.045084971875,
            +0.500000000000, -0.309016994375, +0.809016994375, +3.000000000000};
        const std::vector<double> gGLJ = {
            -2.500000000000, -0.445902906223, +0.160188620509, -0.250000000000,
            +4.108757210636, -0.676776695297, -0.930801637161, +1.203427124746,
            -2.608757210636, +1.645087351447, -0.323223304703, -4.453427124746,
            +1.000000000000, -0.522407749927, +1.093836321356, +3.500000000000};
#elif _NPOL == 4
        // nPol = 4
        const std::vector<double> pGLL = {
            -1.000000000000, -0.654653670708, +0.000000000000, +0.654653670708,
            +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.507787629558, +0.132300820777, +0.708820142114,
            +1.000000000000};
        const std::vector<double> wGLL = {
            +0.100000000000, +0.544444444444, +0.711111111111, +0.544444444444,
            +0.100000000000};
        const std::vector<double> wGLJ = {
            +0.013333333333, +0.289656694561, +0.736004369482, +0.794338935957,
            +0.166666666667};
        const std::vector<double> gGLL = {
            -5.000000000000, -1.240990253031, +0.375000000000, -0.259009746969,
            +0.500000000000, +6.756502488724, -0.000000000000, -1.336584577695,
            +0.763762615826, -1.410164177942, -2.666666666667, +1.745743121888,
            +0.000000000000, -1.745743121888, +2.666666666667, +1.410164177942,
            -0.763762615826, +1.336584577695, +0.000000000000, -6.756502488724,
            -0.500000000000, +0.259009746969, -0.375000000000, +1.240990253031,
            +5.000000000000};
        const std::vector<double> gGLJ = {
            -4.000000000000, -0.616438928116, +0.168105664494, -0.107222291934,
            +0.200000000000, +6.695837336887, -1.015821685975, -0.980080160442,
            +0.496350282495, -0.874333733232, -4.639743885094, +2.490338715010,
            -0.441578766725, -1.669642166429, +2.421846619237, +3.193906548208,
            -1.361164311623, +1.801975418767, -0.292599547300, -7.497512886006,
            -1.250000000000, +0.503086210704, -0.548422156095, +1.573113723168,
            +5.750000000000};
#elif _NPOL == 5
        // nPol = 5
        const std::vector<double> pGLL = {
            -1.000000000000, -0.765055323929, -0.285231516481, +0.285231516481,
            +0.765055323929, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.650778856692, -0.156370431808, +0.373489378736,
            +0.797296273400, +1.000000000000};
        const std::vector<double> wGLL = {
            +0.066666666667, +0.378474956298, +0.554858377035, +0.554858377035,
            +0.378474956298, +0.066666666667};
        const std::vector<double> wGLJ = {
            +0.006349206349, +0.150329375443, +0.452292685042, +0.685821572112,
            +0.590921446768, +0.114285714286};
        const std::vector<double> gGLL = {
            -7.500000000000, -1.786364948339, +0.484951047854, -0.269700610832,
            +0.237781177984, -0.500000000000, +10.14141593632, +0.000000000000,
            -1.721256952830, +0.786356672223, -0.653547507430, +1.349913314190,
            -4.036187270305, +2.523426777429, -0.000000000000, -1.752961966368,
            +1.152828158536, -2.244684648176, +2.244684648176, -1.152828158536,
            +1.752961966368, +0.000000000000, -2.523426777429, +4.036187270305,
            -1.349913314190, +0.653547507430, -0.786356672223, +1.721256952830,
            -0.000000000000, -10.14141593632, +0.500000000000, -0.237781177984,
            +0.269700610832, -0.484951047854, +1.786364948339, +7.500000000000};
        const std::vector<double> gGLJ = {
            -5.833333333333, -0.832247046049, +0.198615502767, -0.099070383168,
            +0.081562432334, -0.166666666667, +9.852505318461, -1.431757525514,
            -1.166074456360, +0.457090980013, -0.348310009836, +0.694763597562,
            -7.074296699204, +3.508342566249, -0.592677187775, -1.532650853411,
            +0.917377924310, -1.720350072793, +5.350637717282, -2.085306704601,
            +2.323992964175, -0.364036306171, -2.541983204845, +3.910039630469,
            -3.795513003205, +1.369152597947, -1.198556395436, +2.190237890068,
            -0.278195647206, -11.21778648857, +1.500000000000, -0.528183888031,
            +0.434699572628, -0.651571327331, +2.169548505244, +8.500000000000};
#elif _NPOL == 6
        // nPol = 6
        const std::vector<double> pGLL = {
            -1.000000000000, -0.830223896279, -0.468848793471, +0.000000000000,
            +0.468848793471, +0.830223896279, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.740123648580, -0.353852634128, +0.098902793151,
            +0.528842304451, +0.850846569722, +1.000000000000};
        const std::vector<double> wGLL = {
            +0.047619047619, +0.276826047362, +0.431745381210, +0.487619047619,
            +0.431745381210, +0.276826047362, +0.047619047619};
        const std::vector<double> wGLJ = {
            +0.003401360544, +0.084736552956, +0.280303211939, +0.501646961918,
            +0.594575445079, +0.452003134230, +0.083333333333};
        const std::vector<double> gGLL = {
            -10.50000000000, -2.442926014244, +0.625256665515, -0.312500000000,
            +0.226099400943, -0.226611870395, +0.500000000000, +14.20157660292,
            +0.000000000000, -2.215804283170, +0.907544471269, -0.616390835518,
            +0.602247179636, -1.317373435702, -5.668985225546, +3.455828214294,
            -0.000000000000, -2.006969240589, +1.066441904006, -0.961339797289,
            +2.049964813077, +3.200000000000, -1.598606688098, +2.266698087086,
            +0.000000000000, -2.266698087086, +1.598606688098, -3.200000000000,
            -2.049964813077, +0.961339797289, -1.066441904006, +2.006969240589,
            +0.000000000000, -3.455828214294, +5.668985225546, +1.317373435702,
            -0.602247179636, +0.616390835518, -0.907544471269, +2.215804283170,
            -0.000000000000, -14.20157660292, -0.500000000000, +0.226611870395,
            -0.226099400943, +0.312500000000, -0.625256665515, +2.442926014244,
            +10.50000000000};
        const std::vector<double> gGLJ = {
            -8.000000000000, -1.090282836280, +0.241099171111, -0.105970006790,
            +0.069964084602, -0.066282700029, +0.142857142857, +13.58086096613,
            -1.923991918725, -1.423407465389, +0.489847228970, -0.297496555578,
            +0.272146418542, -0.579489952568, -9.934388192566, +4.708542777910,
            -0.773817284430, -1.651015988227, +0.777857271666, -0.653679237176,
            +1.354669228574, +7.814451198196, -2.899933566384, +2.954754420560,
            -0.454999298497, -2.136427944510, +1.401015695890, -2.722814612918,
            -6.115042231034, +2.087459788858, -1.649980499013, +2.532194336691,
            -0.327044848605, -3.561813185444, +5.669278851179, +4.404118259271,
            -1.451687965340, +1.054090896573, -1.262368824535, +2.707731603596,
            -0.270146649744, -15.61450065713, -1.750000000000, +0.569893719960,
            -0.402739239413, +0.452312552389, -0.794583611171, +2.878759657960,
            +11.75000000000};
#elif _NPOL == 7
        // nPol = 7
        const std::vector<double> pGLL = {
            -1.000000000000, -0.871740148510, -0.591700181433, -0.209299217902,
            +0.209299217902, +0.591700181433, +0.871740148510, +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.799381854549, -0.491905791339, -0.111733935400,
            +0.283539772448, +0.633793327056, +0.885688481784, +1.000000000000};
        const std::vector<double> wGLL = {
            +0.035714285714, +0.210704227144, +0.341122692484, +0.412458794659,
            +0.412458794659, +0.341122692484, +0.210704227144, +0.035714285714};
        const std::vector<double> wGLJ = {
            +0.001984126984, +0.051068915198, +0.179218780464, +0.353399617690,
            +0.490974910515, +0.504783970569, +0.355077615089, +0.063492063492};
        const std::vector<double> gGLL = {
            -14.00000000000, -3.209915703003, +0.792476681321, -0.372150435729,
            +0.243330712724, -0.203284568901, +0.219957514771, -0.500000000000,
            +18.93759860712, -0.000000000000, -2.806475794736, +1.078944688790,
            -0.661157350900, +0.537039586158, -0.573565414940, +1.297687388320,
            -7.569289819348, +4.543585064567, +0.000000000000, -2.378187233516,
            +1.135358016881, -0.845022556507, +0.869448098331, -1.941659425544,
            +4.297908164265, -2.112061214315, +2.875517405973, +0.000000000000,
            -2.388924359158, +1.372785831806, -1.294232050914, +2.810188989258,
            -2.810188989258, +1.294232050914, -1.372785831806, +2.388924359158,
            -0.000000000000, -2.875517405973, +2.112061214315, -4.297908164265,
            +1.941659425544, -0.869448098331, +0.845022556507, -1.135358016881,
            +2.378187233516, +0.000000000000, -4.543585064567, +7.569289819348,
            -1.297687388320, +0.573565414940, -0.537039586158, +0.661157350900,
            -1.078944688790, +2.806475794736, +0.000000000000, -18.93759860712,
            +0.500000000000, -0.219957514771, +0.203284568901, -0.243330712724,
            +0.372150435729, -0.792476681321, +3.209915703003, +14.00000000000};
        const std::vector<double> gGLJ = {
            -10.50000000000, -1.389476486679, +0.292862473417, -0.119295379021,
            +0.070042356542, -0.054268751950, +0.056061941539, -0.125000000000,
            +17.88168233075, -2.492296989769, -1.736102032676, +0.552814143253,
            -0.297818477729, +0.221935236223, -0.225060166491, +0.498419768109,
            -13.22658674299, +6.092592487024, -0.984069472702, -1.873176748170,
            +0.779131554161, -0.529317990644, +0.515713835870, -1.126135471067,
            +10.62405321724, -3.825503363859, +3.693697418072, -0.562894407347,
            -2.146374303081, +1.122320073345, -1.000212472296, +2.122132292829,
            -8.666038013272, +2.863217279736, -2.134452896586, +2.981938515596,
            -0.389547726321, -2.815751477198, +1.952829121665, -3.881305630039,
            +6.903287013962, -2.193689631258, +1.490866282647, -1.603083745748,
            +2.894946728143, -0.306036260352, -4.733379404537, +7.699575875750,
            -5.016397805686, +1.564823275746, -1.021759206458, +1.004961639649,
            -1.412304157192, +3.329576944332, -0.265155143509, -20.68768683558,
            +2.000000000000, -0.619666570940, +0.398957434286, -0.381264018212,
            +0.501924025479, -0.968457773757, +3.699202287759, +15.50000000000};
#elif _NPOL == 8
        // nPol = 8
        const std::vector<double> pGLL = {
            -1.000000000000, -0.899757995411, -0.677186279511, -0.363117463826,
            +0.000000000000, +0.363117463826, +0.677186279511, +0.899757995411,
            +1.000000000000};
        const std::vector<double> pGLJ = {
            -1.000000000000, -0.840590144594, -0.591181348791, -0.271893788434,
            +0.078983005512, +0.419173799923, +0.707657051057, +0.909616131210,
            +1.000000000000};
        const std::vector<double> wGLL = {
            +0.027777777778, +0.165495361561, +0.274538712500, +0.346428510973,
            +0.371519274376, +0.346428510973, +0.274538712500, +0.165495361561,
            +0.027777777778};
        const std::vector<double> wGLJ = {
            +0.001234567901, +0.032489354604, +0.118611429028, +0.248991726162,
            +0.380501475574, +0.455288537175, +0.427395174812, +0.285487734745,
            +0.050000000000};
        const std::vector<double> gGLL = {
            -18.00000000000, -4.087013702034, +0.985360090075, -0.444613449281,
            +0.273437500000, -0.207734512036, +0.189655591978, -0.215654018702,
            +0.500000000000, +24.34974517159, +0.000000000000, -3.488358753434,
            +1.287960750064, -0.741782397916, +0.547300160534, -0.492350938316,
            +0.555704981284, -1.284830632700, -9.738701657212, +5.786805816637,
            -0.000000000000, -2.834458912079, +1.269413086358, -0.855726185093,
            +0.738349277190, -0.816756381741, +1.874440873447, +5.544963906949,
            -2.696065440314, +3.576680940126, +0.000000000000, -2.659310217574,
            +1.376964893761, -1.079803811283, +1.145653738455, -2.590745676559,
            -3.657142857143, +1.665221645005, -1.717832157195, +2.851915968463,
            +0.000000000000, -2.851915968463, +1.717832157195, -1.665221645005,
            +3.657142857143, +2.590745676559, -1.145653738455, +1.079803811283,
            -1.376964893761, +2.659310217574, -0.000000000000, -3.576680940126,
            +2.696065440314, -5.544963906949, -1.874440873447, +0.816756381741,
            -0.738349277190, +0.855726185093, -1.269413086358, +2.834458912079,
            +0.000000000000, -5.786805816637, +9.738701657212, +1.284830632700,
            -0.555704981284, +0.492350938316, -0.547300160534, +0.741782397916,
            -1.287960750064, +3.488358753434, -0.000000000000, -24.34974517159,
            -0.500000000000, +0.215654018702, -0.189655591978, +0.207734512036,
            -0.273437500000, +0.444613449281, -0.985360090075, +4.087013702034,
            +18.00000000000};
        const std::vector<double> gGLJ = {
            -13.33333333333, -1.729365450870, +0.352921826495, -0.136768251355,
            +0.074658577581, -0.051891179862, +0.044509941686, -0.048700417732,
            +0.111111111111, +22.75531678581, -3.136568932487, -2.098434933196,
            +0.635181204545, -0.317765019655, +0.212049869585, -0.178080117148,
            +0.192746978165, -0.437953949982, -16.95352768113, +7.660920605690,
            -1.223036176363, -2.161667099532, +0.833112733327, -0.505179655154,
            +0.405595699733, -0.429484455670, +0.967962622390, +13.79193601055,
            -4.867898007680, +4.537819221242, -0.686712998814, -2.305470870321,
            +1.070110376369, -0.779203096012, +0.790426256874, -1.754513076857,
            -11.50511806819, +3.721528492641, -2.672597632026, +3.523149469949,
            -0.463399328299, -2.687275899034, +1.500851881983, -1.389873256750,
            +2.995197928248, +9.568311045531, -2.971554101632, +1.939126002370,
            -1.956727620561, +3.215456421583, -0.352317665410, -3.577733165687,
            +2.574907798813, -5.195318264270, -7.704450394678, +2.342631416589,
            -1.461491918712, +1.337504858365, -1.685819618702, +3.358542477854,
            -0.292798837852, -6.058396544276, +10.00085617749, +5.630865635436,
            -1.693690097776, +1.033732966230, -0.906283132521, +1.042812690946,
            -1.614590604943, +4.046835358801, -0.261832727441, -26.43734254813,
            -2.250000000000, +0.673996075525, -0.408039356041, +0.352323569924,
            -0.393585586459, +0.570552280595, -1.169977665504, +4.630206368017,
            +19.75000000000};
#endif
    }
    
    
    ////////////////// global constants //////////////////
    // point positions
    // GLL
    const eigen::DColP gPositionsGLL =
    Eigen::Map<const eigen::DColP>(internal::pGLL.data());
    // GLJ
    const eigen::DColP gPositionsGLJ =
    Eigen::Map<const eigen::DColP>(internal::pGLJ.data());
    
    // weights
    // GLL
    const eigen::DColP gWeightsGLL =
    Eigen::Map<const eigen::DColP>(internal::wGLL.data());
    // GLJ
    const eigen::DColP gWeightsGLJ =
    Eigen::Map<const eigen::DColP>(internal::wGLJ.data());
    
    // G matrix
    // GLL
    const eigen::DMatPP_RM gGMatrixGLL =
    Eigen::Map<const eigen::DMatPP_RM>(internal::gGLL.data());
    // GLJ
    const eigen::DMatPP_RM gGMatrixGLJ =
    Eigen::Map<const eigen::DMatPP_RM>(internal::gGLJ.data());
    
    
    ///////////////////////////// methods /////////////////////////////
    // xieta on element
    inline const eigen::DMat2N &getXiEtaElement(bool axial) {
        static eigen::DMat2N xietaAX, xietaNA;
        static bool formed = false;
        if (!formed) {
            for (int ipol = 0; ipol < spectral::nPED; ipol++) {
                for (int jpol = 0; jpol < spectral::nPED; jpol++) {
                    int ipnt = ipol * spectral::nPED + jpol;
                    // axial
                    xietaAX(0, ipnt) = gPositionsGLJ(ipol);
                    xietaAX(1, ipnt) = gPositionsGLL(jpol);
                    // non-axial
                    xietaNA(0, ipnt) = gPositionsGLL(ipol);
                    xietaNA(1, ipnt) = gPositionsGLL(jpol);
                }
            }
            formed = true;
        }
        return axial ? xietaAX : xietaNA;
    }
    
    // weights on element
    inline eigen::DRowN getWeightsElement(bool axial) {
        if (axial) {
            eigen::DMatPP_RM wmat = (gWeightsGLJ * gWeightsGLL.transpose());
            return Eigen::Map<eigen::DRowN>(wmat.data());
        } else {
            eigen::DMatPP_RM wmat = (gWeightsGLL * gWeightsGLL.transpose());
            return Eigen::Map<eigen::DRowN>(wmat.data());
        }
    }
    
    // interpolate a nodal field to a GLL field
    template <int ROW>
    Eigen::Matrix<double, ROW, spectral::nPEM>
    interpolateGLL(const Eigen::Matrix<double, ROW, 4> &nodal, bool axial) {
        // form linear shape functions
        static Eigen::Matrix<double, 4, spectral::nPEM> lsfAX, lsfNA;
        static bool formed = false;
        if (!formed) {
            {
                // axial
                const eigen::DMat2N &xieta = getXiEtaElement(true);
                const auto &xip = 1. + xieta.row(0).array();
                const auto &xim = 1. - xieta.row(0).array();
                const auto &etp = 1. + xieta.row(1).array();
                const auto &etm = 1. - xieta.row(1).array();
                lsfAX.row(0) = xim.cwiseProduct(etm);
                lsfAX.row(1) = xip.cwiseProduct(etm);
                lsfAX.row(2) = xip.cwiseProduct(etp);
                lsfAX.row(3) = xim.cwiseProduct(etp);
                lsfAX /= 4.;
            }
            {
                // non-axial
                const eigen::DMat2N &xieta = getXiEtaElement(false);
                const auto &xip = 1. + xieta.row(0).array();
                const auto &xim = 1. - xieta.row(0).array();
                const auto &etp = 1. + xieta.row(1).array();
                const auto &etm = 1. - xieta.row(1).array();
                lsfNA.row(0) = xim.cwiseProduct(etm);
                lsfNA.row(1) = xip.cwiseProduct(etm);
                lsfNA.row(2) = xip.cwiseProduct(etp);
                lsfNA.row(3) = xim.cwiseProduct(etp);
                lsfNA /= 4.;
            }
            formed = true;
        }
        return axial ? nodal * lsfAX : nodal * lsfNA;
    }
}

#endif /* spectrals_hpp */
