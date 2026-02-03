/*
 * modelsolver01-06.cpp
 * 文件作用: 压裂水平井复合页岩油模型核心计算类实现
 * 修改记录:
 * 1. [多裂缝支持] 修改 PWD_composite 逻辑，支持计算多条横向裂缝(Transverse Fractures)。
 * 2. [参数重构] "nf" 现在代表裂缝条数，原微元段数参数重命名为 "n_seg"。
 * 3. [几何计算] 引入 spacingD 参数，微元距离计算升级为二维欧氏距离。
 * 4. [产量约束] 矩阵方程约束条件保持为 Sum(q_all) = 1/z，即输入 q 代表总产量。
 */

#include "modelsolver01-06.h"
#include "pressurederivativecalculator.h"

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <algorithm>
#include <QDebug>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 简单二维点结构体
struct Point2D {
    double x;
    double y;
};

// ---------------------- 辅助数学函数 ----------------------

static double safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    return boost::math::cyl_bessel_k(v, x);
}

static double safe_bessel_i_scaled(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x);
    try {
        return boost::math::cyl_bessel_i(v, x) * std::exp(-x);
    } catch (...) {
        return 0.0;
    }
}

// ---------------------- 类实现 ----------------------

ModelSolver01_06::ModelSolver01_06(ModelType type)
    : m_type(type)
    , m_highPrecision(true)
{
}

ModelSolver01_06::~ModelSolver01_06()
{
}

void ModelSolver01_06::setHighPrecision(bool high)
{
    m_highPrecision = high;
}

QString ModelSolver01_06::getModelName(ModelType type)
{
    switch(type) {
    case Model_1: return "模型1: 变井储+无限大边界";
    case Model_2: return "模型2: 恒定井储+无限大边界";
    case Model_3: return "模型3: 变井储+封闭边界";
    case Model_4: return "模型4: 恒定井储+封闭边界";
    case Model_5: return "模型5: 变井储+定压边界";
    case Model_6: return "模型6: 恒定井储+定压边界";
    default: return "未知模型";
    }
}

QVector<double> ModelSolver01_06::generateLogTimeSteps(int count, double startExp, double endExp)
{
    QVector<double> t;
    if (count <= 0) return t;
    t.reserve(count);
    for (int i = 0; i < count; ++i) {
        double exponent = startExp + (endExp - startExp) * i / (count - 1);
        t.append(pow(10.0, exponent));
    }
    return t;
}

ModelCurveData ModelSolver01_06::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    QVector<double> tPoints = providedTime;
    if (tPoints.isEmpty()) {
        tPoints = generateLogTimeSteps(100, -3.0, 3.0);
    }

    double phi = params.value("phi", 0.05);
    double mu = params.value("mu", 0.5);
    double B = params.value("B", 1.05);
    double Ct = params.value("Ct", 5e-4);
    double q = params.value("q", 5.0);
    double h = params.value("h", 20.0);
    double kf = params.value("kf", 1e-3);
    double L = params.value("L", 1000.0);
    if (L < 1e-9) L = 1000.0;

    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tPoints, QVector<double>(tPoints.size(), 0.0), QVector<double>(tPoints.size(), 0.0));
    }

    // tD = 14.4 * kf * t / (phi * mu * Ct * L^2)
    double td_coeff = 14.4 * kf / (phi * mu * Ct * pow(L, 2));

    QVector<double> tD_vec;
    tD_vec.reserve(tPoints.size());
    for(double t : tPoints) {
        tD_vec.append(td_coeff * t);
    }

    QVector<double> PD_vec, Deriv_vec;
    auto func = std::bind(&ModelSolver01_06::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);

    QMap<QString, double> calcParams = params;

    // 默认 Stehfest N=10
    if (!calcParams.contains("N") || calcParams["N"] < 4) {
        calcParams["N"] = 10;
    }

    // [参数默认值设置]
    // 裂缝条数 nf (默认1)
    if (!calcParams.contains("nf") || calcParams["nf"] < 1) {
        calcParams["nf"] = 1;
    }
    // 单缝段数 n_seg (原 nf, 默认3)
    if (!calcParams.contains("n_seg")) {
        calcParams["n_seg"] = 10;
    }
    if (calcParams["n_seg"] < 1) calcParams["n_seg"] = 1;

    calculatePDandDeriv(tD_vec, calcParams, func, PD_vec, Deriv_vec);

    // 转换为有因次压力: Dp = 1.842e-3*q*mu*B*PD/(kf*h);
    // 注意: q 是总产量，PD 是基于总产量归一化计算的无因次压力
    double p_coeff = 1.842e-3 * q * mu * B / (kf * h);

    QVector<double> finalP(tPoints.size()), finalDP(tPoints.size());
    for(int i=0; i<tPoints.size(); ++i) {
        finalP[i] = p_coeff * PD_vec[i];
        finalDP[i] = p_coeff * Deriv_vec[i];
    }

    return std::make_tuple(tPoints, finalP, finalDP);
}

void ModelSolver01_06::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                                           std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                                           QVector<double>& outPD, QVector<double>& outDeriv)
{
    int numPoints = tD.size();
    outPD.resize(numPoints);
    outDeriv.resize(numPoints);

    int N = (int)params.value("N", 10);
    if (N > 18) N = 18;
    if (N % 2 != 0) N = 10;
    double ln2 = log(2.0);
    double gamaD = params.value("gamaD", 0.0);

    for (int k = 0; k < numPoints; ++k) {
        double t = tD[k];
        if (t <= 1e-10) { outPD[k] = 0.0; continue; }

        double pd_val = 0.0;
        for (int m = 1; m <= N; ++m) {
            double z = m * ln2 / t;
            double pf = laplaceFunc(z, params);
            if (std::isnan(pf) || std::isinf(pf)) pf = 0.0;
            pd_val += stefestCoefficient(m, N) * pf;
        }
        outPD[k] = pd_val * ln2 / t;

        if (std::abs(gamaD) > 1e-9) {
            double arg = 1.0 - gamaD * outPD[k];
            if (arg > 1e-12) {
                outPD[k] = -1.0 / gamaD * std::log(arg);
            }
        }
    }

    if (numPoints > 2) {
        outDeriv = PressureDerivativeCalculator::calculateBourdetDerivative(tD, outPD, 0.1);
    } else {
        outDeriv.fill(0.0);
    }
}

// 核心 Laplace 函数
double ModelSolver01_06::flaplace_composite(double z, const QMap<QString, double>& p) {
    double M12 = 1.0;
    if (p.contains("M12")) {
        M12 = p.value("M12");
    } else {
        double kf = p.value("kf", 1.0);
        double km = p.value("km", 0.01);
        if(km < 1e-12) km = 1e-12;
        M12 = kf / km;
    }

    double L = p.value("L", 1000.0);
    double Lf = p.value("Lf", 100.0);
    double rm = p.value("rm", 500.0);
    double re = p.value("re", 20000.0);

    double LfD = (L > 1e-9) ? Lf / L : 0.1;
    double rmD = (L > 1e-9) ? rm / L : 0.5;
    double reD = (L > 1e-9) ? re / L : 20.0;

    double omga1 = p.value("omega1", 0.4);
    double omga2 = p.value("omega2", 0.08);
    double remda1 = p.contains("lambda1") ? p.value("lambda1") : p.value("remda1", 1e-3);
    double remda2 = p.contains("lambda2") ? p.value("lambda2") : p.value("remda2", 1e-4);
    double eta12 = p.value("eta", 0.2);
    if (p.contains("eta12")) eta12 = p.value("eta12");

    // [参数变更]
    // nf: 裂缝条数 (Input)
    // n_seg: 单条裂缝段数 (Input, 原 nf)
    // frac_spacing: 裂缝间距 (Input)
    int n_fracs = (int)p.value("nf", 1);
    int n_seg = (int)p.value("n_seg", 10);
    double spacing = p.value("frac_spacing", 100.0); // 默认间距
    double spacingD = (L > 1e-9) ? spacing / L : 0.1;

    // 安全保护
    if (n_fracs < 1) n_fracs = 1;
    if (n_seg < 1) n_seg = 1;

    double one_minus_omega1 = 1.0 - omga1;
    double den_fs1 = one_minus_omega1 * z + remda1;
    double fs1 = 1.0;
    if (std::abs(den_fs1) > 1e-20) {
        fs1 = (omga1 * one_minus_omega1 * z + remda1) / den_fs1;
    }

    double one_minus_omega2 = 1.0 - omga2;
    double den_fs2 = one_minus_omega2 * eta12 * z + remda2;
    double fs2 = 0.0;
    if (std::abs(den_fs2) > 1e-20) {
        fs2 = eta12 * (omga2 * one_minus_omega2 * eta12 * z + remda2) / den_fs2;
    }

    if (z * fs1 < 0) fs1 = 0;
    if (z * fs2 < 0) fs2 = 0;

    // 调用多裂缝 PWD
    double pf = PWD_composite(z, fs1, fs2, M12, LfD, rmD, reD, n_seg, n_fracs, spacingD, m_type);

    // 井储和表皮
    bool hasStorage = (m_type == Model_1 || m_type == Model_3 || m_type == Model_5);
    if (hasStorage) {
        double CD = p.value("cD", 0.0);
        double S = p.value("S", 0.0);
        if (CD > 1e-12 || std::abs(S) > 1e-12) {
            double num = z * pf + S;
            double den = z + CD * z * z * num;
            if (std::abs(den) > 1e-100) {
                pf = num / den;
            }
        }
    }

    return pf;
}

double ModelSolver01_06::PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                                       int n_seg, int n_fracs, double spacingD, ModelType type) {
    // -------------------------------------------------------------
    // [多裂缝离散化逻辑]
    // 1. 假设所有裂缝为横向裂缝(Transverse)，垂直于井筒(X轴)。
    // 2. 裂缝沿 X 轴分布，间距为 spacingD，对称中心在 X=0。
    // 3. 裂缝长度方向为 Y 轴，范围 [-LfD, LfD]。
    // -------------------------------------------------------------

    int total_segments = n_fracs * n_seg;
    double segLen = 2.0 * LfD / n_seg; // 单个微元长度 (Y方向)

    QVector<Point2D> segmentCenters;
    segmentCenters.reserve(total_segments);

    // 计算第一条裂缝的 X 坐标 (使裂缝簇关于原点对称)
    // X_k = startX + k * spacingD
    double startX = -(n_fracs - 1) * spacingD / 2.0;

    for (int k = 0; k < n_fracs; ++k) {
        double currentX = startX + k * spacingD;

        // 生成当前裂缝的 Y 方向微元
        for (int i = 0; i < n_seg; ++i) {
            double currentY = -LfD + (i + 0.5) * segLen;
            segmentCenters.append({currentX, currentY});
        }
    }

    // -------------------------------------------------------------

    double gama1 = sqrt(z * fs1);
    double gama2 = sqrt(z * fs2);
    double arg_g1_rm = gama1 * rmD;
    double arg_g2_rm = gama2 * rmD;

    double term_mAB_i0 = 0.0;
    double term_mAB_i1 = 0.0;

    double k0_g2_rm = safe_bessel_k(0, arg_g2_rm);
    double k1_g2_rm = safe_bessel_k(1, arg_g2_rm);
    double k0_g1_rm = safe_bessel_k(0, arg_g1_rm);
    double k1_g1_rm = safe_bessel_k(1, arg_g1_rm);

    bool isInfinite = (type == Model_1 || type == Model_2);
    bool isClosed = (type == Model_3 || type == Model_4);
    bool isConstP = (type == Model_5 || type == Model_6);

    // 边界条件处理 (假设边界是以井簇中心为圆心，半径 reD)
    if (!isInfinite && reD > 1e-5) {
        double arg_re = gama2 * reD;
        double i0_re_s = safe_bessel_i_scaled(0, arg_re);
        double i1_re_s = safe_bessel_i_scaled(1, arg_re);
        double k0_re = safe_bessel_k(0, arg_re);
        double k1_re = safe_bessel_k(1, arg_re);

        double i0_g2_rm_s = safe_bessel_i_scaled(0, arg_g2_rm);
        double i1_g2_rm_s = safe_bessel_i_scaled(1, arg_g2_rm);

        double exp_factor = 0.0;
        if ((arg_g2_rm - arg_re) > -700.0) {
            exp_factor = std::exp(arg_g2_rm - arg_re);
        }

        if (isClosed && i1_re_s > 1e-100) {
            term_mAB_i0 = (k1_re / i1_re_s) * i0_g2_rm_s * exp_factor;
            term_mAB_i1 = (k1_re / i1_re_s) * i1_g2_rm_s * exp_factor;
        } else if (isConstP && i0_re_s > 1e-100) {
            term_mAB_i0 = -(k0_re / i0_re_s) * i0_g2_rm_s * exp_factor;
            term_mAB_i1 = -(k0_re / i0_re_s) * i1_g2_rm_s * exp_factor;
        }
    }

    double term1 = term_mAB_i0 + k0_g2_rm;
    double term2 = term_mAB_i1 - k1_g2_rm;
    double Acup = M12 * gama1 * k1_g1_rm * term1 + gama2 * k0_g1_rm * term2;

    double i1_g1_rm_s = safe_bessel_i_scaled(1, arg_g1_rm);
    double i0_g1_rm_s = safe_bessel_i_scaled(0, arg_g1_rm);
    double Acdown_scaled = M12 * gama1 * i1_g1_rm_s * term1 - gama2 * i0_g1_rm_s * term2;

    if (std::abs(Acdown_scaled) < 1e-100) Acdown_scaled = (Acdown_scaled >= 0 ? 1e-100 : -1e-100);
    double Ac_prefactor = Acup / Acdown_scaled;

    // 建立矩阵
    int size = total_segments + 1;
    Eigen::MatrixXd A_mat(size, size);
    Eigen::VectorXd b_vec(size);
    b_vec.setZero();
    b_vec(total_segments) = 1.0; // 流量约束右端项 (Sum(q) = 1)

    for (int i = 0; i < total_segments; ++i) {
        for (int j = 0; j < total_segments; ++j) {
            Point2D pi = segmentCenters[i];
            Point2D pj = segmentCenters[j];

            // 积分变量 a 是沿裂缝方向 (Y轴) 的局部坐标，范围 [-segLen/2, segLen/2]
            // 被积距离 dist = sqrt( (xi - xj)^2 + (yi - (yj + a))^2 )

            double dx_sq = (pi.x - pj.x) * (pi.x - pj.x); // X方向距离平方 (常量)

            auto integrand = [&](double a) -> double {
                double dy = pi.y - (pj.y + a);
                double dist_val = std::sqrt(dx_sq + dy * dy); // 二维距离
                double arg_dist = gama1 * dist_val;

                double k0_val = safe_bessel_k(0, arg_dist);
                double term2_val = 0.0;
                double exponent = arg_dist - arg_g1_rm;

                if (exponent > -700.0) {
                    term2_val = Ac_prefactor * safe_bessel_i_scaled(0, arg_dist) * std::exp(exponent);
                }
                return k0_val + term2_val;
            };

            double val = 0.0;
            double halfLen = segLen / 2.0;

            // 自感应项 (同一微元 i==j): 奇异点在 a=0
            if (i == j) {
                val = 2.0 * adaptiveGauss(integrand, 0.0, halfLen, 1e-6, 0, 8);
            }
            // 同一裂缝的不同微元 (X相同，Y不同): 距离平滑，但接近
            else if (std::abs(pi.x - pj.x) < 1e-9) {
                val = adaptiveGauss(integrand, -halfLen, halfLen, 1e-6, 0, 5);
            }
            // 不同裂缝: 距离远，平滑
            else {
                val = adaptiveGauss(integrand, -halfLen, halfLen, 1e-6, 0, 3);
            }

            A_mat(i, j) = val / (M12 * 2.0 * LfD); // 这里的系数 2*LfD 是裂缝总长，需确认无因次定义
            // 修正：如果是多裂缝，单个裂缝长是 2*LfD?
            // 假设 LfD 是半长，无因次化通常针对单缝。此处保持原逻辑，但需注意分母物理意义。
            // 边界元系数通常为 1/(2*pi*T) * ...
            // 假设 MATLAB 原型中分母是几何因子。此处保持一致。
        }
    }

    // 流量约束方程: sum(q_i) = 1/z  => sum(q_i * z) = 1
    // 这意味着所有裂缝的总产量之和恒定为 q (无因次为1)
    for (int i = 0; i < total_segments; ++i) {
        A_mat(i, total_segments) = -1.0;
        A_mat(total_segments, i) = z;
    }
    A_mat(total_segments, total_segments) = 0.0;

    // 求解得到井底压力
    Eigen::VectorXd x_sol = A_mat.fullPivLu().solve(b_vec);

    // 返回解向量的最后一个元素，即 PWD
    return x_sol(total_segments);
}

double ModelSolver01_06::scaled_besseli(int v, double x) {
    return safe_bessel_i_scaled(v, x);
}

double ModelSolver01_06::gauss15(std::function<double(double)> f, double a, double b) {
    static const double X[] = { 0.0, 0.201194, 0.394151, 0.570972, 0.724418, 0.848207, 0.937299, 0.987993 };
    static const double W[] = { 0.202578, 0.198431, 0.186161, 0.166269, 0.139571, 0.107159, 0.070366, 0.030753 };
    double h = 0.5 * (b - a); double c = 0.5 * (a + b); double s = W[0] * f(c);
    for (int i = 1; i < 8; ++i) { double dx = h * X[i]; s += W[i] * (f(c - dx) + f(c + dx)); }
    return s * h;
}

double ModelSolver01_06::adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth) {
    double c = (a + b) / 2.0; double v1 = gauss15(f, a, b); double v2 = gauss15(f, a, c) + gauss15(f, c, b);
    if (depth >= maxDepth || std::abs(v1 - v2) < eps * (std::abs(v2) + 1.0)) return v2;
    return adaptiveGauss(f, a, c, eps/2, depth+1, maxDepth) + adaptiveGauss(f, c, b, eps/2, depth+1, maxDepth);
}

double ModelSolver01_06::stefestCoefficient(int i, int N) {
    double s = 0.0; int k1 = (i + 1) / 2; int k2 = std::min(i, N / 2);
    for (int k = k1; k <= k2; ++k) {
        double num = pow(k, N / 2.0) * factorial(2 * k);
        double den = factorial(N / 2 - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i);
        if(den!=0) s += num/den;
    }
    return ((i + N / 2) % 2 == 0 ? 1.0 : -1.0) * s;
}

double ModelSolver01_06::factorial(int n) {
    if(n<=1)return 1;
    double r=1;
    for(int i=2;i<=n;++i) r*=i;
    return r;
}
