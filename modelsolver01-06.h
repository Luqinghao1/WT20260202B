/*
 * modelsolver01-06.h
 * 文件作用: 压裂水平井复合页岩油模型核心计算类头文件
 * 功能描述:
 * 1. 定义模型类型枚举 (ModelType) 和曲线数据类型 (ModelCurveData)。
 * 2. 声明纯数学计算逻辑，包括拉普拉斯变换、贝塞尔函数计算、Stehfest 数值反演等。
 * 3. 实现了计算结果的容器定义，不依赖任何 UI 控件，仅负责数据输入与结果输出。
 * 4. 提供6种理论模型的解算接口，算法逻辑已根据MATLAB原型（modelwidget1A-6A）进行严格对齐。
 */

#ifndef MODELSOLVER01_06_H
#define MODELSOLVER01_06_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>

// 类型定义: <时间序列, 压力序列, 导数序列>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver01_06
{
public:
    // 模型类型枚举：对应不同的边界条件和井储模式
    // 对应关系:
    // Model 1/2 -> 无限大边界 (Model 1 变井储/有井储, Model 2 恒定井储/无井储)
    // Model 3/4 -> 封闭边界
    // Model 5/6 -> 定压边界
    enum ModelType {
        Model_1 = 0, // 模型1: 变井储 + 无限大边界 (对应 MATLAB modelwidget1A)
        Model_2,     // 模型2: 恒定井储 + 无限大边界 (对应 MATLAB modelwidget2A)
        Model_3,     // 模型3: 变井储 + 封闭边界
        Model_4,     // 模型4: 恒定井储 + 封闭边界
        Model_5,     // 模型5: 变井储 + 定压边界
        Model_6      // 模型6: 恒定井储 + 定压边界
    };

    // 构造函数：初始化模型类型
    explicit ModelSolver01_06(ModelType type);
    virtual ~ModelSolver01_06();

    // 设置计算精度（高精度模式下Stehfest项数N取值更大）
    void setHighPrecision(bool high);

    // 核心计算接口：根据输入的物理参数和时间序列计算理论压力和导数曲线
    // 参数 params: 包含 kf, M12, phi, mu, L, Lf, omega1/2, lambda1/2, eta 等物理参数的字典
    // 参数 providedTime: 如果为空，则自动生成对数时间步长
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    // 静态辅助函数：获取模型对应的中文名称，用于UI显示
    static QString getModelName(ModelType type);

    // 静态辅助函数：生成对数分布的时间步长序列（用于生成平滑的对数坐标曲线）
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    // 内部函数：通过Stehfest数值反演算法，计算无因次压力(PD)和无因次导数(Deriv)
    // 根据 MATLAB 逻辑，默认 N=10 (MATLAB文件示例中为4，但为了稳定性建议保持较高精度，参数可控)
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    // 内部函数：拉普拉斯空间下的复合模型总函数 (包含双重介质、井储和表皮效应)
    // 修正：此处逻辑已更新为匹配 Composite_shale_oil_reservoir_fitfun 中的 fs1/fs2 算法
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 内部函数：计算点源解的拉普拉斯变换值 (求解裂缝流量分布矩阵)
    // 实现了 PWD_inf / PWD_composite 的核心积分方程求解
    double PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD, int nf, const QVector<double>& xwD, ModelType type);

    // 数学辅助函数：计算缩放的第一类修正贝塞尔函数 I_v(x) * exp(-|x|)
    double scaled_besseli(int v, double x);

    // 数学辅助函数：15点高斯积分公式
    double gauss15(std::function<double(double)> f, double a, double b);

    // 数学辅助函数：自适应高斯积分，用于处理裂缝沿线的积分计算
    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);

    // 数学辅助函数：计算Stehfest反演系数
    double stefestCoefficient(int i, int N);

    // 数学辅助函数：计算阶乘
    double factorial(int n);

private:
    ModelType m_type;       // 当前选择的模型类型
    bool m_highPrecision;   // 高精度计算标志
};

#endif // MODELSOLVER01_06_H
