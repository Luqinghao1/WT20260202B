/*
 * modelsolver01-06.h
 * 文件作用: 压裂水平井复合页岩油模型核心计算类头文件
 * 修改记录:
 * 1. 更新了关于裂缝参数的注释说明 (nf 为裂缝条数, n_seg 为单缝段数)。
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
    // 模型类型枚举
    enum ModelType {
        Model_1 = 0, // 模型1: 变井储 + 无限大边界
        Model_2,     // 模型2: 恒定井储 + 无限大边界
        Model_3,     // 模型3: 变井储 + 封闭边界
        Model_4,     // 模型4: 恒定井储 + 封闭边界
        Model_5,     // 模型5: 变井储 + 定压边界
        Model_6      // 模型6: 恒定井储 + 定压边界
    };

    explicit ModelSolver01_06(ModelType type);
    virtual ~ModelSolver01_06();

    void setHighPrecision(bool high);

    // 核心计算接口
    // params 关键参数变更:
    // "nf"      -> 裂缝条数 (Number of Fractures)
    // "n_seg"   -> 单条裂缝的离散微元段数 (原 nf)
    // "frac_spacing" -> 裂缝间距 (m)
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    static QString getModelName(ModelType type);
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 内部函数：PWD_composite 签名更新，增加了 spacingD 和 n_fracs 参数
    // n_seg: 单条裂缝的段数
    // n_fracs: 裂缝条数
    // spacingD: 无因次裂缝间距
    double PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                         int n_seg, int n_fracs, double spacingD, ModelType type);

    double scaled_besseli(int v, double x);
    double gauss15(std::function<double(double)> f, double a, double b);
    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);
    double stefestCoefficient(int i, int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
};

#endif // MODELSOLVER01_06_H
