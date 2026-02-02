/*
 * 文件名: fittingparameterchart.h
 * 文件作用: 拟合参数图表管理类头文件
 * 修改记录:
 * 1. [修改] resetParams 增加 preserveStates 参数，用于控制是否保留拟合/显示状态。
 */

#ifndef FITTINGPARAMETERCHART_H
#define FITTINGPARAMETERCHART_H

#include <QObject>
#include <QTableWidget>
#include <QList>
#include <QMap>
#include <QTimer>
#include "modelmanager.h"

// 拟合参数结构体
struct FitParameter {
    QString name;       // 英文名 (如 kf)
    QString displayName;// 中文显示名 (如 渗透率)
    double value;       // 当前值
    double min;         // 拟合下限
    double max;         // 拟合上限
    double step;        // 滚轮步长
    bool isFit;         // 是否参与拟合
    bool isVisible;     // 是否在表格中显示
};

class FittingParameterChart : public QObject
{
    Q_OBJECT
public:
    explicit FittingParameterChart(QTableWidget* parentTable, QObject *parent = nullptr);

    void setModelManager(ModelManager* m);

    // 重置为指定模型的默认参数，并自动计算上下限
    // [修改] preserveStates: 若为 true，则仅重置数值，保留 isFit 和 isVisible 状态
    void resetParams(ModelManager::ModelType type, bool preserveStates = false);

    // 切换模型（尽可能保留同名参数值）
    void switchModel(ModelManager::ModelType newType);

    // 获取/设置参数
    QList<FitParameter> getParameters() const;
    void setParameters(const QList<FitParameter>& params);

    // 更新内部参数列表（从界面表格读取）
    void updateParamsFromTable();

    // 根据当前值自动调整所有参数的上下限和步长
    void autoAdjustLimits();

    // 辅助：获取静态翻译信息
    static void getParamDisplayInfo(const QString& name, QString& chName, QString& symbol, QString& uniSym, QString& unit);
    // 辅助：获取该模型默认显示的参数列表
    static QStringList getDefaultFitKeys(ModelManager::ModelType type);

    // 获取表格中的原始文本（用于支持敏感性分析的逗号分隔输入）
    QMap<QString, QString> getRawParamTexts() const;

signals:
    void parameterChangedByWheel(); // 滚轮改变参数信号

protected:
    bool eventFilter(QObject *watched, QEvent *event) override;

private:
    QTableWidget* m_table;
    ModelManager* m_modelManager;
    QList<FitParameter> m_params;
    QTimer* m_wheelTimer; // 滚轮防抖定时器

    void refreshParamTable();
    void addRowToTable(const FitParameter& p, int& serialNo, bool highlight);

private slots:
    void onTableItemChanged(QTableWidgetItem* item);
    void onWheelDebounceTimeout();
};

#endif // FITTINGPARAMETERCHART_H
