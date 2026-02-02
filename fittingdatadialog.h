/*
 * 文件名: fittingdatadialog.h
 * 文件作用: 拟合数据加载配置窗口头文件
 * 功能描述:
 * 1. 定义数据加载设置结构体 FittingDataSettings。
 * 2. 声明数据加载对话框类，支持文件选择、列映射及试井参数设置。
 */

#ifndef FITTINGDATADIALOG_H
#define FITTINGDATADIALOG_H

#include <QDialog>
#include <QMap>
#include <QStandardItemModel>

namespace Ui {
class FittingDataDialog;
}

// 试井类型枚举
enum WellTestType {
    Test_Drawdown, // 压力降落
    Test_Buildup   // 压力恢复
};

// 数据设置结构体
struct FittingDataSettings {
    bool isFromProject;         // 是否来自项目内部数据
    QString projectFileName;    // 项目文件名 (如果是项目数据)
    QString filePath;           // 外部文件路径 (如果是外部数据)

    int timeColIndex;           // 时间列索引
    int pressureColIndex;       // 压力列索引
    int derivColIndex;          // 导数列索引 (-1表示自动计算)
    int skipRows;               // 跳过表头行数

    WellTestType testType;      // 试井类型
    double initialPressure;     // 初始地层压力 (仅降落试井有效)
    double producingTime;       // [新增] 关井前生产时间 (仅恢复试井有效)

    double lSpacing;            // 导数计算步长
    bool enableSmoothing;       // 是否启用平滑
    int smoothingSpan;          // 平滑窗口
};

class FittingDataDialog : public QDialog
{
    Q_OBJECT

public:
    explicit FittingDataDialog(const QMap<QString, QStandardItemModel*>& projectModels, QWidget *parent = nullptr);
    ~FittingDataDialog();

    // 获取设置结果
    FittingDataSettings getSettings() const;

    // 获取预览用的数据模型
    QStandardItemModel* getPreviewModel() const;

private slots:
    void onSourceChanged();
    void onProjectFileSelectionChanged(int index);
    void onBrowseFile();
    void onAccepted();
    void onDerivColumnChanged(int index);
    void onTestTypeChanged();
    void onSmoothingToggled(bool checked);

private:
    Ui::FittingDataDialog *ui;
    QMap<QString, QStandardItemModel*> m_projectDataMap;
    QStandardItemModel* m_fileModel;

    // 解析辅助函数
    bool parseTextFile(const QString& filePath);
    bool parseExcelFile(const QString& filePath);
    QStandardItemModel* getCurrentProjectModel() const;
    void updateColumnComboBoxes(const QStringList& headers);
};

#endif // FITTINGDATADIALOG_H
