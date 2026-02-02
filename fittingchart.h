/*
 * 文件名: fittingchart.h
 * 文件作用: 拟合绘图管理类头文件
 * 修改记录:
 * 1. [修复] 文本框拖拽逻辑重写，解决乱跑问题。
 * 2. [优化] 接管右键菜单，实现原有功能+压力求解功能的融合。
 * 3. [优化] 拟合线支持样式编辑，但禁止移动/旋转。
 * 4. [新增] 完善交互逻辑，支持拟合线伸缩、标注添加、删除功能。
 * 5. [修复] 将 ChartAnnotation 重命名为 FittingChartAnnotation 以避免与 ChartWidget 冲突。
 */

#ifndef FITTINGCHART_H
#define FITTINGCHART_H

#include <QObject>
#include <QVector>
#include <QPointer>
#include <QAction>
#include <QPointF>
#include <QMap>
#include "mousezoom.h"
#include "fittingdatadialog.h"
#include "fittingpressuredialog.h"
#include "styleselectordialog.h"

// 标注结构体 (复用 ChartWidget 的设计理念)
// [修改] 重命名为 FittingChartAnnotation 避免与 chartwidget.h 中的定义冲突
struct FittingChartAnnotation {
    QCPItemText* textItem = nullptr;
    QCPItemLine* arrowItem = nullptr;
};

class FittingChart : public QObject
{
    Q_OBJECT
public:
    explicit FittingChart(QObject *parent = nullptr);

    // 初始化图表指针
    void initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian);

    void setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                         const QVector<double>& deriv, const QVector<double>& rawP);

    void setSettings(const FittingDataSettings& settings);

    // 绘制所有曲线
    void plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid, bool autoScale = false);

    // 绘制抽样点
    void plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

    double getCalculatedInitialPressure() const;

private slots:
    // 半对数图右键菜单
    void onSemiLogContextMenu(const QPoint &pos);

    // 交互槽函数
    void onShowPressureSolver();

    // 鼠标事件处理
    void onPlotMousePress(QMouseEvent *event);
    void onPlotMouseMove(QMouseEvent *event);
    void onPlotMouseRelease(QMouseEvent *event);

    // 弹窗交互
    void onPickStart();
    void onPickEnd();
    void onCalculatePressure();

    // 功能槽函数
    void onLineStyleRequested(QCPItemLine* line);
    void onAddAnnotationRequested(QCPItemLine* line);
    void onDeleteSelectedRequested();
    void onEditItemRequested(QCPAbstractItem* item);

private:
    MouseZoom* m_plotLogLog;
    MouseZoom* m_plotSemiLog;
    MouseZoom* m_plotCartesian;

    QVector<double> m_obsT;
    QVector<double> m_obsDeltaP;
    QVector<double> m_obsDeriv;
    QVector<double> m_obsRawP;

    FittingDataSettings m_settings;
    double m_calculatedPi;

    // --- 原始地层压力求解相关 ---
    QPointer<FittingPressureDialog> m_pressureDialog;
    bool m_isPickingStart;
    bool m_isPickingEnd;

    // 结果缓存
    bool m_hasManualPressure;
    double m_manualSlope;
    double m_manualIntercept;
    double m_manualStartX;
    double m_manualEndX;

    // 绘图项指针
    QPointer<QCPItemLine> m_manualFitLine;    // 拟合红线 (P*)
    QPointer<QCPItemLine> m_manualZeroLine;   // X=0 蓝线
    QPointer<QCPItemText> m_manualResultText; // 结果文本框

    // --- 交互状态管理 ---
    enum InteractionMode {
        Mode_None,
        Mode_Dragging_Start,       // 拖拽线段起点 (伸缩)
        Mode_Dragging_End,         // 拖拽线段终点 (伸缩)
        Mode_Dragging_Text,        // 拖拽文本
        Mode_Dragging_ArrowStart,  // 拖拽箭头起点 (标注)
        Mode_Dragging_ArrowEnd     // 拖拽箭头终点 (标注)
    };

    InteractionMode m_interMode;
    QPointF m_lastMousePos;

    // 当前选中的活动项
    QCPItemLine* m_activeLine;
    QCPItemText* m_activeText;
    QCPItemLine* m_activeArrow;

    // 标注管理
    QMap<QCPItemLine*, FittingChartAnnotation> m_annotations;

    // 内部绘图函数
    void plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);
    void plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);
    void plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);

    void showResultOnLogPlot();
    void drawPressureFitResult();

    // 辅助函数
    double distToSegment(const QPointF& p, const QPointF& s, const QPointF& e);
    // 限制线段端点移动 (只允许沿线伸缩)
    void constrainLinePoint(QCPItemLine* line, bool isMovingStart, double mouseX, double mouseY);
    // 更新标注箭头位置
    void updateAnnotationArrow(QCPItemLine* line);
    // 添加标注到线
    void addAnnotationToLine(QCPItemLine* line);
};

#endif // FITTINGCHART_H
