/*
 * 文件名: fittingchart.cpp
 * 文件作用: 拟合绘图管理类实现文件
 * 修改记录:
 * 1. [修复] 文本框拖拽逻辑重写，解决乱跑问题。
 * 2. [交互] 实现了拟合线的伸缩交互，限制了平移和旋转。
 * 3. [功能] 增加了标注添加、删除、样式修改功能，并在删除拟合线时同步删除结果文本。
 * 4. [修复] 使用 FittingChartAnnotation 替代 ChartAnnotation。
 * 5. [优化] 结果文本框改为 ptPlotCoords 模式，拖动交互与标注一致。
 */

#include "fittingchart.h"
#include <cmath>
#include <algorithm>
#include <QDebug>
#include <QMenu>
#include <QAction>
#include <QMouseEvent>
#include <QMetaObject>
#include <QInputDialog>

FittingChart::FittingChart(QObject *parent)
    : QObject(parent)
    , m_plotLogLog(nullptr), m_plotSemiLog(nullptr), m_plotCartesian(nullptr)
    , m_calculatedPi(0.0)
    , m_isPickingStart(false)
    , m_isPickingEnd(false)
    , m_hasManualPressure(false)
    , m_manualSlope(0)
    , m_manualIntercept(0)
    , m_manualStartX(0)
    , m_manualEndX(0)
    , m_interMode(Mode_None)
    , m_activeLine(nullptr)
    , m_activeText(nullptr)
    , m_activeArrow(nullptr)
{
}

void FittingChart::initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian)
{
    m_plotLogLog = logLog;
    m_plotSemiLog = semiLog;
    m_plotCartesian = cartesian;

    // --- 1. 初始化双对数图 (Log-Log) ---
    if (m_plotLogLog) {
        m_plotLogLog->xAxis->setLabel("时间 Time (h)");
        m_plotLogLog->yAxis->setLabel("压差 & 导数 (MPa)");

        QSharedPointer<QCPAxisTickerLog> logTickerX(new QCPAxisTickerLog);
        logTickerX->setLogBase(10.0);
        m_plotLogLog->xAxis->setTicker(logTickerX);
        m_plotLogLog->xAxis->setScaleType(QCPAxis::stLogarithmic);
        m_plotLogLog->xAxis->setNumberFormat("gb");

        QSharedPointer<QCPAxisTickerLog> logTickerY(new QCPAxisTickerLog);
        logTickerY->setLogBase(10.0);
        m_plotLogLog->yAxis->setTicker(logTickerY);
        m_plotLogLog->yAxis->setScaleType(QCPAxis::stLogarithmic);
        m_plotLogLog->yAxis->setNumberFormat("gb");

        if(m_plotLogLog->axisRect() && m_plotLogLog->axisRect()->insetLayout())
            m_plotLogLog->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
    }

    // --- 2. 初始化笛卡尔图 (Cartesian) ---
    if (m_plotCartesian) {
        m_plotCartesian->xAxis->setLabel("时间 Time (h)");
        m_plotCartesian->yAxis->setLabel("压差 Delta P (MPa)");

        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        m_plotCartesian->xAxis->setTicker(linearTicker);
        m_plotCartesian->xAxis->setScaleType(QCPAxis::stLinear);
        m_plotCartesian->xAxis->setNumberFormat("gb");

        m_plotCartesian->yAxis->setTicker(linearTicker);
        m_plotCartesian->yAxis->setScaleType(QCPAxis::stLinear);
        m_plotCartesian->yAxis->setNumberFormat("gb");

        if(m_plotCartesian->axisRect() && m_plotCartesian->axisRect()->insetLayout())
            m_plotCartesian->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom | Qt::AlignRight);
    }

    // --- 3. 半对数图 (Semi-Log / Horner) ---
    if (m_plotSemiLog) {
        if(m_plotSemiLog->axisRect() && m_plotSemiLog->axisRect()->insetLayout())
            m_plotSemiLog->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignRight);

        // 断开 MouseZoom 自带的右键菜单连接
        m_plotSemiLog->disconnect(SIGNAL(customContextMenuRequested(QPoint)));

        // 开启自定义右键菜单策略
        m_plotSemiLog->setContextMenuPolicy(Qt::CustomContextMenu);
        connect(m_plotSemiLog, &QWidget::customContextMenuRequested,
                this, &FittingChart::onSemiLogContextMenu);

        // 连接鼠标事件（用于交互）
        connect(m_plotSemiLog, &QCustomPlot::mousePress,
                this, &FittingChart::onPlotMousePress);
        connect(m_plotSemiLog, &QCustomPlot::mouseMove,
                this, &FittingChart::onPlotMouseMove);
        connect(m_plotSemiLog, &QCustomPlot::mouseRelease,
                this, &FittingChart::onPlotMouseRelease);

        // 增加双击编辑文本功能
        connect(m_plotSemiLog, &QCustomPlot::mouseDoubleClick, [=](QMouseEvent* event){
            if (event->button() != Qt::LeftButton) return;
            for (int i = 0; i < m_plotSemiLog->itemCount(); ++i) {
                if (auto text = qobject_cast<QCPItemText*>(m_plotSemiLog->item(i))) {
                    if (text->selectTest(event->pos(), false) < 10.0) { onEditItemRequested(text); return; }
                }
            }
        });
    }
}

// 辅助函数：计算点到线段的距离
double FittingChart::distToSegment(const QPointF& p, const QPointF& s, const QPointF& e)
{
    double l2 = (s.x()-e.x())*(s.x()-e.x()) + (s.y()-e.y())*(s.y()-e.y());
    if (l2 == 0) return std::sqrt((p.x()-s.x())*(p.x()-s.x()) + (p.y()-s.y())*(p.y()-s.y()));
    double t = ((p.x()-s.x())*(e.x()-s.x()) + (p.y()-s.y())*(e.y()-s.y())) / l2;
    t = std::max(0.0, std::min(1.0, t));
    QPointF proj = s + t * (e - s);
    return std::sqrt((p.x()-proj.x())*(p.x()-proj.x()) + (p.y()-proj.y())*(p.y()-proj.y()));
}

// 右键菜单处理
void FittingChart::onSemiLogContextMenu(const QPoint &pos)
{
    if (!m_plotSemiLog) return;

    QMenu *menu = new QMenu(m_plotSemiLog);
    QPointF pMouse = pos;
    double tolerance = 8.0;

    // 1. 查找是否击中线条 (拟合线 或 零线)
    QCPItemLine* hitLine = nullptr;
    QList<QCPItemLine*> candidates;
    if (m_manualFitLine) candidates << m_manualFitLine;
    if (m_manualZeroLine) candidates << m_manualZeroLine;

    // 同时也检查是否有普通标识线被击中 (如果有的话)
    for(int i=0; i<m_plotSemiLog->itemCount(); ++i) {
        if(auto l = qobject_cast<QCPItemLine*>(m_plotSemiLog->item(i))) {
            if(!candidates.contains(l) && !l->property("isArrow").toBool()) candidates << l;
        }
    }

    for (auto line : candidates) {
        if (!line->visible()) continue;
        double x1 = m_plotSemiLog->xAxis->coordToPixel(line->start->coords().x());
        double y1 = m_plotSemiLog->yAxis->coordToPixel(line->start->coords().y());
        double x2 = m_plotSemiLog->xAxis->coordToPixel(line->end->coords().x());
        double y2 = m_plotSemiLog->yAxis->coordToPixel(line->end->coords().y());

        // 检查端点或线身
        double dStart = std::sqrt(pow(pMouse.x()-x1,2) + pow(pMouse.y()-y1,2));
        double dEnd   = std::sqrt(pow(pMouse.x()-x2,2) + pow(pMouse.y()-y2,2));
        double dLine  = distToSegment(pMouse, QPointF(x1, y1), QPointF(x2, y2));

        if (dStart < tolerance || dEnd < tolerance || dLine < tolerance) {
            hitLine = line;
            break;
        }
    }

    // 2. 查找是否击中文本
    QCPItemText* hitText = nullptr;
    for (int i = 0; i < m_plotSemiLog->itemCount(); ++i) {
        if (auto text = qobject_cast<QCPItemText*>(m_plotSemiLog->item(i))) {
            if (text->selectTest(pMouse, false) < tolerance) {
                hitText = text;
                break;
            }
        }
    }

    // 3. 构建菜单内容
    if (hitLine) {
        m_activeLine = hitLine; // 标记当前操作线

        QAction *actStyle = menu->addAction("线条设置 (颜色/线型)");
        connect(actStyle, &QAction::triggered, [=](){ this->onLineStyleRequested(hitLine); });

        QAction *actAnnotate = menu->addAction("添加/修改标注");
        connect(actAnnotate, &QAction::triggered, [=](){ this->onAddAnnotationRequested(hitLine); });

        menu->addSeparator();
        QAction *actDelete = menu->addAction("删除线段");
        connect(actDelete, &QAction::triggered, this, &FittingChart::onDeleteSelectedRequested);
    }
    else if (hitText) {
        QAction *actEdit = menu->addAction("编辑文本");
        connect(actEdit, &QAction::triggered, [=](){ this->onEditItemRequested(hitText); });

        if (hitText == m_manualResultText) {
            QAction *actResetPos = menu->addAction("重置位置");
            connect(actResetPos, &QAction::triggered, [=](){
                if(m_manualResultText && m_plotSemiLog) {
                    // [修改] 重置到当前视图的中心位置 (因为现在是 ptPlotCoords)
                    double pixelX = m_plotSemiLog->axisRect()->left() + m_plotSemiLog->axisRect()->width() * 0.05;
                    double pixelY = m_plotSemiLog->axisRect()->top() + m_plotSemiLog->axisRect()->height() * 0.35;
                    m_manualResultText->position->setCoords(
                        m_plotSemiLog->xAxis->pixelToCoord(pixelX),
                        m_plotSemiLog->yAxis->pixelToCoord(pixelY)
                        );
                    m_plotSemiLog->replot();
                }
            });
        } else {
            QAction *actDelText = menu->addAction("删除文本");
            connect(actDelText, &QAction::triggered, [=](){
                m_plotSemiLog->removeItem(hitText);
                m_plotSemiLog->replot();
            });
        }
    }

    // 4. 添加通用菜单项
    if (!hitLine && !hitText) {
        QAction *actSave = menu->addAction("导出图片");
        connect(actSave, &QAction::triggered, [=](){ QMetaObject::invokeMethod(m_plotSemiLog, "saveImageRequested"); });

        QAction *actExport = menu->addAction("导出数据");
        connect(actExport, &QAction::triggered, [=](){ QMetaObject::invokeMethod(m_plotSemiLog, "exportDataRequested"); });

        QMenu* subMenuLine = menu->addMenu("标识线绘制");
        auto addDrawAction = [&](QString name, double k){
            QAction* act = subMenuLine->addAction(name);
            connect(act, &QAction::triggered, [=](){ QMetaObject::invokeMethod(m_plotSemiLog, "drawLineRequested", Q_ARG(double, k)); });
        };
        addDrawAction("斜率 k=1", 1.0);
        addDrawAction("斜率 k=1/2", 0.5);
        addDrawAction("斜率 k=1/4", 0.25);
        addDrawAction("水平线", 0.0);

        QAction *actSetting = menu->addAction("图表设置");
        connect(actSetting, &QAction::triggered, [=](){ QMetaObject::invokeMethod(m_plotSemiLog, "settingsRequested"); });

        menu->addSeparator();

        QAction *actReset = menu->addAction("重置视图");
        connect(actReset, &QAction::triggered, [=](){ QMetaObject::invokeMethod(m_plotSemiLog, "resetViewRequested"); });

        menu->addSeparator();
        QAction *actionSolve = menu->addAction("原始地层压力");
        connect(actionSolve, &QAction::triggered, this, &FittingChart::onShowPressureSolver);
    }

    menu->exec(m_plotSemiLog->mapToGlobal(pos));
    delete menu;
}

void FittingChart::onLineStyleRequested(QCPItemLine* line)
{
    if (!line) return;
    StyleSelectorDialog dlg(StyleSelectorDialog::ModeLine, m_plotSemiLog);
    dlg.setPen(line->pen());
    dlg.setWindowTitle("样式设置");
    if (dlg.exec() == QDialog::Accepted) {
        line->setPen(dlg.getPen());
        m_plotSemiLog->replot();
    }
}

void FittingChart::onAddAnnotationRequested(QCPItemLine* line)
{
    if (!line) return;

    // 如果已有标注，先清除旧的
    if (m_annotations.contains(line)) {
        FittingChartAnnotation old = m_annotations.take(line);
        if(old.textItem) m_plotSemiLog->removeItem(old.textItem);
        if(old.arrowItem) m_plotSemiLog->removeItem(old.arrowItem);
    }

    bool ok;
    // 默认显示斜率信息
    QString defaultText = "Annotation";
    if (line == m_manualFitLine) {
        defaultText = QString("k=%1").arg(m_manualSlope, 0, 'f', 4);
    }

    QString text = QInputDialog::getText(m_plotSemiLog, "添加标注", "输入标注内容:", QLineEdit::Normal, defaultText, &ok);
    if (!ok || text.isEmpty()) return;

    QCPItemText* txt = new QCPItemText(m_plotSemiLog);
    txt->setText(text);
    txt->position->setType(QCPItemPosition::ptPlotCoords);

    // 计算位置：线段中点上方
    double midX = (line->start->coords().x() + line->end->coords().x()) / 2.0;
    double midY = (line->start->coords().y() + line->end->coords().y()) / 2.0;

    // 简单偏移
    double yOffset = (m_plotSemiLog->yAxis->range().upper - m_plotSemiLog->yAxis->range().lower) * 0.05;
    txt->position->setCoords(midX, midY + yOffset);
    txt->setFont(QFont("Microsoft YaHei", 9));

    QCPItemLine* arr = new QCPItemLine(m_plotSemiLog);
    arr->setHead(QCPLineEnding::esSpikeArrow);
    arr->start->setParentAnchor(txt->bottom); // 箭头起点跟随文本底部
    arr->end->setCoords(midX, midY);          // 箭头终点指向线中点
    arr->setProperty("isArrow", true); // 标记为箭头，避免被误选为操作线

    FittingChartAnnotation note;
    note.textItem = txt;
    note.arrowItem = arr;
    m_annotations.insert(line, note);

    m_plotSemiLog->replot();
}

void FittingChart::onDeleteSelectedRequested()
{
    if (!m_activeLine) return;

    // 如果删除的是拟合线，需要清除结果文本和状态
    if (m_activeLine == m_manualFitLine) {
        if (m_manualResultText) {
            m_plotSemiLog->removeItem(m_manualResultText);
            m_manualResultText = nullptr;
        }
        m_hasManualPressure = false;
        m_manualFitLine = nullptr;
    }
    else if (m_activeLine == m_manualZeroLine) {
        m_manualZeroLine = nullptr;
    }

    // 清除关联的标注
    if (m_annotations.contains(m_activeLine)) {
        FittingChartAnnotation note = m_annotations.take(m_activeLine);
        if(note.textItem) m_plotSemiLog->removeItem(note.textItem);
        if(note.arrowItem) m_plotSemiLog->removeItem(note.arrowItem);
    }

    m_plotSemiLog->removeItem(m_activeLine);
    m_activeLine = nullptr;
    m_plotSemiLog->replot();
}

void FittingChart::onEditItemRequested(QCPAbstractItem* item) {
    if (auto text = qobject_cast<QCPItemText*>(item)) {
        bool ok;
        QString newContent = QInputDialog::getText(m_plotSemiLog, "修改标注", "内容:", QLineEdit::Normal, text->text(), &ok);
        if (ok && !newContent.isEmpty()) { text->setText(newContent); m_plotSemiLog->replot(); }
    }
}

void FittingChart::onShowPressureSolver()
{
    if (!m_pressureDialog) {
        m_pressureDialog = new FittingPressureDialog(nullptr);
        connect(m_pressureDialog, &FittingPressureDialog::requestPickStart, this, &FittingChart::onPickStart);
        connect(m_pressureDialog, &FittingPressureDialog::requestPickEnd, this, &FittingChart::onPickEnd);
        connect(m_pressureDialog, &FittingPressureDialog::requestCalculate, this, &FittingChart::onCalculatePressure);
    }
    m_pressureDialog->show();
    m_pressureDialog->raise();
    m_pressureDialog->activateWindow();
}

void FittingChart::onPickStart()
{
    m_isPickingStart = true;
    m_isPickingEnd = false;
    if (m_plotSemiLog) m_plotSemiLog->setCursor(Qt::CrossCursor);
}

void FittingChart::onPickEnd()
{
    m_isPickingStart = false;
    m_isPickingEnd = true;
    if (m_plotSemiLog) m_plotSemiLog->setCursor(Qt::CrossCursor);
}

// ============================================================================
// 鼠标交互核心逻辑
// ============================================================================

void FittingChart::onPlotMousePress(QMouseEvent *event)
{
    if (!m_plotSemiLog) return;

    // 1. 优先处理取点逻辑
    if (m_isPickingStart || m_isPickingEnd) {
        double x = m_plotSemiLog->xAxis->pixelToCoord(event->pos().x());
        double y = m_plotSemiLog->yAxis->pixelToCoord(event->pos().y());

        if (m_pressureDialog) {
            if (m_isPickingStart) {
                m_pressureDialog->setStartCoordinate(x, y);
                m_isPickingStart = false;
            } else if (m_isPickingEnd) {
                m_pressureDialog->setEndCoordinate(x, y);
                m_isPickingEnd = false;
            }
        }
        m_plotSemiLog->setCursor(Qt::ArrowCursor);
        return;
    }

    // 重置状态
    m_interMode = Mode_None;
    m_activeLine = nullptr;
    m_activeText = nullptr;
    m_activeArrow = nullptr;
    m_lastMousePos = event->pos();
    double tolerance = 8.0;

    // 2. 检测文本 (优先选中文字以便拖拽)
    for (int i = 0; i < m_plotSemiLog->itemCount(); ++i) {
        if (auto text = qobject_cast<QCPItemText*>(m_plotSemiLog->item(i))) {
            if (text->selectTest(event->pos(), false) < tolerance) {
                m_interMode = Mode_Dragging_Text;
                m_activeText = text;
                // 选中高亮
                m_plotSemiLog->deselectAll();
                text->setSelected(true);
                m_plotSemiLog->setInteractions(QCP::Interactions()); // 禁用默认交互
                m_plotSemiLog->replot();
                return;
            }
        }
    }

    // 3. 检测线条端点 (用于伸缩)
    QList<QCPItemLine*> draggableLines;
    if (m_manualFitLine) draggableLines << m_manualFitLine;
    if (m_manualZeroLine) draggableLines << m_manualZeroLine;
    // 以及所有有 FittingChartAnnotation 的线条 (视为用户添加的标识线)
    // 这里为了统一，也可以遍历所有 line

    for (auto line : draggableLines) {
        if (!line->visible()) continue;
        double x1 = m_plotSemiLog->xAxis->coordToPixel(line->start->coords().x());
        double y1 = m_plotSemiLog->yAxis->coordToPixel(line->start->coords().y());
        double x2 = m_plotSemiLog->xAxis->coordToPixel(line->end->coords().x());
        double y2 = m_plotSemiLog->yAxis->coordToPixel(line->end->coords().y());

        QPointF p(event->pos());
        double dStart = std::sqrt(pow(p.x()-x1,2) + pow(p.y()-y1,2));
        double dEnd   = std::sqrt(pow(p.x()-x2,2) + pow(p.y()-y2,2));
        double dLine  = distToSegment(p, QPointF(x1, y1), QPointF(x2, y2));

        // 选中判断
        if (dStart < tolerance) {
            m_interMode = Mode_Dragging_Start;
            m_activeLine = line;
        } else if (dEnd < tolerance) {
            m_interMode = Mode_Dragging_End;
            m_activeLine = line;
        } else if (dLine < tolerance) {
            // 仅选中，不拖拽 (禁止平移)
            m_activeLine = line;
        }

        if (m_activeLine) {
            m_plotSemiLog->deselectAll();
            m_activeLine->setSelected(true);
            m_plotSemiLog->setInteractions(QCP::Interactions()); // 锁定默认交互
            m_plotSemiLog->replot();
            return;
        }
    }

    // 4. 点击空白处
    m_plotSemiLog->deselectAll();
    m_plotSemiLog->replot();
}

void FittingChart::onPlotMouseMove(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton) {
        if (!m_plotSemiLog) return;

        QPointF currentPos = event->pos();
        QPointF delta = currentPos - m_lastMousePos;
        double mouseX = m_plotSemiLog->xAxis->pixelToCoord(currentPos.x());
        double mouseY = m_plotSemiLog->yAxis->pixelToCoord(currentPos.y());

        // --- 文本拖拽处理 ---
        if (m_interMode == Mode_Dragging_Text && m_activeText) {
            // 区分坐标类型
            if (m_activeText->position->type() == QCPItemPosition::ptPlotCoords) {
                // 坐标系模式：通过像素偏差计算新坐标，确保 1:1 拖拽手感
                double px = m_plotSemiLog->xAxis->coordToPixel(m_activeText->position->coords().x()) + delta.x();
                double py = m_plotSemiLog->yAxis->coordToPixel(m_activeText->position->coords().y()) + delta.y();
                m_activeText->position->setCoords(m_plotSemiLog->xAxis->pixelToCoord(px), m_plotSemiLog->yAxis->pixelToCoord(py));
            }
            else if (m_activeText->position->type() == QCPItemPosition::ptAxisRectRatio) {
                // 比例模式：将 pixel delta 转换为 ratio delta
                double w = m_plotSemiLog->axisRect()->width();
                double h = m_plotSemiLog->axisRect()->height();
                if (w > 0 && h > 0) {
                    double currentRX = m_activeText->position->coords().x();
                    double currentRY = m_activeText->position->coords().y();
                    m_activeText->position->setCoords(currentRX + delta.x()/w, currentRY + delta.y()/h);
                }
            }
            // 联动箭头起点
            // 注意：FittingChartAnnotation 中箭头起点是 setParentAnchor 到文本的，会自动跟随，无需手动更新
        }
        // --- 线条伸缩处理 ---
        else if ((m_interMode == Mode_Dragging_Start || m_interMode == Mode_Dragging_End) && m_activeLine) {
            constrainLinePoint(m_activeLine, m_interMode == Mode_Dragging_Start, mouseX, mouseY);
            updateAnnotationArrow(m_activeLine);
        }

        m_plotSemiLog->replot();
        m_lastMousePos = currentPos;
    }
}

void FittingChart::constrainLinePoint(QCPItemLine* line, bool isMovingStart, double mouseX, double mouseY)
{
    // 1. 如果是 Zero Line (x=0)
    if (line == m_manualZeroLine) {
        // 强制 X=0，只允许 Y 变化
        if (isMovingStart) line->start->setCoords(0, mouseY);
        else line->end->setCoords(0, mouseY);
        return;
    }

    // 2. 如果是 Fit Line (P* 线)，必须沿 y = kx + b 移动
    if (line == m_manualFitLine) {
        // 由于是半对数图(Horner)，X 轴实际上是 log 值，Y 轴是 P
        // 拟合方程 y = kx + b 是在 PlotCoords 空间成立的
        // 为方便交互，我们使用 X 轴主导：根据鼠标 X 计算对应的 Y
        double newY = m_manualSlope * mouseX + m_manualIntercept;

        if (isMovingStart) line->start->setCoords(mouseX, newY);
        else line->end->setCoords(mouseX, newY);
        return;
    }

    // 3. 其他普通线 (如果有)，允许自由移动端点
    if (isMovingStart) line->start->setCoords(mouseX, mouseY);
    else line->end->setCoords(mouseX, mouseY);
}

void FittingChart::updateAnnotationArrow(QCPItemLine* line) {
    if (m_annotations.contains(line)) {
        FittingChartAnnotation note = m_annotations[line];
        double midX = (line->start->coords().x() + line->end->coords().x()) / 2.0;
        double midY = (line->start->coords().y() + line->end->coords().y()) / 2.0;
        if(note.arrowItem) note.arrowItem->end->setCoords(midX, midY);
    }
}

void FittingChart::onPlotMouseRelease(QMouseEvent *event)
{
    Q_UNUSED(event);
    if (m_interMode != Mode_None) {
        // 恢复默认交互
        if (m_plotSemiLog) {
            m_plotSemiLog->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectItems);
        }
    }
    m_interMode = Mode_None;
}

void FittingChart::onCalculatePressure()
{
    if (!m_pressureDialog) return;

    double x1 = m_pressureDialog->getStartX();
    double x2 = m_pressureDialog->getEndX();
    if (x1 > x2) std::swap(x1, x2);

    QVector<double> X, Y;
    double tp = m_settings.producingTime;
    bool useHorner = (tp > 1e-5);

    for(int i=0; i<m_obsT.size(); ++i) {
        double dt = m_obsT[i];
        if (dt < 1e-6) continue;

        double x_val;
        if (useHorner) {
            double ratio = (tp + dt) / dt;
            if (ratio <= 0) continue;
            x_val = log10(ratio);
        } else {
            x_val = dt;
        }

        if (x_val >= x1 && x_val <= x2) {
            if (i < m_obsRawP.size()) {
                X.append(x_val);
                Y.append(m_obsRawP[i]);
            }
        }
    }

    if (X.size() < 2) return;

    double sumX=0, sumY=0, sumXY=0, sumXX=0;
    int n = X.size();
    for(int i=0; i<n; ++i) {
        sumX += X[i];
        sumY += Y[i];
        sumXY += X[i] * Y[i];
        sumXX += X[i] * X[i];
    }

    double denominator = n * sumXX - sumX * sumX;
    if (std::abs(denominator) < 1e-9) return;

    double k = (n * sumXY - sumX * sumY) / denominator;
    double b = (sumY - k * sumX) / n;

    m_hasManualPressure = true;
    m_manualSlope = k;
    m_manualIntercept = b;
    m_manualStartX = x1;
    m_manualEndX = x2;
    m_calculatedPi = b;

    drawPressureFitResult();
    m_plotSemiLog->replot();
}

void FittingChart::drawPressureFitResult()
{
    if (!m_plotSemiLog || !m_hasManualPressure) return;

    if (m_manualFitLine) m_plotSemiLog->removeItem(m_manualFitLine);
    if (m_manualZeroLine) m_plotSemiLog->removeItem(m_manualZeroLine);
    if (m_manualResultText) m_plotSemiLog->removeItem(m_manualResultText);

    double x_start = std::max(m_manualStartX, m_manualEndX);
    if (x_start < m_plotSemiLog->xAxis->range().upper) x_start = m_plotSemiLog->xAxis->range().upper;

    // 1. 拟合线
    m_manualFitLine = new QCPItemLine(m_plotSemiLog);
    m_manualFitLine->start->setCoords(x_start, m_manualSlope * x_start + m_manualIntercept);
    m_manualFitLine->end->setCoords(0, m_manualIntercept);
    m_manualFitLine->setPen(QPen(Qt::red, 2, Qt::DashLine));
    m_manualFitLine->setSelectable(true); // 允许选中，以便触发点击事件

    // 2. 零线
    m_manualZeroLine = new QCPItemLine(m_plotSemiLog);
    m_manualZeroLine->start->setCoords(0, m_plotSemiLog->yAxis->range().lower);
    m_manualZeroLine->end->setCoords(0, m_plotSemiLog->yAxis->range().upper);
    m_manualZeroLine->setPen(QPen(Qt::blue, 1, Qt::DotLine));
    m_manualZeroLine->setSelectable(true);

    // 3. 结果文本
    m_manualResultText = new QCPItemText(m_plotSemiLog);
    // [修改] 改为 PlotCoords 模式，与 Annotation 交互一致，支持 1:1 拖拽
    m_manualResultText->position->setType(QCPItemPosition::ptPlotCoords);

    // 计算初始位置：基于轴矩形左上角相对位置转换为坐标，确保初始可见
    double initPixelX = m_plotSemiLog->axisRect()->left() + m_plotSemiLog->axisRect()->width() * 0.05;
    double initPixelY = m_plotSemiLog->axisRect()->top() + m_plotSemiLog->axisRect()->height() * 0.35;
    m_manualResultText->position->setCoords(
        m_plotSemiLog->xAxis->pixelToCoord(initPixelX),
        m_plotSemiLog->yAxis->pixelToCoord(initPixelY)
        );

    m_manualResultText->setPositionAlignment(Qt::AlignTop | Qt::AlignLeft);
    QString slopeStr = QString::number(m_manualSlope, 'f', 4);
    QString interceptStr = QString::number(m_manualIntercept, 'f', 4);
    QString sign = (m_manualIntercept >= 0) ? "+" : "";
    QString text = QString("原始地层压力 P* = %1 MPa\n拟合方程: y = %2x %3 %4")
                       .arg(interceptStr).arg(slopeStr).arg(sign).arg(interceptStr);

    m_manualResultText->setText(text);
    m_manualResultText->setFont(QFont("Microsoft YaHei", 9));
    m_manualResultText->setColor(Qt::black);
    m_manualResultText->setBrush(QBrush(QColor(255, 255, 255, 220)));
    m_manualResultText->setPen(QPen(Qt::black));
    m_manualResultText->setPadding(QMargins(4, 4, 4, 4));
    m_manualResultText->setSelectable(true);
}

void FittingChart::setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                                   const QVector<double>& deriv, const QVector<double>& rawP)
{
    m_obsT = t;
    m_obsDeltaP = deltaP;
    m_obsDeriv = deriv;
    m_obsRawP = rawP;
}

void FittingChart::setSettings(const FittingDataSettings& settings)
{
    m_settings = settings;
}

void FittingChart::plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid, bool autoScale)
{
    if (!m_plotLogLog || !m_plotSemiLog || !m_plotCartesian) return;
    plotLogLog(t_model, p_model, d_model, isModelValid, autoScale);
    plotSemiLog(t_model, p_model, d_model, isModelValid, autoScale);
    plotCartesian(t_model, p_model, d_model, isModelValid, autoScale);
}

void FittingChart::plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotLogLog;
    plot->clearGraphs();
    plot->clearItems();

    QVector<double> vt, vp, vd;
    for(int i=0; i<m_obsT.size(); ++i) {
        if (m_obsT[i] > 1e-10 && m_obsDeltaP[i] > 1e-10) {
            vt << m_obsT[i];
            vp << m_obsDeltaP[i];
            if (i < m_obsDeriv.size() && m_obsDeriv[i] > 1e-10) vd << m_obsDeriv[i];
            else vd << 1e-10;
        }
    }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    plot->addGraph();
    plot->graph(1)->setData(vt, vd);
    plot->graph(1)->setPen(Qt::NoPen);
    plot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, Qt::magenta, 6));
    plot->graph(1)->setName("实测导数");

    if (hasModel) {
        QVector<double> vtm, vpm, vdm;
        for(int i=0; i<tm.size(); ++i) {
            if(tm[i] > 1e-10) {
                vtm << tm[i];
                vpm << (pm[i] > 1e-10 ? pm[i] : 1e-10);
                vdm << (dm[i] > 1e-10 ? dm[i] : 1e-10);
            }
        }
        plot->addGraph();
        plot->graph(2)->setData(vtm, vpm);
        plot->graph(2)->setPen(QPen(Qt::red, 2));
        plot->graph(2)->setName("理论压差");

        plot->addGraph();
        plot->graph(3)->setData(vtm, vdm);
        plot->graph(3)->setPen(QPen(Qt::blue, 2));
        plot->graph(3)->setName("理论导数");
    }

    if (autoScale) {
        plot->rescaleAxes();
        plot->xAxis->scaleRange(1.1, plot->xAxis->range().center());
        plot->yAxis->scaleRange(1.1, plot->yAxis->range().center());
    }

    if (m_settings.testType == Test_Buildup && m_calculatedPi > 1e-6) {
        showResultOnLogPlot();
    }
}

void FittingChart::plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotSemiLog;
    plot->clearGraphs();

    // 注意：不要清空 items，否则拟合线会被清除。
    // 我们手动管理 items，或者只清除自动生成的。
    // 为了简单，我们只保留 manualFitLine 和相关 item，清除其他的？
    // 现有逻辑 clearItems() 会清除所有，需要重建 m_manualFitLine
    // 方案：重建
    plot->clearItems();
    m_annotations.clear(); // 清空标注记录
    m_manualFitLine = nullptr;
    m_manualZeroLine = nullptr;
    m_manualResultText = nullptr;

    Q_UNUSED(dm);

    double tp = m_settings.producingTime;
    bool useHorner = (tp > 1e-5);

    QVector<double> vt, vp;
    for(int i=0; i<m_obsT.size(); ++i) {
        double dt = m_obsT[i];
        if (dt < 1e-6) continue;

        double x_val;
        if (useHorner) {
            double ratio = (tp + dt) / dt;
            if (ratio <= 0) continue;
            x_val = log10(ratio);
        } else {
            x_val = dt;
        }

        double y_val = 0.0;
        if (i < m_obsRawP.size()) y_val = m_obsRawP[i];

        vt << x_val;
        vp << y_val;
    }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 0, 180), 5));
    plot->graph(0)->setName("实测压力");

    if (hasModel) {
        QVector<double> vtm, vpm;
        double p_start = 0.0;
        if (!m_obsRawP.isEmpty()) p_start = m_obsRawP.first();
        else p_start = m_settings.initialPressure;

        for(int i=0; i<tm.size(); ++i) {
            double dt = tm[i];
            if (dt < 1e-10) continue;

            double x_val;
            if (useHorner) {
                double ratio = (tp + dt) / dt;
                if (ratio <= 0) continue;
                x_val = log10(ratio);
            } else {
                x_val = dt;
            }

            double y_val;
            if (m_settings.testType == Test_Drawdown) {
                y_val = p_start - pm[i];
            } else {
                y_val = p_start + pm[i];
            }

            vtm << x_val;
            vpm << y_val;
        }

        plot->addGraph();
        plot->graph(1)->setData(vtm, vpm);
        plot->graph(1)->setPen(QPen(Qt::red, 2));
        plot->graph(1)->setName("理论压力");
    }

    if (m_hasManualPressure) {
        drawPressureFitResult();
    }

    if (useHorner) {
        plot->xAxis->setLabel("Horner 时间比 lg((tp+dt)/dt)");
        plot->yAxis->setLabel("实测压力 Pressure (MPa)");
        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        plot->xAxis->setTicker(linearTicker);
        plot->xAxis->setScaleType(QCPAxis::stLinear);
        plot->xAxis->setRangeReversed(false);

        plot->yAxis->setTicker(linearTicker);
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    } else {
        plot->xAxis->setLabel("时间 Time (h)");
        plot->yAxis->setLabel("实测压力 Pressure (MPa)");
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->xAxis->setTicker(logTicker);
        plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        plot->xAxis->setRangeReversed(false);

        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        plot->yAxis->setTicker(linearTicker);
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    }

    if (autoScale) {
        plot->rescaleAxes();
        if(useHorner) {
            plot->xAxis->setRangeUpper(plot->xAxis->range().upper);
            plot->xAxis->setRangeLower(0.0);
        }
    }
}

void FittingChart::plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotCartesian;
    plot->clearGraphs();
    Q_UNUSED(dm);

    QVector<double> vt, vp;
    for(int i=0; i<m_obsT.size(); ++i) { vt << m_obsT[i]; vp << m_obsDeltaP[i]; }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    if (hasModel) {
        QVector<double> vtm, vpm;
        for(int i=0; i<tm.size(); ++i) { vtm << tm[i]; vpm << pm[i]; }
        plot->addGraph();
        plot->graph(1)->setData(vtm, vpm);
        plot->graph(1)->setPen(QPen(Qt::red, 2));
        plot->graph(1)->setName("理论压差");
    }

    if (autoScale) plot->rescaleAxes();
}

void FittingChart::plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d)
{
    if (!m_plotLogLog) return;
    QCPGraph* gP = m_plotLogLog->addGraph();
    gP->setData(t, p);
    gP->setPen(Qt::NoPen);
    gP->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor(0, 100, 0)), QBrush(QColor(0, 100, 0)), 6));
    gP->setName("抽样压差");

    QCPGraph* gD = m_plotLogLog->addGraph();
    gD->setData(t, d);
    gD->setPen(Qt::NoPen);
    gD->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, QPen(Qt::magenta), QBrush(Qt::magenta), 6));
    gD->setName("抽样导数");
}

double FittingChart::getCalculatedInitialPressure() const { return m_calculatedPi; }

void FittingChart::showResultOnLogPlot()
{
    if(!m_plotLogLog) return;
    QCPItemText *textLabel = new QCPItemText(m_plotLogLog);
    textLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignRight);
    textLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.95, 0.05);
    textLabel->setText(QString("Horner推算Pi: %1 MPa").arg(m_calculatedPi, 0, 'f', 2));
    textLabel->setFont(QFont("Microsoft YaHei", 10, QFont::Bold));
    textLabel->setColor(Qt::red);
    textLabel->setBrush(QBrush(QColor(255, 255, 255, 200)));
    textLabel->setPadding(QMargins(5, 5, 5, 5));
    textLabel->setPen(QPen(Qt::black));
}
