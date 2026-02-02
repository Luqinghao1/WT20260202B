/*
 * 文件名: wt_fittingwidget.cpp
 * 文件作用: 试井拟合分析主界面实现文件
 * 修改记录:
 * 1. 连接 m_chartSemiLog 信号。
 * 2. 实现 onSemiLogLineMoved。
 * 3. [修复] 移除对不存在的 UI 控件 lineEdit_Equation/Pi 的访问，仅更新参数表。
 */

#include "wt_fittingwidget.h"
#include "ui_wt_fittingwidget.h"
#include "modelparameter.h"
#include "modelselect.h"
#include "fittingdatadialog.h"
#include "pressurederivativecalculator.h"
#include "pressurederivativecalculator1.h"
#include "paramselectdialog.h"
#include "fittingreport.h"
#include "fittingchart.h"

#include <QMessageBox>
#include <QDebug>
#include <cmath>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QBuffer>
#include <QFileInfo>
#include <QDateTime>

FittingWidget::FittingWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::FittingWidget),
    m_modelManager(nullptr),
    m_core(new FittingCore(this)),
    m_chartManager(new FittingChart(this)),
    m_mdiArea(nullptr),
    m_chartLogLog(nullptr), m_chartSemiLog(nullptr), m_chartCartesian(nullptr),
    m_subWinLogLog(nullptr), m_subWinSemiLog(nullptr), m_subWinCartesian(nullptr),
    m_plotLogLog(nullptr), m_plotSemiLog(nullptr), m_plotCartesian(nullptr),
    m_currentModelType(ModelManager::Model_1),
    m_isFitting(false),
    m_isCustomSamplingEnabled(false)
{
    ui->setupUi(this);

    if (ui->plotContainer->layout()) {
        QLayoutItem* item;
        while ((item = ui->plotContainer->layout()->takeAt(0)) != nullptr) {
            delete item->widget(); delete item;
        }
        delete ui->plotContainer->layout();
    }
    QVBoxLayout* containerLayout = new QVBoxLayout(ui->plotContainer);
    containerLayout->setContentsMargins(0,0,0,0);
    containerLayout->setSpacing(0);

    m_mdiArea = new QMdiArea(this);
    m_mdiArea->setViewMode(QMdiArea::SubWindowView);
    m_mdiArea->setBackground(QBrush(QColor(240, 240, 240)));
    containerLayout->addWidget(m_mdiArea);

    m_chartLogLog = new FittingChart1(this);
    m_chartSemiLog = new FittingChart2(this);
    m_chartCartesian = new FittingChart3(this);

    m_plotLogLog = m_chartLogLog->getPlot();
    m_plotSemiLog = m_chartSemiLog->getPlot();
    m_plotCartesian = m_chartCartesian->getPlot();

    m_chartLogLog->setTitle("双对数曲线 (Log-Log)");
    m_chartSemiLog->setTitle("半对数曲线 (Semi-Log)");
    m_chartCartesian->setTitle("历史拟合曲线 (History Plot)");

    m_subWinLogLog = m_mdiArea->addSubWindow(m_chartLogLog);
    m_subWinSemiLog = m_mdiArea->addSubWindow(m_chartSemiLog);
    m_subWinCartesian = m_mdiArea->addSubWindow(m_chartCartesian);

    m_subWinLogLog->setWindowTitle("双对数图");
    m_subWinSemiLog->setWindowTitle("半对数图");
    m_subWinCartesian->setWindowTitle("标准坐标系");

    connect(m_chartLogLog, &FittingChart1::exportDataTriggered, this, &FittingWidget::onExportCurveData);
    connect(m_chartSemiLog, &FittingChart2::exportDataTriggered, this, &FittingWidget::onExportCurveData);
    connect(m_chartCartesian, &FittingChart3::exportDataTriggered, this, &FittingWidget::onExportCurveData);

    // 连接信号
    connect(m_chartSemiLog, &FittingChart2::sigLineMoved, this, &FittingWidget::onSemiLogLineMoved);

    ui->splitter->setSizes(QList<int>() << 350 << 650);
    ui->splitter->setCollapsible(0, false);

    m_paramChart = new FittingParameterChart(ui->tableParams, this);

    connect(m_paramChart, &FittingParameterChart::parameterChangedByWheel, this, [this](){
        updateModelCurve(nullptr, false, false);
    });

    setupPlot();
    m_chartManager->initializeCharts(m_plotLogLog, m_plotSemiLog, m_plotCartesian);

    qRegisterMetaType<QMap<QString,double>>("QMap<QString,double>");
    qRegisterMetaType<ModelManager::ModelType>("ModelManager::ModelType");
    qRegisterMetaType<QVector<double>>("QVector<double>");

    connect(m_core, &FittingCore::sigIterationUpdated, this, &FittingWidget::onIterationUpdate, Qt::QueuedConnection);
    connect(m_core, &FittingCore::sigProgress, ui->progressBar, &QProgressBar::setValue);
    connect(m_core, &FittingCore::sigFitFinished, this, &FittingWidget::onFitFinished);

    connect(ui->sliderWeight, &QSlider::valueChanged, this, &FittingWidget::onSliderWeightChanged);
    connect(ui->btnSamplingSettings, &QPushButton::clicked, this, &FittingWidget::onOpenSamplingSettings);

    ui->sliderWeight->setRange(0, 100);
    ui->sliderWeight->setValue(50);
    onSliderWeightChanged(50);
}

FittingWidget::~FittingWidget()
{
    delete ui;
}

// 响应半对数直线移动
void FittingWidget::onSemiLogLineMoved(double k, double b)
{
    // 更新参数表中的 Pi (初始地层压力)
    QList<FitParameter> params = m_paramChart->getParameters();
    bool updated = false;
    for(auto& p : params) {
        if(p.name == "Pi" || p.name == "p*") {
            p.value = b;
            updated = true;
            break;
        }
    }
    if(updated) {
        m_paramChart->setParameters(params);
        // 如果需要，可以在这里调用 updateModelCurve 刷新理论曲线
    }

    // 方程更新逻辑 (移除，因为UI中没有对应控件)
    // 如果您后续添加了 lineEdit_Equation，可以解开注释:
    /*
    QString eq = QString("P = %1 log(t) + %2")
                    .arg(k, 0, 'f', 4)
                    .arg(b, 0, 'f', 4);
    if(ui->lineEdit_Equation) ui->lineEdit_Equation->setText(eq);
    */
}

void FittingWidget::showEvent(QShowEvent *event)
{
    QWidget::showEvent(event);
    layoutCharts();
}

void FittingWidget::setModelManager(ModelManager *m)
{
    m_modelManager = m;
    m_paramChart->setModelManager(m);
    if (m_core) m_core->setModelManager(m);
    initializeDefaultModel();
}

void FittingWidget::setProjectDataModels(const QMap<QString, QStandardItemModel *> &models)
{
    m_dataMap = models;
}

void FittingWidget::updateBasicParameters() {}

void FittingWidget::initializeDefaultModel()
{
    if(!m_modelManager) return;
    m_currentModelType = ModelManager::Model_1;
    ui->btn_modelSelect->setText("当前: " + ModelManager::getModelTypeName(m_currentModelType));
    on_btnResetParams_clicked();
}

void FittingWidget::setupPlot() {
    if(m_plotLogLog) m_plotLogLog->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    if(m_plotSemiLog) m_plotSemiLog->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    if(m_plotCartesian) m_plotCartesian->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void FittingWidget::resizeEvent(QResizeEvent *event)
{
    QWidget::resizeEvent(event);
    layoutCharts();
}

void FittingWidget::layoutCharts()
{
    if (!m_mdiArea || !m_subWinLogLog || !m_subWinSemiLog || !m_subWinCartesian) return;
    QRect rect = m_mdiArea->contentsRect();
    int w = rect.width();
    int h = rect.height();
    if (w <= 0 || h <= 0) return;

    m_subWinLogLog->setGeometry(0, 0, w / 2, h);
    m_subWinCartesian->setGeometry(w / 2, 0, w - (w / 2), h / 2);
    m_subWinSemiLog->setGeometry(w / 2, h / 2, w - (w / 2), h - (h / 2));

    if (m_subWinLogLog->isMinimized()) m_subWinLogLog->showNormal();
    if (m_subWinCartesian->isMinimized()) m_subWinCartesian->showNormal();
    if (m_subWinSemiLog->isMinimized()) m_subWinSemiLog->showNormal();
}

void FittingWidget::setObservedData(const QVector<double>& t, const QVector<double>& deltaP, const QVector<double>& d)
{
    setObservedData(t, deltaP, d, QVector<double>());
}

void FittingWidget::setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                                    const QVector<double>& d, const QVector<double>& rawP)
{
    m_obsTime = t;
    m_obsDeltaP = deltaP;
    m_obsDerivative = d;
    m_obsRawP = rawP;

    if (m_core) m_core->setObservedData(t, deltaP, d);
    if (m_chartManager) m_chartManager->setObservedData(t, deltaP, d, rawP);

    updateModelCurve(nullptr, true);
}

void FittingWidget::on_btnLoadData_clicked() {
    FittingDataDialog dlg(m_dataMap, this);
    if (dlg.exec() != QDialog::Accepted) return;

    FittingDataSettings settings = dlg.getSettings();
    QStandardItemModel* sourceModel = dlg.getPreviewModel();

    if (!sourceModel || sourceModel->rowCount() == 0) {
        QMessageBox::warning(this, "警告", "所选数据源为空，无法加载！");
        return;
    }

    QVector<double> rawTime, rawPressureData, finalDeriv;
    int skip = settings.skipRows;
    int rows = sourceModel->rowCount();

    for (int i = skip; i < rows; ++i) {
        QStandardItem* itemT = sourceModel->item(i, settings.timeColIndex);
        QStandardItem* itemP = sourceModel->item(i, settings.pressureColIndex);

        if (itemT && itemP) {
            bool okT, okP;
            double t = itemT->text().toDouble(&okT);
            double p = itemP->text().toDouble(&okP);

            if (okT && okP && t > 0) {
                rawTime.append(t);
                rawPressureData.append(p);
                if (settings.derivColIndex >= 0) {
                    QStandardItem* itemD = sourceModel->item(i, settings.derivColIndex);
                    if (itemD) finalDeriv.append(itemD->text().toDouble());
                    else finalDeriv.append(0.0);
                }
            }
        }
    }

    if (rawTime.isEmpty()) {
        QMessageBox::warning(this, "警告", "未能提取到有效数据。");
        return;
    }

    QVector<double> finalDeltaP;
    double p_shutin = rawPressureData.first();

    for (double p : rawPressureData) {
        double deltaP = 0.0;
        if (settings.testType == Test_Drawdown) {
            deltaP = std::abs(settings.initialPressure - p);
        } else {
            deltaP = std::abs(p - p_shutin);
        }
        finalDeltaP.append(deltaP);
    }

    if (settings.derivColIndex == -1) {
        finalDeriv = PressureDerivativeCalculator::calculateBourdetDerivative(rawTime, finalDeltaP, settings.lSpacing);
        if (settings.enableSmoothing) {
            finalDeriv = PressureDerivativeCalculator1::smoothData(finalDeriv, settings.smoothingSpan);
        }
    } else {
        if (settings.enableSmoothing) {
            finalDeriv = PressureDerivativeCalculator1::smoothData(finalDeriv, settings.smoothingSpan);
        }
        if (finalDeriv.size() != rawTime.size()) {
            finalDeriv.resize(rawTime.size());
        }
    }

    m_chartManager->setSettings(settings);
    setObservedData(rawTime, finalDeltaP, finalDeriv, rawPressureData);
    QMessageBox::information(this, "成功", "观测数据已成功加载。");
}

void FittingWidget::onSliderWeightChanged(int value)
{
    double wPressure = value / 100.0;
    double wDerivative = 1.0 - wPressure;
    ui->label_ValDerivative->setText(QString("导数权重: %1").arg(wDerivative, 0, 'f', 2));
    ui->label_ValPressure->setText(QString("压差权重: %1").arg(wPressure, 0, 'f', 2));
}

void FittingWidget::on_btnSelectParams_clicked()
{
    m_paramChart->updateParamsFromTable();
    QList<FitParameter> currentParams = m_paramChart->getParameters();
    ParamSelectDialog dlg(currentParams, this);
    if(dlg.exec() == QDialog::Accepted) {
        QList<FitParameter> updatedParams = dlg.getUpdatedParams();
        for(auto& p : updatedParams) {
            if(p.name == "LfD") p.isFit = false;
        }
        m_paramChart->setParameters(updatedParams);
        hideUnwantedParams();

        updateModelCurve(nullptr, false);
    }
}

void FittingWidget::hideUnwantedParams()
{
    for(int i = 0; i < ui->tableParams->rowCount(); ++i) {
        QTableWidgetItem* item = ui->tableParams->item(i, 1);
        if(item) {
            QString name = item->data(Qt::UserRole).toString();
            if(name == "LfD") {
                ui->tableParams->setRowHidden(i, true);
            }
        }
    }
}

void FittingWidget::onOpenSamplingSettings()
{
    if (m_obsTime.isEmpty()) {
        QMessageBox::warning(this, "提示", "请先加载观测数据，以便确定时间范围。");
        return;
    }
    double tMin = m_obsTime.first();
    double tMax = m_obsTime.last();

    SamplingSettingsDialog dlg(m_customIntervals, m_isCustomSamplingEnabled, tMin, tMax, this);
    if (dlg.exec() == QDialog::Accepted) {
        m_customIntervals = dlg.getIntervals();
        m_isCustomSamplingEnabled = dlg.isCustomSamplingEnabled();
        if(m_core) m_core->setSamplingSettings(m_customIntervals, m_isCustomSamplingEnabled);

        updateModelCurve(nullptr, false);
    }
}

void FittingWidget::on_btnUpdateLimits_clicked()
{
    m_paramChart->updateParamsFromTable();
    m_paramChart->autoAdjustLimits();
    QMessageBox::information(this, "提示", "参数上下限及滚轮步长已根据当前值更新。\n\n范围: 0.1倍 ~ 10倍\n步长: 范围的1/20");
}

void FittingWidget::on_btnRunFit_clicked() {
    if(m_isFitting) return;
    if(m_obsTime.isEmpty()) {
        QMessageBox::warning(this,"错误","请先加载观测数据。");
        return;
    }

    m_paramChart->updateParamsFromTable();
    m_isFitting = true;
    ui->btnRunFit->setEnabled(false);

    ui->btnSelectParams->setEnabled(false);
    ui->btnUpdateLimits->setEnabled(false);

    ModelManager::ModelType modelType = m_currentModelType;
    QList<FitParameter> paramsCopy = m_paramChart->getParameters();
    double w = ui->sliderWeight->value() / 100.0;

    if(m_core) m_core->startFit(modelType, paramsCopy, w);
}

void FittingWidget::on_btnStop_clicked() {
    if(m_core) m_core->stopFit();
}

void FittingWidget::on_btnImportModel_clicked() {
    updateModelCurve(nullptr, false, false);
}

void FittingWidget::on_btnResetParams_clicked() {
    if(!m_modelManager) return;
    m_paramChart->resetParams(m_currentModelType, true);
    loadProjectParams();
    hideUnwantedParams();
    updateModelCurve(nullptr, true, true);
}

void FittingWidget::on_btn_modelSelect_clicked() {
    ModelSelect dlg(this);
    if (dlg.exec() == QDialog::Accepted) {
        QString code = dlg.getSelectedModelCode();
        QString name = dlg.getSelectedModelName();
        bool found = false;
        ModelManager::ModelType newType = ModelManager::Model_1;

        if (code == "modelwidget1") newType = ModelManager::Model_1;
        else if (code == "modelwidget2") newType = ModelManager::Model_2;
        else if (code == "modelwidget3") newType = ModelManager::Model_3;
        else if (code == "modelwidget4") newType = ModelManager::Model_4;
        else if (code == "modelwidget5") newType = ModelManager::Model_5;
        else if (code == "modelwidget6") newType = ModelManager::Model_6;
        else if (!code.isEmpty()) found = true;

        if (code.startsWith("modelwidget")) found = true;

        if (found) {
            m_paramChart->switchModel(newType);
            m_currentModelType = newType;
            ui->btn_modelSelect->setText("当前: " + name);
            loadProjectParams();
            hideUnwantedParams();
            updateModelCurve(nullptr, true);
        } else {
            QMessageBox::warning(this, "提示", "所选组合暂无对应的模型。\nCode: " + code);
        }
    }
}

void FittingWidget::loadProjectParams()
{
    ModelParameter* mp = ModelParameter::instance();
    QList<FitParameter> params = m_paramChart->getParameters();
    bool changed = false;
    for(auto& p : params) {
        if(p.name == "phi") { p.value = mp->getPhi(); changed = true; }
        else if(p.name == "h") { p.value = mp->getH(); changed = true; }
        else if(p.name == "rw") { p.value = mp->getRw(); changed = true; }
        else if(p.name == "mu") { p.value = mp->getMu(); changed = true; }
        else if(p.name == "Ct") { p.value = mp->getCt(); changed = true; }
        else if(p.name == "B") { p.value = mp->getB(); changed = true; }
        else if(p.name == "q") { p.value = mp->getQ(); changed = true; }
    }
    if(changed) m_paramChart->setParameters(params);
}

void FittingWidget::updateModelCurve(const QMap<QString, double>* explicitParams, bool autoScale, bool calcError) {
    if(!m_modelManager) {
        QMessageBox::critical(this, "错误", "ModelManager 未初始化！");
        return;
    }

    if(m_obsTime.isEmpty()) {
        m_chartLogLog->clearGraphs();
        m_chartSemiLog->clearGraphs();
        m_chartCartesian->clearGraphs();
        return;
    }

    ui->tableParams->clearFocus();

    QMap<QString, double> rawParams;
    QString sensitivityKey = "";
    QVector<double> sensitivityValues;

    if (explicitParams) {
        rawParams = *explicitParams;
    } else {
        QList<FitParameter> allParams = m_paramChart->getParameters();
        for(const auto& p : allParams) rawParams.insert(p.name, p.value);

        QMap<QString, QString> rawTexts = m_paramChart->getRawParamTexts();
        for(auto it = rawTexts.begin(); it != rawTexts.end(); ++it) {
            QVector<double> vals = parseSensitivityValues(it.value());
            if (!vals.isEmpty()) {
                rawParams.insert(it.key(), vals.first());
                if (vals.size() > 1 && sensitivityKey.isEmpty()) {
                    sensitivityKey = it.key();
                    sensitivityValues = vals;
                }
            } else {
                rawParams.insert(it.key(), 0.0);
            }
        }
    }

    QMap<QString, double> solverParams = FittingCore::preprocessParams(rawParams, m_currentModelType);

    QVector<double> targetT;
    if (m_obsTime.size() > 300) {
        double tMin = m_obsTime.first() > 1e-5 ? m_obsTime.first() : 1e-5;
        double tMax = m_obsTime.last();
        targetT = ModelManager::generateLogTimeSteps(300, log10(tMin), log10(tMax));
    } else if (!m_obsTime.isEmpty()) {
        targetT = m_obsTime;
    } else {
        for(double e = -4; e <= 4; e += 0.1) targetT.append(pow(10, e));
    }

    bool isSensitivityMode = !sensitivityKey.isEmpty();
    ui->btnRunFit->setEnabled(!isSensitivityMode);

    if (isSensitivityMode) {
        ui->label_Error->setText(QString("敏感性分析模式: %1 (%2 个值)").arg(sensitivityKey).arg(sensitivityValues.size()));
        m_chartLogLog->clearGraphs();

        m_chartManager->plotAll(QVector<double>(), QVector<double>(), QVector<double>(), false, autoScale);

        QList<QColor> colors = { Qt::red, Qt::blue, QColor(0,180,0), Qt::magenta, QColor(255,140,0), Qt::cyan, Qt::darkRed, Qt::darkBlue };
        for(int i = 0; i < sensitivityValues.size(); ++i) {
            double val = sensitivityValues[i];

            QMap<QString, double> currentParams = rawParams;
            currentParams[sensitivityKey] = val;
            QMap<QString, double> currentSolverParams = FittingCore::preprocessParams(currentParams, m_currentModelType);

            ModelCurveData res = m_modelManager->calculateTheoreticalCurve(m_currentModelType, currentSolverParams, targetT);
            QColor c = colors[i % colors.size()];
            QString suffix = QString("%1=%2").arg(sensitivityKey).arg(val);

            QCPGraph* gP = m_plotLogLog->addGraph();
            gP->setData(std::get<0>(res), std::get<1>(res));
            gP->setPen(QPen(c, 2)); gP->setName("P: "+suffix);

            QCPGraph* gD = m_plotLogLog->addGraph();
            gD->setData(std::get<0>(res), std::get<2>(res));
            gD->setPen(QPen(c, 2, Qt::DashLine)); gD->setName("P': "+suffix);
        }

    } else {
        ModelCurveData res = m_modelManager->calculateTheoreticalCurve(m_currentModelType, solverParams, targetT);

        m_chartManager->plotAll(std::get<0>(res), std::get<1>(res), std::get<2>(res), true, autoScale);

        if (!m_obsTime.isEmpty() && m_core && calcError) {
            QVector<double> sampleT, sampleP, sampleD;
            m_core->getLogSampledData(m_obsTime, m_obsDeltaP, m_obsDerivative, sampleT, sampleP, sampleD);

            QVector<double> residuals = m_core->calculateResiduals(rawParams, m_currentModelType, ui->sliderWeight->value()/100.0, sampleT, sampleP, sampleD);
            double sse = m_core->calculateSumSquaredError(residuals);
            ui->label_Error->setText(QString("误差(MSE): %1").arg(sse/residuals.size(), 0, 'e', 3));

            if (m_isCustomSamplingEnabled) {
                m_chartManager->plotSampledPoints(sampleT, sampleP, sampleD);
            }
        }
    }

    m_plotLogLog->replot();
    m_plotSemiLog->replot();
    m_plotCartesian->replot();
}

void FittingWidget::onIterationUpdate(double err, const QMap<QString,double>& p,
                                      const QVector<double>& t, const QVector<double>& p_curve, const QVector<double>& d_curve) {
    ui->label_Error->setText(QString("误差(MSE): %1").arg(err, 0, 'e', 3));

    ui->tableParams->blockSignals(true);
    for(int i=0; i<ui->tableParams->rowCount(); ++i) {
        QString key = ui->tableParams->item(i, 1)->data(Qt::UserRole).toString();
        if(p.contains(key)) {
            double val = p[key];
            ui->tableParams->item(i, 2)->setText(QString::number(val, 'g', 5));
        }
    }
    ui->tableParams->blockSignals(false);

    m_chartManager->plotAll(t, p_curve, d_curve, true, false);

    if (m_isCustomSamplingEnabled && m_core) {
        QVector<double> sampleT, sampleP, sampleD;
        m_core->getLogSampledData(m_obsTime, m_obsDeltaP, m_obsDerivative, sampleT, sampleP, sampleD);
        m_chartManager->plotSampledPoints(sampleT, sampleP, sampleD);
    }

    if(m_plotLogLog) m_plotLogLog->replot();
    if(m_plotSemiLog) m_plotSemiLog->replot();
    if(m_plotCartesian) m_plotCartesian->replot();
}

void FittingWidget::onFitFinished() {
    m_isFitting = false;
    ui->btnRunFit->setEnabled(true);

    ui->btnSelectParams->setEnabled(true);
    ui->btnUpdateLimits->setEnabled(true);

    QMessageBox::information(this, "完成", "拟合完成。");
}

void FittingWidget::on_btnExportData_clicked() {
    m_paramChart->updateParamsFromTable();
    QList<FitParameter> params = m_paramChart->getParameters();

    QString defaultDir = ModelParameter::instance()->getProjectPath();
    if(defaultDir.isEmpty()) defaultDir = ".";

    QString fileName = QFileDialog::getSaveFileName(this, "导出拟合参数", defaultDir + "/FittingParameters.csv", "CSV Files (*.csv);;Text Files (*.txt)");
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
    QTextStream out(&file);

    if(fileName.endsWith(".csv", Qt::CaseInsensitive)) {
        file.write("\xEF\xBB\xBF");
        out << QString("参数中文名,参数英文名,拟合值,单位\n");
        for(const auto& param : params) {
            QString htmlSym, uniSym, unitStr, dummyName;
            FittingParameterChart::getParamDisplayInfo(param.name, dummyName, htmlSym, uniSym, unitStr);
            if(unitStr == "无因次" || unitStr == "小数") unitStr = "";
            out << QString("%1,%2,%3,%4\n").arg(param.displayName).arg(uniSym).arg(param.value, 0, 'g', 10).arg(unitStr);
        }
    } else {
        for(const auto& param : params) {
            QString htmlSym, uniSym, unitStr, dummyName;
            FittingParameterChart::getParamDisplayInfo(param.name, dummyName, htmlSym, uniSym, unitStr);
            if(unitStr == "无因次" || unitStr == "小数") unitStr = "";
            out << QString("%1 (%2): %3 %4").arg(param.displayName).arg(uniSym).arg(param.value, 0, 'g', 10).arg(unitStr) << "\n";
        }
    }
    file.close();
    QMessageBox::information(this, "完成", "参数数据已成功导出。");
}

void FittingWidget::onExportCurveData() {
    QString defaultDir = ModelParameter::instance()->getProjectPath();
    if(defaultDir.isEmpty()) defaultDir = ".";

    QString path = QFileDialog::getSaveFileName(this, "导出拟合曲线数据", defaultDir + "/FittingCurves.csv", "CSV Files (*.csv)");
    if (path.isEmpty()) return;

    auto graphObsP = m_plotLogLog->graph(0);
    auto graphObsD = m_plotLogLog->graph(1);

    if (!graphObsP) return;

    QCPGraph *graphModP = (m_plotLogLog->graphCount() > 2) ? m_plotLogLog->graph(2) : nullptr;
    QCPGraph *graphModD = (m_plotLogLog->graphCount() > 3) ? m_plotLogLog->graph(3) : nullptr;

    QFile f(path);
    if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&f);
        out << "Obs_Time,Obs_DP,Obs_Deriv,Model_Time,Model_DP,Model_Deriv\n";

        auto itObsP = graphObsP->data()->begin();
        auto itObsD = graphObsD->data()->begin();
        auto endObsP = graphObsP->data()->end();

        QCPGraphDataContainer::const_iterator itModP, endModP, itModD;
        bool hasModel = (graphModP != nullptr && graphModD != nullptr);
        if(hasModel) {
            itModP = graphModP->data()->begin();
            endModP = graphModP->data()->end();
            itModD = graphModD->data()->begin();
        }

        while (itObsP != endObsP || (hasModel && itModP != endModP)) {
            QStringList line;
            if (itObsP != endObsP) {
                line << QString::number(itObsP->key, 'g', 10);
                line << QString::number(itObsP->value, 'g', 10);
                if (itObsD != graphObsD->data()->end()) {
                    line << QString::number(itObsD->value, 'g', 10);
                    ++itObsD;
                } else { line << ""; }
                ++itObsP;
            } else { line << "" << "" << ""; }

            if (hasModel && itModP != endModP) {
                line << QString::number(itModP->key, 'g', 10);
                line << QString::number(itModP->value, 'g', 10);
                if (itModD != graphModD->data()->end()) {
                    line << QString::number(itModD->value, 'g', 10);
                    ++itModD;
                } else { line << ""; }
                ++itModP;
            } else { line << "" << "" << ""; }
            out << line.join(",") << "\n";
        }
        f.close();
        QMessageBox::information(this, "导出成功", "拟合曲线数据已保存。");
    }
}

void FittingWidget::on_btnExportReport_clicked()
{
    QString wellName = "未命名井";
    QString projectFilePath = ModelParameter::instance()->getProjectFilePath();
    QFile pwtFile(projectFilePath);
    if (pwtFile.exists() && pwtFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QByteArray data = pwtFile.readAll();
        QJsonDocument doc = QJsonDocument::fromJson(data);
        if (!doc.isNull() && doc.isObject()) {
            QJsonObject root = doc.object();
            if (root.contains("wellName")) wellName = root["wellName"].toString();
            else if (root.contains("basicParams")) {
                QJsonObject basic = root["basicParams"].toObject();
                if (basic.contains("wellName")) wellName = basic["wellName"].toString();
            }
        }
        pwtFile.close();
    }
    if (wellName == "未命名井" || wellName.isEmpty()) {
        QFileInfo fi(projectFilePath);
        wellName = fi.completeBaseName();
    }

    FittingReportData reportData;
    reportData.wellName = wellName;
    reportData.modelType = m_currentModelType;

    QString mseText = ui->label_Error->text().remove("误差(MSE): ");
    reportData.mse = mseText.toDouble();

    reportData.t = m_obsTime;
    reportData.p = m_obsDeltaP;
    reportData.d = m_obsDerivative;

    m_paramChart->updateParamsFromTable();
    reportData.params = m_paramChart->getParameters();

    reportData.imgLogLog = getPlotImageBase64(m_plotLogLog);
    reportData.imgSemiLog = getPlotImageBase64(m_plotSemiLog);
    reportData.imgCartesian = getPlotImageBase64(m_plotCartesian);

    QString reportFileName = QString("%1试井解释报告.doc").arg(wellName);
    QString defaultDir = QFileInfo(projectFilePath).absolutePath();
    if(defaultDir.isEmpty() || defaultDir == ".") defaultDir = ModelParameter::instance()->getProjectPath();
    if(defaultDir.isEmpty()) defaultDir = ".";

    QString fileName = QFileDialog::getSaveFileName(this, "导出报告", defaultDir + "/" + reportFileName, "Word 文档 (*.doc);;HTML 文件 (*.html)");
    if(fileName.isEmpty()) return;

    QString errorMsg;
    if (FittingReportGenerator::generate(fileName, reportData, &errorMsg)) {
        QMessageBox::information(this, "成功", QString("报告及数据已导出！\n\n文件路径: %1").arg(fileName));
    } else {
        QMessageBox::critical(this, "错误", "报告导出失败:\n" + errorMsg);
    }
}

QString FittingWidget::getPlotImageBase64(MouseZoom* plot) {
    if(!plot) return "";
    QPixmap pixmap = plot->toPixmap(800, 600);
    QByteArray byteArray;
    QBuffer buffer(&byteArray);
    buffer.open(QIODevice::WriteOnly);
    pixmap.save(&buffer, "PNG");
    return QString::fromLatin1(byteArray.toBase64().data());
}

void FittingWidget::on_btnSaveFit_clicked()
{
    emit sigRequestSave();
}

QJsonObject FittingWidget::getJsonState() const
{
    const_cast<FittingWidget*>(this)->m_paramChart->updateParamsFromTable();
    QList<FitParameter> params = m_paramChart->getParameters();

    QJsonObject root;
    root["modelType"] = (int)m_currentModelType;
    root["modelName"] = ModelManager::getModelTypeName(m_currentModelType);
    root["fitWeightVal"] = ui->sliderWeight->value();

    QJsonObject plotRange;
    plotRange["xMin"] = m_plotLogLog->xAxis->range().lower;
    plotRange["xMax"] = m_plotLogLog->xAxis->range().upper;
    plotRange["yMin"] = m_plotLogLog->yAxis->range().lower;
    plotRange["yMax"] = m_plotLogLog->yAxis->range().upper;
    root["plotView"] = plotRange;

    QJsonArray paramsArray;
    for(const auto& p : params) {
        QJsonObject pObj;
        pObj["name"] = p.name;
        pObj["value"] = p.value;
        pObj["isFit"] = p.isFit;
        pObj["min"] = p.min;
        pObj["max"] = p.max;
        pObj["isVisible"] = p.isVisible;
        pObj["step"] = p.step;
        paramsArray.append(pObj);
    }
    root["parameters"] = paramsArray;

    QJsonArray timeArr, pressArr, derivArr, rawPArr;
    for(double v : m_obsTime) timeArr.append(v);
    for(double v : m_obsDeltaP) pressArr.append(v);
    for(double v : m_obsDerivative) derivArr.append(v);
    for(double v : m_obsRawP) rawPArr.append(v);

    QJsonObject obsData;
    obsData["time"] = timeArr;
    obsData["pressure"] = pressArr;
    obsData["derivative"] = derivArr;
    obsData["rawPressure"] = rawPArr;
    root["observedData"] = obsData;

    root["useCustomSampling"] = m_isCustomSamplingEnabled;
    QJsonArray intervalArr;
    for(const auto& item : m_customIntervals) {
        QJsonObject obj;
        obj["start"] = item.tStart;
        obj["end"] = item.tEnd;
        obj["count"] = item.count;
        intervalArr.append(obj);
    }
    root["customIntervals"] = intervalArr;

    return root;
}

void FittingWidget::loadFittingState(const QJsonObject& root)
{
    if (root.isEmpty()) return;

    if (root.contains("modelType")) {
        int type = root["modelType"].toInt();
        m_currentModelType = (ModelManager::ModelType)type;
        ui->btn_modelSelect->setText("当前: " + ModelManager::getModelTypeName(m_currentModelType));
    }

    m_paramChart->resetParams(m_currentModelType);

    QMap<QString, double> explicitParamsMap;
    if (root.contains("parameters")) {
        QJsonArray arr = root["parameters"].toArray();
        QList<FitParameter> currentParams = m_paramChart->getParameters();
        for(int i=0; i<arr.size(); ++i) {
            QJsonObject pObj = arr[i].toObject();
            QString name = pObj["name"].toString();
            for(auto& p : currentParams) {
                if(p.name == name) {
                    p.value = pObj["value"].toDouble();
                    p.isFit = pObj["isFit"].toBool();
                    p.min = pObj["min"].toDouble();
                    p.max = pObj["max"].toDouble();
                    if(pObj.contains("isVisible")) p.isVisible = pObj["isVisible"].toBool();
                    else p.isVisible = true;
                    if(pObj.contains("step")) p.step = pObj["step"].toDouble();
                    explicitParamsMap.insert(p.name, p.value);
                    break;
                }
            }
        }
        m_paramChart->setParameters(currentParams);
    }

    if (root.contains("fitWeightVal")) ui->sliderWeight->setValue(root["fitWeightVal"].toInt());

    if (root.contains("observedData")) {
        QJsonObject obs = root["observedData"].toObject();
        QJsonArray tArr = obs["time"].toArray();
        QJsonArray pArr = obs["pressure"].toArray();
        QJsonArray dArr = obs["derivative"].toArray();
        QJsonArray rawPArr;
        if (obs.contains("rawPressure")) rawPArr = obs["rawPressure"].toArray();

        QVector<double> t, p, d, rawP;
        for(auto v : tArr) t.append(v.toDouble());
        for(auto v : pArr) p.append(v.toDouble());
        for(auto v : dArr) d.append(v.toDouble());
        for(auto v : rawPArr) rawP.append(v.toDouble());

        setObservedData(t, p, d, rawP);
    }

    if (root.contains("useCustomSampling")) m_isCustomSamplingEnabled = root["useCustomSampling"].toBool();
    if (root.contains("customIntervals")) {
        m_customIntervals.clear();
        QJsonArray arr = root["customIntervals"].toArray();
        for(auto v : arr) {
            QJsonObject obj = v.toObject();
            SamplingInterval item;
            item.tStart = obj["start"].toDouble();
            item.tEnd = obj["end"].toDouble();
            item.count = obj["count"].toInt();
            m_customIntervals.append(item);
        }
        if(m_core) m_core->setSamplingSettings(m_customIntervals, m_isCustomSamplingEnabled);
    }

    hideUnwantedParams();

    updateModelCurve(&explicitParamsMap);

    if (root.contains("plotView")) {
        QJsonObject range = root["plotView"].toObject();
        if (range.contains("xMin") && range.contains("xMax")) {
            double xMin = range["xMin"].toDouble();
            double xMax = range["xMax"].toDouble();
            double yMin = range["yMin"].toDouble();
            double yMax = range["yMax"].toDouble();
            if (xMax > xMin && yMax > yMin && xMin > 0 && yMin > 0) {
                m_plotLogLog->xAxis->setRange(xMin, xMax);
                m_plotLogLog->yAxis->setRange(yMin, yMax);
                m_plotLogLog->replot();
            }
        }
    }
}

QVector<double> FittingWidget::parseSensitivityValues(const QString& text) {
    QVector<double> values;
    QString cleanText = text;
    cleanText.replace(QChar(0xFF0C), ",");
    QStringList parts = cleanText.split(',', Qt::SkipEmptyParts);
    for (const QString& part : parts) {
        bool ok;
        double v = part.trimmed().toDouble(&ok);
        if (ok) values.append(v);
    }
    return values;
}
