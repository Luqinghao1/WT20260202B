/*
 * 文件名: wt_modelwidget.cpp
 * 文件作用: 压裂水平井复合页岩油模型界面类实现
 */

#include "wt_modelwidget.h"
#include "ui_wt_modelwidget.h"
#include "modelmanager.h"
#include "modelparameter.h"

#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QDateTime>
#include <QCoreApplication>
#include <QSplitter>
#include <QLabel>
#include <QLineEdit>
#include <QGridLayout>
#include <cmath>

WT_ModelWidget::WT_ModelWidget(ModelType type, QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::WT_ModelWidget)
    , m_type(type)
    , m_highPrecision(true)
{
    ui->setupUi(this);

    // 初始化求解器
    m_solver = new ModelSolver01_06(m_type);

    // 初始化曲线颜色列表
    m_colorList = { Qt::red, Qt::blue, QColor(0,180,0), Qt::magenta, QColor(255,140,0), Qt::cyan };

    // 设置 Splitter 初始比例 (左 20% : 右 80%)
    QList<int> sizes;
    sizes << 240 << 960;
    ui->splitter->setSizes(sizes);
    ui->splitter->setCollapsible(0, false); // 左侧不可折叠

    // 设置当前模型名称到选择按钮
    ui->btnSelectModel->setText(getModelName() + "  (点击切换)");

    initUi();
    initChart();
    setupConnections();

    // 初始化参数显示
    onResetParameters();
}

WT_ModelWidget::~WT_ModelWidget()
{
    delete m_solver; // 清理求解器资源
    delete ui;
}

QString WT_ModelWidget::getModelName() const {
    return ModelSolver01_06::getModelName(m_type);
}

// 转发给 Solver 进行计算
WT_ModelWidget::ModelCurveData WT_ModelWidget::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    if (m_solver) {
        return m_solver->calculateTheoreticalCurve(params, providedTime);
    }
    return ModelCurveData();
}

void WT_ModelWidget::setHighPrecision(bool high)
{
    m_highPrecision = high;
    if (m_solver) m_solver->setHighPrecision(high);
}

void WT_ModelWidget::initUi() {
    using MT = ModelSolver01_06::ModelType;

    // 1. 标签文本更新 (复用旧控件名，但在界面上显示新含义)
    if(ui->label_km) ui->label_km->setText("流度比 M12"); // kmEdit -> M12
    if(ui->label_rmD) ui->label_rmD->setText("复合半径 rm (m)"); // rmDEdit -> rm
    if(ui->label_reD) ui->label_reD->setText("外区半径 re (m)"); // reDEdit -> re
    if(ui->label_cD) ui->label_cD->setText("井筒储集 C (m³/MPa)"); // cDEdit -> C (有因次)

    // 2. 控制边界参数可见性 (无限大模型隐藏外边界)
    bool isInfinite = (m_type == MT::Model_1 || m_type == MT::Model_2);
    if (isInfinite) {
        ui->label_reD->setVisible(false);
        ui->reDEdit->setVisible(false);
    } else {
        ui->label_reD->setVisible(true);
        ui->reDEdit->setVisible(true);
    }

    // 3. 控制井储表皮可见性
    bool hasStorage = (m_type == MT::Model_1 || m_type == MT::Model_3 || m_type == MT::Model_5);
    ui->label_cD->setVisible(hasStorage);
    ui->cDEdit->setVisible(hasStorage);
    ui->label_s->setVisible(hasStorage);
    ui->sEdit->setVisible(hasStorage);

    // 4. 动态添加新参数输入框 (remda2, eta12)
    QWidget* parentWidget = ui->remda1Edit->parentWidget();
    QGridLayout* layout = qobject_cast<QGridLayout*>(parentWidget->layout());

    if (layout) {
        if (!parentWidget->findChild<QLineEdit*>("remda2Edit")) {
            QLabel* labelRemda2 = new QLabel("外区窜流系数 λ<sub>2</sub>:", parentWidget);
            QLineEdit* editRemda2 = new QLineEdit(parentWidget);
            editRemda2->setObjectName("remda2Edit");

            QLabel* labelEta12 = new QLabel("导压系数比 η<sub>12</sub>:", parentWidget);
            QLineEdit* editEta12 = new QLineEdit(parentWidget);
            editEta12->setObjectName("eta12Edit");

            int row = layout->rowCount();
            layout->addWidget(labelRemda2, row, 0);
            layout->addWidget(editRemda2, row, 1);

            layout->addWidget(labelEta12, row + 1, 0);
            layout->addWidget(editEta12, row + 1, 1);
        }
    }
}

void WT_ModelWidget::initChart() {
    MouseZoom* plot = ui->chartWidget->getPlot();

    plot->setBackground(Qt::white);
    plot->axisRect()->setBackground(Qt::white);

    // 设置对数坐标轴
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic); plot->xAxis->setTicker(logTicker);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic); plot->yAxis->setTicker(logTicker);
    plot->xAxis->setNumberFormat("eb"); plot->xAxis->setNumberPrecision(0);
    plot->yAxis->setNumberFormat("eb"); plot->yAxis->setNumberPrecision(0);

    QFont labelFont("Microsoft YaHei", 10, QFont::Bold);
    QFont tickFont("Microsoft YaHei", 9);
    plot->xAxis->setLabel("时间 Time (h)");
    plot->yAxis->setLabel("压力 & 导数 Pressure & Derivative (MPa)");
    plot->xAxis->setLabelFont(labelFont); plot->yAxis->setLabelFont(labelFont);
    plot->xAxis->setTickLabelFont(tickFont); plot->yAxis->setTickLabelFont(tickFont);

    plot->xAxis2->setVisible(true); plot->yAxis2->setVisible(true);
    plot->xAxis2->setTickLabels(false); plot->yAxis2->setTickLabels(false);
    connect(plot->xAxis, SIGNAL(rangeChanged(QCPRange)), plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(plot->yAxis, SIGNAL(rangeChanged(QCPRange)), plot->yAxis2, SLOT(setRange(QCPRange)));
    plot->xAxis2->setScaleType(QCPAxis::stLogarithmic); plot->yAxis2->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis2->setTicker(logTicker); plot->yAxis2->setTicker(logTicker);

    plot->xAxis->grid()->setVisible(true); plot->yAxis->grid()->setVisible(true);
    plot->xAxis->grid()->setSubGridVisible(true); plot->yAxis->grid()->setSubGridVisible(true);
    plot->xAxis->grid()->setPen(QPen(QColor(220, 220, 220), 1, Qt::SolidLine));
    plot->yAxis->grid()->setPen(QPen(QColor(220, 220, 220), 1, Qt::SolidLine));
    plot->xAxis->grid()->setSubGridPen(QPen(QColor(240, 240, 240), 1, Qt::DotLine));
    plot->yAxis->grid()->setSubGridPen(QPen(QColor(240, 240, 240), 1, Qt::DotLine));

    plot->xAxis->setRange(1e-3, 1e3); plot->yAxis->setRange(1e-3, 1e2);

    plot->legend->setVisible(true);
    plot->legend->setFont(QFont("Microsoft YaHei", 9));
    plot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

    ui->chartWidget->setTitle("复合页岩油储层试井曲线");
}

void WT_ModelWidget::setupConnections() {
    connect(ui->calculateButton, &QPushButton::clicked, this, &WT_ModelWidget::onCalculateClicked);
    connect(ui->resetButton, &QPushButton::clicked, this, &WT_ModelWidget::onResetParameters);
    connect(ui->chartWidget, &ChartWidget::exportDataTriggered, this, &WT_ModelWidget::onExportData);
    connect(ui->btnExportDataTab, &QPushButton::clicked, this, &WT_ModelWidget::onExportData);

    connect(ui->checkShowPoints, &QCheckBox::toggled, this, &WT_ModelWidget::onShowPointsToggled);

    // 转发模型选择按钮信号
    connect(ui->btnSelectModel, &QPushButton::clicked, this, &WT_ModelWidget::requestModelSelection);
}

QVector<double> WT_ModelWidget::parseInput(const QString& text) {
    QVector<double> values;
    QString cleanText = text;
    cleanText.replace("，", ",");
    QStringList parts = cleanText.split(",", Qt::SkipEmptyParts);
    for(const QString& part : parts) {
        bool ok;
        double v = part.trimmed().toDouble(&ok);
        if(ok) values.append(v);
    }
    if(values.isEmpty()) values.append(0.0);
    return values;
}

void WT_ModelWidget::setInputText(QLineEdit* edit, double value) {
    if(!edit) return;
    edit->setText(QString::number(value, 'g', 8));
}

// 重置参数函数
void WT_ModelWidget::onResetParameters() {
    using MT = ModelSolver01_06::ModelType;
    ModelParameter* mp = ModelParameter::instance();

    // 1. 基础参数
    setInputText(ui->phiEdit, mp->getPhi());
    setInputText(ui->hEdit, mp->getH());
    setInputText(ui->rwEdit, 0.1); // [新增] 井筒半径默认 0.1m
    setInputText(ui->muEdit, mp->getMu());
    setInputText(ui->BEdit, mp->getB());
    setInputText(ui->CtEdit, mp->getCt());
    setInputText(ui->qEdit, mp->getQ());

    setInputText(ui->tEdit, 1000.0);
    setInputText(ui->pointsEdit, 100);

    // 2. 模型参数默认值
    setInputText(ui->kfEdit, 1e-2); // 内区渗透率
    setInputText(ui->kmEdit, 10.0); // 流度比 M12

    double valL = 1000.0;
    setInputText(ui->LEdit, valL);  // 水平井长 L
    setInputText(ui->LfEdit, 20.0); // 裂缝半长 Lf
    setInputText(ui->nfEdit, 4);    // 裂缝条数

    // [逻辑更新] 复合半径 rm 默认值设为与 L 一致
    setInputText(ui->rmDEdit, valL);

    setInputText(ui->omga1Edit, 0.4);
    setInputText(ui->omga2Edit, 0.08);
    setInputText(ui->remda1Edit, 1e-3);

    // 动态控件参数
    QLineEdit* editRemda2 = this->findChild<QLineEdit*>("remda2Edit");
    if(editRemda2) setInputText(editRemda2, 1e-4);

    QLineEdit* editEta12 = this->findChild<QLineEdit*>("eta12Edit");
    if(editEta12) setInputText(editEta12, 0.2);

    setInputText(ui->gamaDEdit, 0.02);

    // 边界参数
    bool isInfinite = (m_type == MT::Model_1 || m_type == MT::Model_2);
    if (!isInfinite) {
        setInputText(ui->reDEdit, 20000.0); // 外边界 re
    }

    // 井储表皮
    bool hasStorage = (m_type == MT::Model_1 || m_type == MT::Model_3 || m_type == MT::Model_5);
    if (hasStorage) {
        setInputText(ui->cDEdit, 0.1); // 有因次 C (m3/MPa)
        setInputText(ui->sEdit, 0.01);
    }
}

void WT_ModelWidget::onDependentParamsChanged() {
    // 控件已移除，保留空实现
}

void WT_ModelWidget::onShowPointsToggled(bool checked) {
    MouseZoom* plot = ui->chartWidget->getPlot();
    for(int i = 0; i < plot->graphCount(); ++i) {
        if (checked) plot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
        else plot->graph(i)->setScatterStyle(QCPScatterStyle::ssNone);
    }
    plot->replot();
}

void WT_ModelWidget::onCalculateClicked() {
    ui->calculateButton->setEnabled(false);
    ui->calculateButton->setText("计算中...");
    QCoreApplication::processEvents();
    runCalculation();
    ui->calculateButton->setEnabled(true);
    ui->calculateButton->setText("开始计算");
}

void WT_ModelWidget::runCalculation() {
    MouseZoom* plot = ui->chartWidget->getPlot();
    plot->clearGraphs();

    // 收集界面输入参数
    QMap<QString, QVector<double>> rawParams;

    // 基础参数
    rawParams["phi"] = parseInput(ui->phiEdit->text());
    rawParams["h"] = parseInput(ui->hEdit->text());
    rawParams["rw"] = parseInput(ui->rwEdit->text()); // [新增] 读取 rw
    rawParams["mu"] = parseInput(ui->muEdit->text());
    rawParams["B"] = parseInput(ui->BEdit->text());
    rawParams["Ct"] = parseInput(ui->CtEdit->text());
    rawParams["q"] = parseInput(ui->qEdit->text());
    rawParams["t"] = parseInput(ui->tEdit->text());

    // 模型参数
    rawParams["kf"] = parseInput(ui->kfEdit->text());
    rawParams["M12"] = parseInput(ui->kmEdit->text());
    rawParams["L"] = parseInput(ui->LEdit->text());
    rawParams["Lf"] = parseInput(ui->LfEdit->text());
    rawParams["nf"] = parseInput(ui->nfEdit->text());
    rawParams["rm"] = parseInput(ui->rmDEdit->text());
    rawParams["omega1"] = parseInput(ui->omga1Edit->text());
    rawParams["omega2"] = parseInput(ui->omga2Edit->text());
    rawParams["lambda1"] = parseInput(ui->remda1Edit->text());
    rawParams["gamaD"] = parseInput(ui->gamaDEdit->text());

    QLineEdit* editRemda2 = this->findChild<QLineEdit*>("remda2Edit");
    if(editRemda2) rawParams["lambda2"] = parseInput(editRemda2->text());
    else rawParams["lambda2"] = {1e-4};

    QLineEdit* editEta12 = this->findChild<QLineEdit*>("eta12Edit");
    if(editEta12) rawParams["eta12"] = parseInput(editEta12->text());
    else rawParams["eta12"] = {0.2};

    if (ui->reDEdit->isVisible()) rawParams["re"] = parseInput(ui->reDEdit->text());
    else rawParams["re"] = {20000.0};

    // 井储参数处理 (C -> cD 转换)
    if (ui->cDEdit->isVisible()) {
        QVector<double> C_vals = parseInput(ui->cDEdit->text());
        QVector<double> cD_vals;

        // [重要修正] CD = 0.159 * C / (phi * h * Ct * L^2)
        // 使用界面输入的 L

        double phi = rawParams["phi"].isEmpty() ? 0.05 : rawParams["phi"].first();
        double Ct = rawParams["Ct"].isEmpty() ? 5e-4 : rawParams["Ct"].first();
        double h = rawParams["h"].isEmpty() ? 20.0 : rawParams["h"].first();
        double L = rawParams["L"].isEmpty() ? 1000.0 : rawParams["L"].first(); // [修正] 使用 L

        double denom = phi * h * Ct * L * L; // [修正] L^2
        double factor = 0.0;
        if (denom > 1e-20) factor = 0.159 / denom;

        for(double valC : C_vals) {
            cD_vals.append(valC * factor);
        }

        rawParams["cD"] = cD_vals;
        rawParams["S"] = parseInput(ui->sEdit->text());
    } else {
        rawParams["cD"] = {0.0};
        rawParams["S"] = {0.0};
    }

    // 敏感性分析检查
    QString sensitivityKey = "";
    QVector<double> sensitivityValues;
    for(auto it = rawParams.begin(); it != rawParams.end(); ++it) {
        if(it.key() == "t") continue;
        if(it.value().size() > 1) {
            sensitivityKey = it.key();
            sensitivityValues = it.value();
            break;
        }
    }
    bool isSensitivity = !sensitivityKey.isEmpty();

    // 构建基础参数字典 (取第一个值)
    QMap<QString, double> baseParams;
    for(auto it = rawParams.begin(); it != rawParams.end(); ++it) {
        baseParams[it.key()] = it.value().isEmpty() ? 0.0 : it.value().first();
    }

    baseParams["N"] = m_highPrecision ? 10.0 : 4.0;

    // [逻辑] 内部计算 LfD
    if(baseParams["L"] > 1e-9) baseParams["LfD"] = baseParams["Lf"] / baseParams["L"];
    else baseParams["LfD"] = 0.0;

    // 生成时间序列
    int nPoints = ui->pointsEdit->text().toInt();
    if(nPoints < 5) nPoints = 5;
    double maxTime = baseParams.value("t", 1000.0);
    if(maxTime < 1e-3) maxTime = 1000.0;
    QVector<double> t = ModelSolver01_06::generateLogTimeSteps(nPoints, -3.0, log10(maxTime));

    int iterations = isSensitivity ? sensitivityValues.size() : 1;
    iterations = qMin(iterations, (int)m_colorList.size());

    QString resultTextHeader = QString("计算完成 (%1)\n").arg(getModelName());
    if(isSensitivity) resultTextHeader += QString("敏感性参数: %1\n").arg(sensitivityKey);

    // 循环计算
    for(int i = 0; i < iterations; ++i) {
        QMap<QString, double> currentParams = baseParams;
        double val = 0;
        if (isSensitivity) {
            val = sensitivityValues[i];
            currentParams[sensitivityKey] = val;

            // 联动更新 LfD
            if (sensitivityKey == "L" || sensitivityKey == "Lf") {
                if(currentParams["L"] > 1e-9) currentParams["LfD"] = currentParams["Lf"] / currentParams["L"];
            }
        }

        ModelCurveData res = calculateTheoreticalCurve(currentParams, t);

        res_tD = std::get<0>(res);
        res_pD = std::get<1>(res);
        res_dpD = std::get<2>(res);

        QColor curveColor = isSensitivity ? m_colorList[i] : Qt::red;
        QString legendName;
        if (isSensitivity) legendName = QString("%1 = %2").arg(sensitivityKey).arg(val);
        else legendName = "理论曲线";

        plotCurve(res, legendName, curveColor, isSensitivity);
    }

    // 更新结果
    QString resultText = resultTextHeader;
    resultText += "t(h)\t\tDp(MPa)\t\tdDp(MPa)\n";
    for(int i=0; i<res_pD.size(); ++i) {
        resultText += QString("%1\t%2\t%3\n").arg(res_tD[i],0,'e',4).arg(res_pD[i],0,'e',4).arg(res_dpD[i],0,'e',4);
    }
    ui->resultTextEdit->setText(resultText);

    ui->chartWidget->getPlot()->rescaleAxes();
    if(plot->xAxis->range().lower <= 0) plot->xAxis->setRangeLower(1e-3);
    if(plot->yAxis->range().lower <= 0) plot->yAxis->setRangeLower(1e-3);
    plot->replot();

    onShowPointsToggled(ui->checkShowPoints->isChecked());
    emit calculationCompleted(getModelName(), baseParams);
}

void WT_ModelWidget::plotCurve(const ModelCurveData& data, const QString& name, QColor color, bool isSensitivity) {
    MouseZoom* plot = ui->chartWidget->getPlot();

    const QVector<double>& t = std::get<0>(data);
    const QVector<double>& p = std::get<1>(data);
    const QVector<double>& d = std::get<2>(data);

    QCPGraph* graphP = plot->addGraph();
    graphP->setData(t, p);
    graphP->setPen(QPen(color, 2, Qt::SolidLine));

    QCPGraph* graphD = plot->addGraph();
    graphD->setData(t, d);

    if (isSensitivity) {
        graphD->setPen(QPen(color, 2, Qt::DashLine));
        graphP->setName(name);
        graphD->removeFromLegend();
    } else {
        graphP->setPen(QPen(Qt::red, 2));
        graphP->setName("压力");
        graphD->setPen(QPen(Qt::blue, 2));
        graphD->setName("压力导数");
    }
}

void WT_ModelWidget::onExportData() {
    if (res_tD.isEmpty()) return;
    QString defaultDir = ModelParameter::instance()->getProjectPath();
    if(defaultDir.isEmpty()) defaultDir = ".";
    QString path = QFileDialog::getSaveFileName(this, "导出CSV数据", defaultDir + "/CalculatedData.csv", "CSV Files (*.csv)");
    if (path.isEmpty()) return;
    QFile f(path);
    if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&f);
        out << "t,Dp,dDp\n";
        for (int i = 0; i < res_tD.size(); ++i) {
            double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
            out << res_tD[i] << "," << res_pD[i] << "," << dp << "\n";
        }
        f.close();
        QMessageBox::information(this, "导出成功", "数据文件已保存");
    }
}

