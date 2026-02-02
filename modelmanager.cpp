/*
 * modelmanager.cpp
 * 文件作用: 模型管理类实现文件
 * 修改记录:
 * 1. [性能优化] initializeModels 不再循环创建所有界面，改为调整容器大小。
 * 2. [逻辑修改] switchToModel 和 calculateTheoreticalCurve 中增加 ensureWidget/ensureSolver 检查，实现按需创建。
 */

#include "modelmanager.h"
#include "modelselect.h"
#include "modelparameter.h"
#include "wt_modelwidget.h"
#include "modelsolver01-06.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QDebug>
#include <cmath>

ModelManager::ModelManager(QWidget* parent)
    : QObject(parent), m_mainWidget(nullptr), m_modelStack(nullptr)
    , m_currentModelType(Model_1)
{
}

ModelManager::~ModelManager()
{
    // 清理求解器内存 (Widget 由 Qt 父子对象机制自动清理)
    // 但为了安全，清理非空指针
    for(auto* s : m_solvers) {
        if(s) delete s;
    }
    m_solvers.clear();
}

void ModelManager::initializeModels(QWidget* parentWidget)
{
    if (!parentWidget) return;
    createMainWidget();

    m_modelStack = new QStackedWidget(m_mainWidget);

    // [修改] 不再立即创建对象，而是初始化容器为 6 个空指针
    m_modelWidgets.clear();
    m_modelWidgets.resize(6);
    m_modelWidgets.fill(nullptr);

    m_solvers.clear();
    m_solvers.resize(6);
    m_solvers.fill(nullptr);

    m_mainWidget->layout()->addWidget(m_modelStack);

    // 默认加载第一个模型，此时才会创建 Model_1 的界面，其他 5 个不会创建
    switchToModel(Model_1);

    if (parentWidget->layout()) parentWidget->layout()->addWidget(m_mainWidget);
    else {
        QVBoxLayout* layout = new QVBoxLayout(parentWidget);
        layout->addWidget(m_mainWidget);
        parentWidget->setLayout(layout);
    }
}

void ModelManager::createMainWidget()
{
    m_mainWidget = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_mainWidget);
    mainLayout->setContentsMargins(0, 0, 0, 0); // 无边距
    mainLayout->setSpacing(0);
    m_mainWidget->setLayout(mainLayout);
}

void ModelManager::connectModelSignals()
{
    // 旧的批量连接逻辑已废弃，改在 ensureWidget 中单独连接
}

// [新增] 确保 Widget 存在的辅助函数
WT_ModelWidget* ModelManager::ensureWidget(ModelType type)
{
    int index = (int)type;
    if (index < 0 || index >= m_modelWidgets.size()) return nullptr;

    if (m_modelWidgets[index] == nullptr) {
        // 创建界面对象
        // qDebug() << "惰性加载界面: Model" << index + 1;
        WT_ModelWidget* widget = new WT_ModelWidget(type, m_modelStack);
        m_modelWidgets[index] = widget;
        m_modelStack->addWidget(widget);

        // 连接信号
        connect(widget, &WT_ModelWidget::requestModelSelection, this, &ModelManager::onSelectModelClicked);
        connect(widget, &WT_ModelWidget::calculationCompleted, this, &ModelManager::onWidgetCalculationCompleted);
    }
    return m_modelWidgets[index];
}

// [新增] 确保 Solver 存在的辅助函数
ModelSolver01_06* ModelManager::ensureSolver(ModelType type)
{
    int index = (int)type;
    if (index < 0 || index >= m_solvers.size()) return nullptr;

    if (m_solvers[index] == nullptr) {
        // 创建求解器对象
        // qDebug() << "惰性加载求解器: Model" << index + 1;
        ModelSolver01_06* solver = new ModelSolver01_06(type);
        m_solvers[index] = solver;
    }
    return m_solvers[index];
}

void ModelManager::switchToModel(ModelType modelType)
{
    if (!m_modelStack) return;
    ModelType old = m_currentModelType;
    m_currentModelType = modelType;

    // [关键] 确保目标模型界面已创建
    WT_ModelWidget* w = ensureWidget(modelType);
    if (w) {
        m_modelStack->setCurrentWidget(w);
        // 如果有缓存的观测数据，可能需要传递给新创建的 Widget
        // (WT_ModelWidget 内部通常通过 ModelParameter 获取全局数据，这里视具体实现而定，
        //  若 Widget 需要主动 push 数据，可在此处补充)
    }

    emit modelSwitched(modelType, old);
}

void ModelManager::onSelectModelClicked()
{
    ModelSelect dlg(m_mainWidget);
    if (dlg.exec() == QDialog::Accepted) {
        QString code = dlg.getSelectedModelCode();
        if (code == "modelwidget1") switchToModel(Model_1);
        else if (code == "modelwidget2") switchToModel(Model_2);
        else if (code == "modelwidget3") switchToModel(Model_3);
        else if (code == "modelwidget4") switchToModel(Model_4);
        else if (code == "modelwidget5") switchToModel(Model_5);
        else if (code == "modelwidget6") switchToModel(Model_6);
        else {
            qDebug() << "未知的模型代码: " << code;
        }
    }
}

QString ModelManager::getModelTypeName(ModelType type)
{
    return ModelSolver01_06::getModelName(type);
}

void ModelManager::onWidgetCalculationCompleted(const QString &t, const QMap<QString, double> &r) {
    emit calculationCompleted(t, r);
}

void ModelManager::setHighPrecision(bool high) {
    // 设置界面里的求解器精度（仅对已创建的设置）
    for(WT_ModelWidget* w : m_modelWidgets) {
        if(w) w->setHighPrecision(high);
    }
    // 设置后台求解器精度（仅对已创建的设置）
    for(ModelSolver01_06* s : m_solvers) {
        if(s) s->setHighPrecision(high);
    }
}

void ModelManager::updateAllModelsBasicParameters()
{
    for(WT_ModelWidget* w : m_modelWidgets) {
        if(w) QMetaObject::invokeMethod(w, "onResetParameters");
    }
    // qDebug() << "所有已加载模型的参数已刷新。";
}

QMap<QString, double> ModelManager::getDefaultParameters(ModelType type)
{
    QMap<QString, double> p;
    ModelParameter* mp = ModelParameter::instance();

    p.insert("phi", mp->getPhi());
    p.insert("h", mp->getH());
    p.insert("mu", mp->getMu());
    p.insert("B", mp->getB());
    p.insert("Ct", mp->getCt());
    p.insert("q", mp->getQ());

    p.insert("nf", 4.0);
    p.insert("kf", 1e-3);
    p.insert("km", 1e-4);
    p.insert("L", 1000.0);
    p.insert("Lf", 100.0);
    p.insert("LfD", 0.1);
    p.insert("rmD", 4.0);
    p.insert("omega1", 0.4);
    p.insert("omega2", 0.08);
    p.insert("lambda1", 1e-3);
    p.insert("gamaD", 0.02);

    if (type == Model_1 || type == Model_3 || type == Model_5) {
        p.insert("cD", 0.01);
        p.insert("S", 1.0);
    } else {
        p.insert("cD", 0.0);
        p.insert("S", 0.0);
    }

    if (type == Model_3 || type == Model_4 || type == Model_5 || type == Model_6) {
        p.insert("reD", 10.0);
    }

    return p;
}

ModelCurveData ModelManager::calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    // [修改] 使用 ensureSolver 确保对象存在
    ModelSolver01_06* solver = ensureSolver(type);
    if (solver) {
        return solver->calculateTheoreticalCurve(params, providedTime);
    }
    return ModelCurveData();
}

QVector<double> ModelManager::generateLogTimeSteps(int count, double startExp, double endExp) {
    return ModelSolver01_06::generateLogTimeSteps(count, startExp, endExp);
}

void ModelManager::setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d)
{
    m_cachedObsTime = t;
    m_cachedObsPressure = p;
    m_cachedObsDerivative = d;
}

void ModelManager::getObservedData(QVector<double>& t, QVector<double>& p, QVector<double>& d) const
{
    t = m_cachedObsTime;
    p = m_cachedObsPressure;
    d = m_cachedObsDerivative;
}

void ModelManager::clearCache()
{
    m_cachedObsTime.clear();
    m_cachedObsPressure.clear();
    m_cachedObsDerivative.clear();
}

bool ModelManager::hasObservedData() const
{
    return !m_cachedObsTime.isEmpty();
}
