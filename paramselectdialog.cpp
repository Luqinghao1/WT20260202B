/*
 * 文件名: paramselectdialog.cpp
 * 文件作用: 参数选择配置对话框的具体实现文件
 * 功能描述:
 * 1. 初始化参数配置表格，显示参数名称、当前值、单位、拟合状态、上下限及调节步长。
 * 2. 提供了对参数上下限和步长的编辑功能。
 * 3. [修复] 修复了 setStripTrailingZeros 编译错误。
 * 4. [优化] 通过自定义 SmartDoubleSpinBox 类，实现了"数值显示只保留有效数字"的需求。
 * 5. 实现了滚轮事件过滤，防止在表格滚动时误触数值输入框导致参数改变。
 */

#include "paramselectdialog.h"
#include "ui_paramselectdialog.h"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QDebug>

// --- [新增] 自定义 SpinBox 以去除末尾多余的 0 ---
class SmartDoubleSpinBox : public QDoubleSpinBox {
public:
    explicit SmartDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent) {}

    // 重写文本格式化函数
    QString textFromValue(double value) const override {
        // 使用 'g' 格式，自动去除末尾多余的 0，保留有效数字
        // prec 参数使用 decimals() 设置的精度
        return QString::number(value, 'g', decimals());
    }
};
// ------------------------------------------------

// 构造函数
ParamSelectDialog::ParamSelectDialog(const QList<FitParameter> &params, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParamSelectDialog),
    m_params(params)
{
    ui->setupUi(this);
    this->setWindowTitle("拟合参数配置");

    // 连接底部按钮信号
    connect(ui->btnOk, &QPushButton::clicked, this, &ParamSelectDialog::onConfirm);
    connect(ui->btnCancel, &QPushButton::clicked, this, &ParamSelectDialog::onCancel);

    // 防止回车键意外触发取消
    ui->btnCancel->setAutoDefault(false);

    // 初始化表格内容
    initTable();
}

// 析构函数
ParamSelectDialog::~ParamSelectDialog()
{
    delete ui;
}

// 事件过滤器
// 作用：拦截 QDoubleSpinBox 的鼠标滚轮事件。
// 防止用户在滚动表格试图浏览时，鼠标经过输入框导致数值被滚轮修改。
bool ParamSelectDialog::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Wheel) {
        if (qobject_cast<QAbstractSpinBox*>(obj)) {
            // 如果是滚轮事件且目标是输入框，则忽略该事件（返回 true 表示已处理）
            return true;
        }
    }
    return QDialog::eventFilter(obj, event);
}

// 初始化表格
// 作用：根据传入的 m_params 列表构建表格行和列
void ParamSelectDialog::initTable()
{
    QStringList headers;
    // 定义列头：显示顺序调整为符合用户习惯
    // 0:显示, 1:当前数值, 2:单位, 3:参数名称, 4:拟合变量(是否拟合), 5:下限, 6:上限, 7:滚轮步长
    headers << "显示" << "当前数值" << "单位" << "参数名称" << "拟合变量" << "下限" << "上限" << "滚轮步长";
    ui->tableWidget->setColumnCount(headers.size());
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->setRowCount(m_params.size());

    // 复选框样式定义
    QString checkBoxStyle =
        "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
        "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
        "QCheckBox::indicator:hover { border-color: #0078d7; }";

    for(int i = 0; i < m_params.size(); ++i) {
        const FitParameter& p = m_params[i];

        // --- Col 0: 显示 (Visible) ---
        QWidget* pWidgetVis = new QWidget();
        QHBoxLayout* pLayoutVis = new QHBoxLayout(pWidgetVis);
        QCheckBox* chkVis = new QCheckBox();
        chkVis->setChecked(p.isVisible);
        chkVis->setStyleSheet(checkBoxStyle);
        pLayoutVis->addWidget(chkVis);
        pLayoutVis->setAlignment(Qt::AlignCenter);
        pLayoutVis->setContentsMargins(0,0,0,0);
        ui->tableWidget->setCellWidget(i, 0, pWidgetVis);

        // --- Col 1: 当前数值 (Value) ---
        // 使用自定义的 SmartDoubleSpinBox
        SmartDoubleSpinBox* spinVal = new SmartDoubleSpinBox();
        spinVal->setRange(-9e9, 9e9); // 极大范围
        spinVal->setDecimals(10);     // 设置允许的最大小数位数 (显示时会自动精简)
        // spinVal->setStripTrailingZeros(true); // [已删除] 该函数不存在
        spinVal->setValue(p.value);
        spinVal->setFrame(false);     // 无边框，融入表格
        spinVal->installEventFilter(this); // 安装过滤器防止误触滚轮
        ui->tableWidget->setCellWidget(i, 1, spinVal);

        // --- Col 2: 单位 (Unit) ---
        QString dummy, dummy2, dummy3, unitStr;
        // 获取参数的静态显示信息
        FittingParameterChart::getParamDisplayInfo(p.name, dummy, dummy2, dummy3, unitStr);
        if(unitStr == "无因次" || unitStr == "小数") unitStr = "-";
        QTableWidgetItem* unitItem = new QTableWidgetItem(unitStr);
        unitItem->setFlags(unitItem->flags() & ~Qt::ItemIsEditable); // 只读
        ui->tableWidget->setItem(i, 2, unitItem);

        // --- Col 3: 参数名称 (Name) ---
        QString displayNameFull = QString("%1 (%2)").arg(p.displayName).arg(p.name);
        QTableWidgetItem* nameItem = new QTableWidgetItem(displayNameFull);
        nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable); // 只读
        nameItem->setData(Qt::UserRole, p.name); // 存储原始键名用于后续识别
        ui->tableWidget->setItem(i, 3, nameItem);

        // --- Col 4: 拟合变量 (Fit) ---
        QWidget* pWidgetFit = new QWidget();
        QHBoxLayout* pLayoutFit = new QHBoxLayout(pWidgetFit);
        QCheckBox* chkFit = new QCheckBox();
        chkFit->setChecked(p.isFit);
        chkFit->setStyleSheet(checkBoxStyle);
        pLayoutFit->addWidget(chkFit);
        pLayoutFit->setAlignment(Qt::AlignCenter);
        pLayoutFit->setContentsMargins(0,0,0,0);

        // 特殊处理：LfD (无因次缝长) 是计算结果，不允许用户选择拟合
        if (p.name == "LfD") {
            chkFit->setEnabled(false);
            chkFit->setChecked(false);
        }
        ui->tableWidget->setCellWidget(i, 4, pWidgetFit);

        // 联动逻辑：如果勾选“拟合”，则该参数必须“显示”
        connect(chkFit, &QCheckBox::checkStateChanged, [chkVis](Qt::CheckState state){
            if (state == Qt::Checked) {
                chkVis->setChecked(true);
                chkVis->setEnabled(false); // 强制显示，不可取消
                // 变灰样式提示不可用
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                      "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
            } else {
                chkVis->setEnabled(true); // 恢复可编辑
                chkVis->setStyleSheet(
                    "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
                    "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
                    "QCheckBox::indicator:hover { border-color: #0078d7; }"
                    );
            }
        });

        // 初始化时的联动状态设置
        if (p.isFit) {
            chkVis->setChecked(true);
            chkVis->setEnabled(false);
            chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                  "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
        }

        // --- Col 5: 下限 (Min) ---
        SmartDoubleSpinBox* spinMin = new SmartDoubleSpinBox();
        spinMin->setRange(-9e9, 9e9);
        spinMin->setDecimals(10);
        // spinMin->setStripTrailingZeros(true); // [已删除]
        spinMin->setValue(p.min);
        spinMin->setFrame(false);
        spinMin->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 5, spinMin);

        // --- Col 6: 上限 (Max) ---
        SmartDoubleSpinBox* spinMax = new SmartDoubleSpinBox();
        spinMax->setRange(-9e9, 9e9);
        spinMax->setDecimals(10);
        // spinMax->setStripTrailingZeros(true); // [已删除]
        spinMax->setValue(p.max);
        spinMax->setFrame(false);
        spinMax->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 6, spinMax);

        // --- Col 7: 滚轮步长 (Step) ---
        SmartDoubleSpinBox* spinStep = new SmartDoubleSpinBox();
        spinStep->setRange(0.0, 10000.0);
        spinStep->setDecimals(10);
        // spinStep->setStripTrailingZeros(true); // [已删除]
        spinStep->setValue(p.step);
        spinStep->setFrame(false);
        spinStep->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 7, spinStep);
    }

    ui->tableWidget->resizeColumnsToContents();
    // 让“参数名称”列自适应拉伸，填满剩余空间
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
}

// 收集数据
// 作用：遍历表格控件，将用户修改后的值回写到 m_params 列表中
void ParamSelectDialog::collectData()
{
    for(int i = 0; i < ui->tableWidget->rowCount(); ++i) {
        if(i >= m_params.size()) break;

        // 0: Visible
        QWidget* wVis = ui->tableWidget->cellWidget(i, 0);
        if (wVis) {
            QCheckBox* chkVis = wVis->findChild<QCheckBox*>();
            if(chkVis) m_params[i].isVisible = chkVis->isChecked();
        }

        // 1: Value
        QDoubleSpinBox* spinVal = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 1));
        if(spinVal) m_params[i].value = spinVal->value();

        // 4: Fit
        QWidget* wFit = ui->tableWidget->cellWidget(i, 4);
        if (wFit) {
            QCheckBox* chkFit = wFit->findChild<QCheckBox*>();
            if(chkFit) m_params[i].isFit = chkFit->isChecked();
        }

        // 5: Min
        QDoubleSpinBox* spinMin = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 5));
        if(spinMin) m_params[i].min = spinMin->value();

        // 6: Max
        QDoubleSpinBox* spinMax = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 6));
        if(spinMax) m_params[i].max = spinMax->value();

        // 7: Step
        QDoubleSpinBox* spinStep = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 7));
        if(spinStep) m_params[i].step = spinStep->value();
    }
}

// 获取更新后的参数列表
QList<FitParameter> ParamSelectDialog::getUpdatedParams() const
{
    return m_params;
}

// 确认按钮槽函数
void ParamSelectDialog::onConfirm()
{
    collectData(); // 收集界面数据
    accept();      // 关闭对话框并返回 Accepted
}

// 取消按钮槽函数
void ParamSelectDialog::onCancel()
{
    reject();      // 关闭对话框并返回 Rejected
}
