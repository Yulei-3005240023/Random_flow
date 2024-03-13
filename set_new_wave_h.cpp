#include "set_new_wave_h.h"
#include "ui_set_new_wave_h.h"

set_new_wave_h::set_new_wave_h(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::set_new_wave_h)
{
    ui->setupUi(this);
}

set_new_wave_h::~set_new_wave_h()
{
    delete ui;
}

void set_new_wave_h::set_text(double h1)
{
    QString str = "振幅不建议超过水头/通量平均值的百分之五：" + QString::number(h1) + "m/无量纲" + "\n计算时长为周期的整数倍，为解析解苛求条件。\n根据香农采样定理采样频率必须大于信号频率的两倍。\n所以水位/通量波动信号的周期确定生成必须大于采样周期的两倍，建议取三倍";
    ui->textBrowser->clear();
    ui->textBrowser->append(str);
}

void set_new_wave_h::on_save_clicked()
{
    int type = 0;
    if(ui->radioButton_sin->isChecked()) type = 0;
    else if(ui->radioButton_cos->isChecked()) type = 1;
    double cycle = ui->doubleSpinBox_cycle->value();
    double amplitude = ui->doubleSpinBox_amplitude->value();
    emit wave_info(type, cycle, amplitude); // 发射信号
}
