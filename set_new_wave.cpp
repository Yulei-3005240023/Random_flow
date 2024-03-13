#include "set_new_wave.h"
#include "ui_set_new_wave.h"

set_new_wave::set_new_wave(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::set_new_wave)
{
    ui->setupUi(this);
}

set_new_wave::~set_new_wave()
{
    delete ui;
}

void set_new_wave::set_text(double we)
{
    we = we / 365; // 年平均值转化为天
    QString str = "振幅不得超过一天的降雨量平均值：" + QString::number(we) + "mm\n计算时长为周期的整数倍，为解析解苛求条件。\n根据香农采样定理采样频率必须大于信号频率的两倍。\n所以降雨波动信号的周期确定生成必须大于采样周期的两倍，建议取三倍";
    ui->textBrowser->clear();
    ui->textBrowser->append(str);
}

void set_new_wave::on_save_clicked()
{
    int type = 0;
    if(ui->radioButton_sin->isChecked()) type = 0;
    else if(ui->radioButton_cos->isChecked()) type = 1;
    double cycle = ui->doubleSpinBox_cycle->value();
    double amplitude = ui->doubleSpinBox_amplitude->value() / 1000;
    emit wave_info(type, cycle, amplitude); // 发射信号
}
