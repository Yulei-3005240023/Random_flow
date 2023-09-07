#include "set_hydrogeological_parameter.h"
#include "ui_set_hydrogeological_parameter.h"

Set_hydrogeological_parameter::Set_hydrogeological_parameter(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Set_hydrogeological_parameter)
{
    ui->setupUi(this);
}

Set_hydrogeological_parameter::~Set_hydrogeological_parameter()
{
    delete ui;
}

void Set_hydrogeological_parameter::set_window(double K_, double Sy_)
{
    K = K_;
    Sy = Sy_;
    ui->doubleSpinBox_hydraulic_conductivity->setValue(K);
    ui->doubleSpinBox_specific_yield->setValue(Sy);
}

double Set_hydrogeological_parameter::hydraulic_conductivity()
{
    return ui->doubleSpinBox_hydraulic_conductivity->value();
}

double Set_hydrogeological_parameter::specific_yield()
{
    return ui->doubleSpinBox_specific_yield->value();
}
