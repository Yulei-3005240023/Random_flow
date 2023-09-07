#include "set_fdm.h"
#include "ui_set_fdm.h"
#include <iostream>

Set_FDM::Set_FDM(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Set_FDM)
{
    ui->setupUi(this);
    sl = 1;
    st = 1;
    solve_method = 1;
}

Set_FDM::~Set_FDM()
{
    delete ui;
}

void Set_FDM::set_window(double sl_, double st_, int solve_method_)
{
    sl = sl_;
    ui->spinBox_step_length->setValue(sl);
    st = st_;
    ui->spinBox_step_time->setValue(st);
    solve_method = solve_method_;
    if(solve_method == 1) ui->LU->click();
    else if(solve_method == 2) ui->QR->click();
    else if(solve_method == 3) ui->Cholesky->click();
}

int Set_FDM::step_length()
{
    return ui->spinBox_step_length->value();
}

int Set_FDM::step_time()
{
    return ui->spinBox_step_time->value();
}

int Set_FDM::solve()
{
    return solve_method;
}

void Set_FDM::how_to_solve()
{
    if(ui->LU->isChecked()) solve_method = 1;
    else if(ui->QR->isChecked()) solve_method = 2;
    else if(ui->Cholesky->isChecked()) solve_method = 3;
}

