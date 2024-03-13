#include "set_new_wave_hh.h"
#include "ui_set_new_wave_hh.h"

set_new_wave_hh::set_new_wave_hh(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::set_new_wave_hh)
{
    ui->setupUi(this);
}

set_new_wave_hh::~set_new_wave_hh()
{
    delete ui;
}
