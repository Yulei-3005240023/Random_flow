#ifndef SET_HYDROGEOLOGICAL_PARAMETER_H
#define SET_HYDROGEOLOGICAL_PARAMETER_H

#include <QDialog>

namespace Ui {
class Set_hydrogeological_parameter;
}

class Set_hydrogeological_parameter : public QDialog
{
    Q_OBJECT

public:
    explicit Set_hydrogeological_parameter(QWidget *parent = nullptr);
    ~Set_hydrogeological_parameter();
    void set_window(double K_, double Sy_);
    double hydraulic_conductivity();
    double specific_yield();

private:
    Ui::Set_hydrogeological_parameter *ui;
    double K;
    double Sy;


};

#endif // SET_HYDROGEOLOGICAL_PARAMETER_H
