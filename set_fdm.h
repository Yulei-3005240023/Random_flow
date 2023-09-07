#ifndef SET_FDM_H
#define SET_FDM_H

#include <QDialog>

namespace Ui {
class Set_FDM;
}

class Set_FDM : public QDialog
{
    Q_OBJECT

public:
    explicit Set_FDM(QWidget *parent = nullptr);
    ~Set_FDM();
    int step_length();
    int step_time();
    int solve();
    void how_to_solve();
    void set_window(double sl_, double st_, int solve_method_);


private:
    Ui::Set_FDM *ui;
    int sl;
    int st;
    int solve_method;
};

#endif // SET_FDM_H
