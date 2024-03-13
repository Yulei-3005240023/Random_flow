#ifndef SET_NEW_WAVE_H_H
#define SET_NEW_WAVE_H_H

#include <QDialog>

namespace Ui {
class set_new_wave_h;
}

class set_new_wave_h : public QDialog
{
    Q_OBJECT

public:
    explicit set_new_wave_h(QWidget *parent = nullptr);
    ~set_new_wave_h();
    void set_text(double we);


signals:
    void wave_info(int type, double cycle, double amplitude);

private slots:
    void on_save_clicked();

private:
    Ui::set_new_wave_h *ui;
};

#endif // SET_NEW_WAVE_H_H
