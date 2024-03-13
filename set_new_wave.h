#ifndef SET_NEW_WAVE_H
#define SET_NEW_WAVE_H

#include <QDialog>

namespace Ui {
class set_new_wave;
}

class set_new_wave : public QDialog
{
    Q_OBJECT

public:
    explicit set_new_wave(QWidget *parent = nullptr);
    ~set_new_wave();
    void set_text(double we);


signals:
    void wave_info(int type, double cycle, double amplitude);

private slots:
    void on_save_clicked();

private:
    Ui::set_new_wave *ui;
};

#endif // SET_NEW_WAVE_H
