#ifndef RANDOM_FLOW_WINDOW_H
#define RANDOM_FLOW_WINDOW_H

#include "Random_flow.h"
#include <QMainWindow>
#include <QtCharts>

QT_BEGIN_NAMESPACE
namespace Ui { class Random_flow_Window; }
QT_END_NAMESPACE

class Random_flow_Window : public QMainWindow
{
    Q_OBJECT

public:
    Random_flow_Window(QWidget *parent = nullptr);
    ~Random_flow_Window();
    void left_boundary();
    void right_boundary();
    void set_time_choose_box();
    void set_location_choose_box();
    void solve(); // 数值解计算函数

public slots:
    void get_wave_info(double cycle, double amplitue);

private slots:
    void do_mouseMovePoint(QPoint point);  // 鼠标移动监测
    void on_actionSet_fdm_triggered();

    void on_solve_FDM_clicked();

    void on_doubleSpinBox_rain_valueChanged(double arg1);

    void on_random_new_wave_clicked();

    void on_draw_solve_line_clicked();

    void on_draw_solve_line_location_clicked();

    void on_time_field_figure_clicked();

    void on_spinBox_t_length_valueChanged(int arg1);

    void on_frequency_field_figure_clicked();

    void on_draw_solve_line_location_fft_clicked();

    void on_actionSet_hydrogeological_parameter_triggered();

    void on_delete_wave_clicked();

    void on_new_wave_clicked();

    void on_doubleSpinBox_rain_x_valueChanged(double arg1);

    void on_random_new_wave_x_clicked();

    void on_delete_wave_x_clicked();

    void on_draw_solve_line_fft_clicked();

    void on_time_field_figure_x_clicked();

    void on_frequency_field_figure_x_clicked();

    void on_move_clicked();

    void on_frame_destroyed();

    void on_graphicsView_destroyed();

    void on_select_clicked();

    void on_reset_1_clicked();

    void on_spinBox_angle_valueChanged(double arg1);

    void on_doubleSpinBox_angle_valueChanged(double arg1);

    void on_use_white_noise_checkBox_time_clicked();

    void on_actual_expectations_clicked();

    void on_power_spectral_density_figure_clicked();

private:
    Ui::Random_flow_Window *ui;
    Random_one_dimension_boussinesq flow; // 实例化具象潜水一维随机流
    double step_length;
    double step_time;
    int how_to_solve;
    double hydraulic_conductivity;
    double specific_yield;
    Eigen::MatrixXd solve_fdm; // 数值解求解后的矩阵存放
    QChart *chart_head;
    QLineSeries *series_head;
    QLineSeries *series_plate;
    QValueAxis *axis_head;
    QValueAxis *axis_x;
    void create_chart_head();
    QChart *chart_W;
    QLineSeries *series_W;
    QValueAxis *axis_W;
    QValueAxis *axis_w;
    void create_chart_W();
    QLabel *lab_chartXY;
    bool move_the_chart;
    bool select_the_chart;
};
#endif // RANDOM_FLOW_WINDOW_H
