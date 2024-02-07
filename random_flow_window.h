#ifndef RANDOM_FLOW_WINDOW_H
#define RANDOM_FLOW_WINDOW_H

#include "Random_flow.h"
#include <QMainWindow>
#include <QtCharts>
#include <mcthread.h>

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
    void clear_chart_head(); // 对于绘制水头图标进行清除上一次绘制曲线的操作
    void save_the_data(QString filename, int model);

public slots:
    void get_wave_info(double cycle, double amplitue); // 主窗口获得波动信息的槽函数

    void get_MC_times(int time); // 主窗口获得蒙特卡洛模拟进行次数的槽函数
    void get_MC_amplitude_complete_fdm(Eigen::VectorXd vector1); // 主窗口获得蒙特卡洛模拟结果的槽函数
    void get_MC_finished(); //获得进程结束的槽函数

    //2号线程对应的槽函数
    void get_MC_times_1(int time); // 主窗口获得蒙特卡洛模拟进行次数的槽函数
    void get_MC_amplitude_complete_fdm_1(Eigen::VectorXd vector1); // 主窗口获得蒙特卡洛模拟结果的槽函数
    void get_MC_finished_1(); //获得进程结束的槽函数

private slots:
    void do_mouseMovePoint(QPoint point); // 鼠标移动监测

    void on_actionSet_fdm_triggered();

    void on_solve_FDM_clicked();

    void on_doubleSpinBox_rain_valueChanged(double arg1);

    void on_random_new_wave_clicked();

    void on_draw_solve_line_clicked();

    void on_draw_solve_line_location_clicked();

    void on_draw_solve_line_as_clicked();

    void on_time_field_figure_clicked();

    void on_spinBox_x_length_valueChanged(int arg1);

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

    void on_amplitude_complete_figure_clicked();

    void on_actionsave_fdm_triggered(); // 储存数值解

    void on_MC_un_wt_amp_start_clicked(); // 关于源汇项随时间成均匀分布的蒙特卡洛线程启动函数

    void on_spinBox_MC_un_wt_amp_valueChanged(int arg1);

    void on_checkBox_use_MC_clicked();

    void on_doubleSpinBox_left_boundary_valueChanged(double arg1);

    void on_actionsave_MC_triggered();

    void on_actionsave_as_triggered();

    void on_actionsave_as_complete_triggered();

    void on_doubleSpinBox_right_boundary_valueChanged(double arg1);

    void on_random_new_wave_h_l_clicked();

    void on_delete_wave_h_l_clicked();

    void on_random_new_wave_h_r_clicked();


    void on_delete_wave_h_r_clicked();

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
    QLineSeries *series_analyze;
    QValueAxis *axis_head;
    QValueAxis *axis_x;
    QLogValueAxis *logAxisX; // 对数横坐标轴
    void create_chart_head();
    QChart *chart_W;
    QLineSeries *series_W;
    QValueAxis *axis_W;
    QValueAxis *axis_w;
    void create_chart_W();
    QLabel *lab_chartXY;
    bool move_the_chart;
    bool select_the_chart;

    Eigen::VectorXd MC_amplitude_complete_fdm; // 蒙特卡洛数值模拟的功率谱结果矩阵存放(总控)
    Eigen::VectorXd MC_amplitude_complete_fdm_1; // 蒙特卡洛数值模拟的功率谱结果矩阵存放(线程1)
    Eigen::VectorXd MC_amplitude_complete_fdm_2; // 蒙特卡洛数值模拟的功率谱结果矩阵存放(线程2)
    int MC_times; // 设定的蒙特卡洛模拟次数
    int MC_act_time = 0; //实际上蒙特卡洛运行的次数
    int MC_act_time_1 = 0; //实际上蒙特卡洛运行的次数(线程2)
    MCThread_uniform_wt_amp *MCThread_uniform_wt_amp1; // 源汇项随时间均匀分布的振幅的蒙特卡洛模拟线程1
    MCThread_uniform_wt_amp *MCThread_uniform_wt_amp2; // 源汇项随时间均匀分布的振幅的蒙特卡洛模拟线程2
    bool MCThread_uniform_wt_amp1_work = false; // 线程1是否工作
    bool MCThread_uniform_wt_amp2_work = false; // 线程2是否工作
    bool use_MC_to_draw = false; // 是否使用蒙特卡洛的数据去绘图
    QMutex mutex; // 互斥量，用于确保线程之间不会冲突
};
#endif // RANDOM_FLOW_WINDOW_H
