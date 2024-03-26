#ifndef MCTHREAD_H
#define MCTHREAD_H

#include <QThread>
#include <Random_flow.h>
//#include <QDeclarativeItem>
#include <QtQuick/QQuickItem>

class MCThread : public QThread
{
    Q_OBJECT
public:
    explicit MCThread(QObject *parent = nullptr);
};

class MCThread_uniform_wt_amp : public MCThread // 是用于对源汇项随时间均匀分布的振幅的蒙特卡洛模拟线程
{
    Q_OBJECT
private:
    //Eigen::VectorXd MC_amplitude_complete_fdm; // 蒙特卡洛数值模拟的功率谱结果矩阵存放
    Random_one_dimension_boussinesq flow; // 实例化具象潜水一维随机流
    int times = 10; // 蒙特卡洛运算次数
    QMutex mutex; // 互斥量
    bool use_preheat = false; // 是否启用预热期模式
protected:
    void run(); // 线程的任务函数
public:
    explicit MCThread_uniform_wt_amp(QObject *parent = nullptr);
    void set_flow(Random_one_dimension_boussinesq flow1);
    void set_times(int t); // 设置运算次数
    void set_preheat(int u); // 设置是否使用预热期
signals:
    void which_time(int time);
    void MC_amp_fdm(Eigen::VectorXd Matrix1);
};

class MCThread_uniform_hl_amp : public MCThread // 是用于对左边界随时间均匀分布的振幅的蒙特卡洛模拟线程
{
    Q_OBJECT
private:
    //Eigen::VectorXd MC_amplitude_complete_fdm; // 蒙特卡洛数值模拟的功率谱结果矩阵存放
    Random_one_dimension_boussinesq flow; // 实例化具象潜水一维随机流
    int times = 10; // 蒙特卡洛运算次数
    QMutex mutex; // 互斥量
protected:
    void run(); // 线程的任务函数
public:
    explicit MCThread_uniform_hl_amp(QObject *parent = nullptr);
    void set_flow(Random_one_dimension_boussinesq flow1);
    void set_times(int t); // 设置运算次数
signals:
    void which_time(int time);
    void MC_amp_fdm(Eigen::VectorXd Matrix1);
};
#endif // MCTHREAD_H
