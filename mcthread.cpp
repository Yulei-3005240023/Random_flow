#include "mcthread.h"

MCThread::MCThread(QObject *parent)
    : QThread{parent}
{

}

MCThread_uniform_wt_amp::MCThread_uniform_wt_amp(QObject *parent)
    :  MCThread{parent}
{

}

void MCThread_uniform_wt_amp::run()
{
    int time = 0;
    if(use_preheat){
        while(time<times)
        {
            flow.set_white_noise_time(1);
            time += 1;
            mutex.lock();
            Eigen::MatrixXd solve_fdm = flow.solve(3);
            Eigen::VectorXd Amplitude_fdm = flow.amplitude_complete_fdm_preheat(solve_fdm, loc);
            mutex.unlock();
            emit which_time(time); // 向主进程发送信号，表示进行到第几次模拟
            //msleep(10);
            emit MC_amp_fdm(Amplitude_fdm);
        }
    }
    else{
        while(time<times)
        {
            flow.set_white_noise_time(1);
            time += 1;
            mutex.lock();
            Eigen::MatrixXd solve_fdm = flow.solve(3);
            Eigen::VectorXd Amplitude_fdm = flow.amplitude_complete_fdm(solve_fdm, loc);
            mutex.unlock();
            emit which_time(time); // 向主进程发送信号，表示进行到第几次模拟
            //msleep(10);
            emit MC_amp_fdm(Amplitude_fdm);
        }
    }
    quit();
}

void MCThread_uniform_wt_amp::set_flow(Random_one_dimension_boussinesq flow1)
{
    flow = flow1;
}

void MCThread_uniform_wt_amp::set_times(int t)
{
    times = t;
}

void MCThread_uniform_wt_amp::set_preheat(int u)
{
    if(u == 1){
        use_preheat = true;
    }
    else{
        use_preheat = false;
    }
}

void MCThread_uniform_wt_amp::set_loc(int l)
{
    loc = l;
}

MCThread_uniform_hl_amp::MCThread_uniform_hl_amp(QObject *parent)
    :  MCThread{parent}
{

}

void MCThread_uniform_hl_amp::run()
{
    int time = 0;
    while(time<times)
    {
        flow.set_white_noise_h_l(1);
        time += 1;
        mutex.lock();
        Eigen::MatrixXd solve_fdm = flow.solve(3);
        Eigen::VectorXd Amplitude_fdm = flow.amplitude_complete_fdm_hl(solve_fdm, loc);
        mutex.unlock();
        emit which_time(time); // 向主进程发送信号，表示进行到第几次模拟
        //msleep(10);
        emit MC_amp_fdm(Amplitude_fdm);
    }
    quit();
}

void MCThread_uniform_hl_amp::set_flow(Random_one_dimension_boussinesq flow1)
{
    flow = flow1;
}

void MCThread_uniform_hl_amp::set_times(int t)
{
    times = t;
}

void MCThread_uniform_hl_amp::set_loc(int l)
{
    loc = l;
}
