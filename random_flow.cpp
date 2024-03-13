#include "Random_flow.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <QtCore/qrandom.h>
#include <random>
#include <QRandomGenerator>
#include <QDateTime>
#include <Eigen/Dense>
#include "fftw3.h"
#include <stdlib.h>
#include <time.h>
#include <complex.h>

Random_flow_lib::Random_flow_lib()
{
    ic = 15.0;
    tl = 10;
    st = 1;
    n = 11;
    xl = 10;
    sl = 1;
    m = 11;
}

Random_flow_lib::~Random_flow_lib()
{

}

void Random_flow_lib::step_length(double step_length)
{
    sl = step_length;
}

void Random_flow_lib::step_time(double step_time)
{
    st = step_time;
}

double Random_flow_lib::show_sl()
{
    return sl;
}

double Random_flow_lib::show_st()
{
    return st;
}

void Random_flow_lib::x_length(double length)
{
    xl = length;
}

void Random_flow_lib::t_length(double length)
{
    tl = length;
}

double Random_flow_lib::show_xl()
{
    return xl;
}

double Random_flow_lib::show_tl()
{
    return tl;
}

void Random_flow_lib::set_n_m()
{
    n = std::floor(tl/st) + 1;
    m = std::floor(xl/sl) + 1;
}

int Random_flow_lib::show_n()
{
    return n;
}

int Random_flow_lib::show_m()
{
    return m;
}

void Random_flow_lib::initial_condition(double initial_condition)
{
    ic = initial_condition;
}

double Random_flow_lib::show_ic()
{
    return ic;
}

Random_one_dimension_boussinesq::Random_one_dimension_boussinesq() : Random_flow_lib()
{
    we = 0;
    we_x = 0;
    Sy = 0.1;
    K = 10;
    a = -1; //设置-1为默认值方便后续程序判别
    a_as = -1;
    ha = 1;
    h_l[0] = 1.0;
    h_l[1] = 18.0;
    h_r[0] = 2.0;
    h_r[1] = 0;
    series_number = 1;
    angle = 0;

    use_white_noise_time = false;
    use_white_noise_length = false;
    //seed = 0;
    //rand = new QRandomGenerator(seed);
}

Random_one_dimension_boussinesq::~Random_one_dimension_boussinesq()
{
    //delete rand;
}

void Random_one_dimension_boussinesq::a_()
{
    a = (show_ic() * K) / Sy;
}

void Random_one_dimension_boussinesq::reference_thickness(double ha_)
{
    ha = ha_;
}

void Random_one_dimension_boussinesq::pressure_diffusion_coefficient(double a_)
{
    a = a_;
}

void Random_one_dimension_boussinesq::source_sink_expectation(double we_)
{
    we = we_;
}

void Random_one_dimension_boussinesq::source_sink_expectation_x(double we_)
{
    we_x = we_;
}

void Random_one_dimension_boussinesq::hydraulic_conductivity(double K_)
{
    K = K_;
}

void Random_one_dimension_boussinesq::specific_yield(double Sy_)
{
    Sy = Sy_;
}

void Random_one_dimension_boussinesq::r_boundary(double h, bool Dirichlet, bool Neumann, bool Robin)
{
    if(Dirichlet){
        h_r[0] = 1;
        h_r[1] = h;
    }
    else if(Neumann){
        h_r[0] = 2;
        h_r[1] = h;
    }
    else if(Robin){
        h_r[0] = 3;
        h_r[1] = h;
    }
}

void Random_one_dimension_boussinesq::l_boundary(double h, bool Dirichlet, bool Neumann, bool Robin)
{
    if(Dirichlet){
        h_l[0] = 1;
        h_l[1] = h;
    }
    else if(Neumann){
        h_l[0] = 2;
        h_l[1] = h;
    }
    else if(Robin){
        h_l[0] = 3;
        h_l[1] = h;
    }
}

void Random_one_dimension_boussinesq::random_h_l()
{
    // 创建一个随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd()); // 使用Mersenne Twister引擎
    // 定义三个随机数分布器
    std::uniform_int_distribution<> dis(1, 51); // 生成1到50之间的随机整数
    std::uniform_real_distribution<> dis_1(0, h_l[1]); // 生成0~h_l[1]的随机小数
    std::uniform_int_distribution<> dis_(0, 1); // 生成0或1的随机整数
    // 随机振幅
    double amplitude = dis_1(gen);
    // 随机周期
    double cycle = 0.0;
    while(true){
        cycle = show_tl() / dis(gen); // 依据香农采样定理采样频率必须大于信号频率的两倍
        if (cycle >= 3 * show_st())  // 所以信号周期的随机生成必须大于采样周期的两倍，本程序取三倍
            break;
    }
    set_list_h_l(dis_(gen), amplitude, cycle);
}

void Random_one_dimension_boussinesq::set_list_h_l(double type_function, double amplitude, double cycle)
{
    std::vector<double> Array1 = {type_function, amplitude, cycle};
    list_h_l.push_back(Array1);
}

void Random_one_dimension_boussinesq::del_last_list_h_l() //  删除源汇项列表的最后一项
{
    std::vector<std::vector<double>> list_h_l_new(list_h_l.size() - 1, std::vector<double>(3));
    for (unsigned long long i = 0; i < list_h_l.size() - 1; i++) {
        for (unsigned long long j = 0; j < list_h_l[i].size(); j++) {
            list_h_l_new[i][j] = list_h_l[i][j];
        }
    }
    list_h_l = list_h_l_new;
}

void Random_one_dimension_boussinesq::clear_list_h_l()
{
    std::vector<std::vector<double>> list_h_l_clear;
    list_h_l = list_h_l_clear;
}

void Random_one_dimension_boussinesq::random_h_r()
{
    // 创建一个随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd()); // 使用Mersenne Twister引擎
    // 定义三个随机数分布器
    std::uniform_int_distribution<> dis(1, 51); // 生成1到50之间的随机整数
    std::uniform_real_distribution<> dis_1(0, h_l[1]); // 生成0~h_l[1]的随机小数
    std::uniform_int_distribution<> dis_(0, 1); // 生成0或1的随机整数
    // 随机振幅
    double amplitude = dis_1(gen);
    // 随机周期
    double cycle = 0.0;
    while(true){
        cycle = show_tl() / dis(gen); // 依据香农采样定理采样频率必须大于信号频率的两倍
        if (cycle >= 3 * show_st())  // 所以信号周期的随机生成必须大于采样周期的两倍，本程序取三倍
            break;
    }
    set_list_h_r(dis_(gen), amplitude, cycle);
}

void Random_one_dimension_boussinesq::set_list_h_r(double type_function, double amplitude, double cycle)
{
    std::vector<double> Array1 = {type_function, amplitude, cycle};
    list_h_r.push_back(Array1);
}

void Random_one_dimension_boussinesq::del_last_list_h_r() //  删除源汇项列表的最后一项
{
    std::vector<std::vector<double>> list_h_r_new(list_h_r.size() - 1, std::vector<double>(3));
    for (unsigned long long i = 0; i < list_h_r.size() - 1; i++) {
        for (unsigned long long j = 0; j < list_h_r[i].size(); j++) {
            list_h_r_new[i][j] = list_h_r[i][j];
        }
    }
    list_h_r = list_h_r_new;
}

void Random_one_dimension_boussinesq::clear_list_h_r()
{
    std::vector<std::vector<double>> list_h_r_clear;
    list_h_r = list_h_r_clear;
}

void Random_one_dimension_boussinesq::random_source_sink_term()
{
    // 创建一个随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd()); // 使用Mersenne Twister引擎
    // 定义三个随机数分布器
    std::uniform_int_distribution<> dis(1, 51); // 生成1到50之间的随机整数
    std::uniform_real_distribution<> dis_1(0, we); // 生成0~we的随机小数
    std::uniform_int_distribution<> dis_(0, 1); // 生成0或1的随机整数
    // 随机振幅
    double amplitude = dis_1(gen);
    // 随机周期
    double cycle = 0.0;
    while(true){
        cycle = show_tl() / dis(gen); // 依据香农采样定理采样频率必须大于信号频率的两倍
        if (cycle >= 3 * show_st())  // 所以信号周期的随机生成必须大于采样周期的两倍，本程序取三倍
            break;
    }
    set_list_source_sink_term(dis_(gen), amplitude, cycle);
}

void Random_one_dimension_boussinesq::set_list_source_sink_term(double type_function, double amplitude, double cycle)
{
    std::vector<double> Array1 = {type_function, amplitude, cycle};
    list_source_sink_term.push_back(Array1);
}

void Random_one_dimension_boussinesq::del_last_list_source_sink_term() //  删除源汇项列表的最后一项
{
    std::vector<std::vector<double>> list_source_sink_term_new(list_source_sink_term.size() - 1, std::vector<double>(3));
    for (unsigned long long i = 0; i < list_source_sink_term.size() - 1; i++) {
        for (unsigned long long j = 0; j < list_source_sink_term[i].size(); j++) {
            list_source_sink_term_new[i][j] = list_source_sink_term[i][j];
        }
    }
    list_source_sink_term = list_source_sink_term_new;
}

void Random_one_dimension_boussinesq::clear_list_source_sink_term()
{
    std::vector<std::vector<double>> list_source_sink_term_clear;
    list_source_sink_term = list_source_sink_term_clear;
}

void Random_one_dimension_boussinesq::random_source_sink_term_x()
{
    // 创建一个随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd()); // 使用Mersenne Twister引擎
    // 定义两个随机数分布器
    std::uniform_int_distribution<> dis(1, 51); // 生成1到50之间的随机整数
    std::uniform_real_distribution<> dis_1(0, we_x); // 生成0~we的随机小数
    std::uniform_int_distribution<> dis_(0, 1); // 生成0或1的随机整数
    // 随机振幅
    double amplitude = dis_1(gen);
    // 随机周期
    double cycle = 0.0;
    while(true){
        cycle = show_xl() / dis(gen); // 依据香农采样定理采样频率必须大于信号频率的两倍
        if (cycle >= 3 * show_sl())  // 所以信号周期的随机生成必须大于采样周期的两倍，本程序取三倍
            break;
    }
    set_list_source_sink_term_x(dis_(gen), amplitude, cycle);
}

void Random_one_dimension_boussinesq::set_list_source_sink_term_x(double type_function, double amplitude, double cycle)
{
    std::vector<double> Array1 = {type_function, amplitude, cycle};
    list_source_sink_term_x.push_back(Array1);
}

void Random_one_dimension_boussinesq::del_last_list_source_sink_term_x() //  删除源汇项列表的最后一项
{
    std::vector<std::vector<double>> list_source_sink_term_new(list_source_sink_term_x.size() - 1, std::vector<double>(3));
    for (unsigned long long i = 0; i < list_source_sink_term_x.size() - 1; i++) {
        for (unsigned long long j = 0; j < list_source_sink_term_x[i].size(); j++) {
            list_source_sink_term_new[i][j] = list_source_sink_term_x[i][j];
        }
    }
    list_source_sink_term_x = list_source_sink_term_new;
}

void Random_one_dimension_boussinesq::clear_list_source_sink_term_x()
{
    std::vector<std::vector<double>> list_source_sink_term_clear;
    list_source_sink_term_x = list_source_sink_term_clear;
}

void Random_one_dimension_boussinesq::show()
{
    for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
        for (unsigned long long j = 0; j < list_source_sink_term[i].size(); j++) {
            std::cout << list_source_sink_term[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

double Random_one_dimension_boussinesq:: source_sink_term(double t)
{
    if(use_white_noise_time == true)
    {
        return numberical_value_w[t/show_st()];
    }
    else
    {
        if(list_source_sink_term.size() == 0) return we; // 在不设随机波动以及，给定波动的情况下，输出源汇项期望值
        double w_ = 0; //  w_ = we;
        for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
            if(list_source_sink_term[i][0] == 0){
                w_ += list_source_sink_term[i][1] * std::sin(2 * 3.1514926 * (1 / list_source_sink_term[i][2]) * t);
            }
            else if(list_source_sink_term[i][0] == 1){
                w_ += list_source_sink_term[i][1] * std::cos(2 * 3.1514926 * (1 / list_source_sink_term[i][2]) * t);
            }
        };
        return w_;
    }
}

double Random_one_dimension_boussinesq:: source_sink_term_x(double x)
{
    double w_ = we_x; // 遵循叠加原理"随空间变换的源汇项初始值为0"
    for (unsigned long long i = 0; i < list_source_sink_term_x.size(); i++) {
        if(list_source_sink_term_x[i][0] == 0){
            w_ += list_source_sink_term_x[i][1] * std::sin(2 * 3.1514926 * (1 / list_source_sink_term_x[i][2]) * x);
        }
        else if(list_source_sink_term_x[i][0] == 1){
            w_ += list_source_sink_term_x[i][1] * std::cos(2 * 3.1514926 * (1 / list_source_sink_term_x[i][2]) * x);
        }
    };
    return w_;
}

double Random_one_dimension_boussinesq:: list_h_l_cal(double t)
{
    if(use_white_noise_h_l == true)
    {
        return h_l_value[t/show_st()];
    }
    else
    {
        double h_ = h_l[1];
        for (unsigned long long i = 0; i < list_h_l.size(); i++) {
            if(list_h_l[i][0] == 0){
                h_ += list_h_l[i][1] * std::sin(2 * 3.1514926 * (1 / list_h_l[i][2]) * t);
            }
            else if(list_h_l[i][0] == 1){
                h_ += list_h_l[i][1] * std::cos(2 * 3.1514926 * (1 / list_h_l[i][2]) * t);
            }
        };
        return h_;
    }
}

double Random_one_dimension_boussinesq:: list_h_r_cal(double t)
{
    double h_ = h_r[1];
    for (unsigned long long i = 0; i < list_h_r.size(); i++) {
        if(list_h_r[i][0] == 0){
            h_ += list_h_r[i][1] * std::sin(2 * 3.1514926 * (1 / list_h_r[i][2]) * t);
        }
        else if(list_h_r[i][0] == 1){
            h_ += list_h_r[i][1] * std::cos(2 * 3.1514926 * (1 / list_h_r[i][2]) * t);
        }
    };
    return h_;
}

std::vector<std::vector<double>> Random_one_dimension_boussinesq::share_list_source_sink_term()
{
    return list_source_sink_term;
}

std::vector<std::vector<double>> Random_one_dimension_boussinesq::share_list_source_sink_term_x()
{
    return list_source_sink_term_x;
}

std::vector<std::vector<double>> Random_one_dimension_boussinesq::share_list_h_l()
{
    return list_h_l;
}

std::vector<std::vector<double>> Random_one_dimension_boussinesq::share_list_h_r()
{
    return list_h_r;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve(int how_to_solve)
{
    // 如果未设定压力扩散系数
    if(a==-1){
        a = K / Sy;
    }
    set_n_m(); // 设置X,T轴差分点数
    create_plate(); // 设置底板高程数组
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个差分位置的水头值
    H_ALL.setZero();
    for(int k = 0; k < show_n(); k++) //对行(时间轴)进行扫描
    {
        int iteration_times = 0;  // 迭代运算次数计数
        Eigen::MatrixXd H_previous_iteration(1, show_m()); // 前次迭代的水头厚度(h-z)
        Eigen::MatrixXd H_previous_iteration_1(1, show_m()); // 前次迭代的水头高h
        H_previous_iteration.setZero();
        Eigen::MatrixXd H(show_m(), 1);
        H.setZero();
        // 迭代运算开始
        while (true) {
            Eigen::MatrixXd H_a(show_m(), show_m()); // 每一时刻的系数矩阵
            H_a.setZero();
            int l_a = 0;
            Eigen::MatrixXd H_b(show_m(), 1); // 每一时刻的常数矩阵
            H_b.setZero();
            if (iteration_times == 0 && k != 0){
                H_previous_iteration.row(0) = H_ALL.row(k - 1);  // 前次迭代的当前时刻水头数值,此处未开始计算，使用上一时刻的水头值进行近似
            }
            for(int i  = 0; i < show_m(); i++)  // 对列(X轴)进行扫描
            {
                // 时间边界赋值(初始条件）
                if(k == 0){
                    H_a(l_a, l_a) = 1;
                    H_b(l_a, 0) = show_ic();
                }
                // 左边界赋值
                else if((i - 1) < 0 && h_l[0] == 1){
                    H_a(l_a, l_a) = 1;
                    H_b(l_a, 0) = list_h_l_cal(k * show_st());
                }
                else if((i - 1) < 0 && h_l[0] == 2){
                    // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) / K - source_sink_term_x(i * show_sl()) / K - Sy / (K * show_st()) * H_ALL(k - 1, i) -
                                  2 * show_sl() * list_h_l_cal(k * show_st()) * (H_previous_iteration(0, i) + list_h_l_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl());
                    // 给位置为(i, k)处的水头赋上系数值
                    H_a(l_a, l_a) = -(H_previous_iteration(0, i + 1) + H_previous_iteration(0, i)) / (2 * show_sl() * show_sl()) -
                                    (H_previous_iteration(0,i) + list_h_l_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl()) - Sy / (K * show_st());
                    // 给位置为(i+1, k)处的水头赋上系数值
                    H_a(l_a, l_a + 1) = (H_previous_iteration(0, i + 1) + H_previous_iteration(0, i)) / (
                                            2 * show_sl() * show_sl()) + (H_previous_iteration(0, i) +list_h_l_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl());
                }
                // 右边界赋值
                else if((i + 1) == show_m() && h_r[0] == 1){
                    H_a(l_a, l_a) = 1;
                    H_b(l_a, 0) = list_h_r_cal(k * show_st());
                }
                else if((i + 1) == show_m() && h_r[0] == 2){
                    // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) / K - source_sink_term_x(i * show_sl()) / K - Sy / (K * show_st()) * H_ALL(k - 1, i) +
                                  2 * show_sl() * list_h_r_cal(k * show_st()) * (H_previous_iteration(0, i) + list_h_r_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl());
                    // 给位置为(i, k)处的水头赋上系数值
                    H_a(l_a, l_a) = - (H_previous_iteration(0, i) + list_h_r_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl()) -
                                    (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (2 * show_sl() * show_sl()) - Sy / (K *show_st());
                    // 给位置为(i-1, k)处的水头赋上系数值
                    H_a(l_a, l_a - 1) = (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (
                                            2 * show_sl() * show_sl()) + (H_previous_iteration(0, i) + list_h_r_cal(k * show_st()) * 0.5 * show_sl()) / (show_sl() * show_sl());
                }
                // 非边界部分赋值
                else{
                    // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) / K - source_sink_term_x(i * show_sl()) / K - Sy / (K * show_st()) * H_ALL(k - 1, i);
                        // 给位置为(i, k)处的水头赋上系数值
                    H_a(l_a, l_a) = -(H_previous_iteration(0 ,i + 1) + H_previous_iteration(0, i)) / (2 * show_sl() * show_sl()) -
                                    (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (2 * show_sl() * show_sl()) - Sy / (K * show_st());
                    // 给位置为(i-1，k)处的水头赋上系数值
                    H_a(l_a, l_a - 1) = (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (2 * show_sl() * show_sl());
                    // 给位置为(i+1, k)处的水头赋上系数值
                    H_a(l_a, l_a + 1) = (H_previous_iteration(0, i + 1) + H_previous_iteration(0, i)) / (2 * show_sl() * show_sl());
                }
                l_a += 1;
            }

            if(how_to_solve == 1){ // 使用LU分解法进行求解
                H = H_a.lu().solve(H_b);
            }
            else if (how_to_solve == 2){ // 使用QR分解法进行求解
                H = H_a.householderQr().solve(H_b);
            }
            else if (how_to_solve == 3){ //使用追赶法进行求解
                H = solve_zhuiganfa(H_a, H_b);
            }
            if(k == 0){  // 第零时刻不参与迭代计算
                break;
            }
            // 判断是否满足精度需求
            int precision = 0;
            if(iteration_times == 0){  //第0次计算不参与比较精度，继续计算
                iteration_times += 1;
                H_previous_iteration_1.row(0) = H.col(0);
                for(int o = 0;o < show_m(); o++){
                    H_previous_iteration.row(0)(o) = H.col(0)(o) - plate_elevation(o);
                }
            }
            else{
                for(int u = 0; u < show_m(); u++){
                    if (abs(H_previous_iteration_1(0, u) - H(u, 0)) > 0.000001){
                        precision = 1;
                    }
                }
                if (precision != 1){
                    break;
                }
                else{
                    iteration_times += 1;
                    H_previous_iteration_1.row(0) = H.col(0);
                    for(int o = 0;o < show_m(); o++){
                        H_previous_iteration.row(0)(o) = H.col(0)(o) - plate_elevation(o);
                    }
                }
            }
            if(iteration_times > 10000){
                //qDebug()<<iteration_times;
                break;
            }
        }

        for(int o = 0; o < show_m(); o++)  // 对空间进行扫描，整合成所有适合的计算水头
            H_ALL(k, o) = H(o, 0);

    }
    return H_ALL;
}

std::complex<double> Random_one_dimension_boussinesq::M(double a, double x, double w, double l){
    double imag_r = w / a;
    std::complex<double>r(0,  imag_r);
    std::complex<double>r1 = sqrt(r) ;
    std::complex<double>r2(-r1.real(), -r1.imag());
    std::complex<double> m;
    m = (r1 *exp(r1 * l)* exp(r2 * x) - r2 *exp(r2 * l)* exp(r1 * x)) / (r1 *exp(r1 * l) -  r2 *exp(r2 * l));
    std::complex<double> m_(m.real(), m.imag());
    return m_;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve_an_wt() // 解析解求解计算  只算右边界得了
{
    a_(); // 重新计算压力扩散系数a
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个离散位置的水头值
    H_ALL.setZero();
    std::vector<double> wt(show_n()); // 离散化源汇项储存数组
    double t = 0.0;
    // 把源汇项做离散傅里叶变换
    std::vector<std::complex<double>> dft_wt(show_n());
    std::vector<std::complex<double>> dft_wt_fftw(show_n());
    // 离散化源汇项赋值
//    // 判定是否唯一常数
//    if(list_source_sink_term.size() == 0 && use_white_noise_time == false){
//        // 按频率索引开始赋值
//        for(int k = 0; k < show_n(); k++){
//            if(k == 0){
//                dft_wt[k] = we * show_n();
//            }
//            else{
//                dft_wt[k] = 0.0001;
//            }
//        }
//    }

    for(int i = 0; i < show_n(); i++){
        wt[i] = source_sink_term(t);
        t += show_st();
    }

    // 创建输出数组，用于存储FFT结果
    fftw_complex *out_w = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * show_n());
    // 创建FFT计划
    fftw_plan plan_w = fftw_plan_dft_r2c_1d(show_n(), wt.data(), out_w, FFTW_ESTIMATE);
    // 执行FFT计算
    fftw_execute(plan_w);
    for(int k = 0; k < show_n(); k++){
        dft_wt_fftw[k] = std::complex<double>(out_w[k][0], out_w[k][1]);
    }
    // 释放资源
    fftw_destroy_plan(plan_w);
    fftw_free(out_w);

    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++){
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 * k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = wt[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_wt[k] = Fw;
    }

    std::complex<double>j(0, 1); // 此式子代表虚数j
    //std::complex<double>one(1,0); // 代表实数1
    // 求水头变化量离散傅里叶变换, 就右边界一条线
    std::vector<std::complex<double>> dft_h(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        if(k == 0){ // 初始条件为0
            std::complex<double>Fh0(0, 0);
            dft_h[k] = Fh0;
        }
        else{
            std::complex<double>Fh(0, 0);
            double wSy_real = k / (show_st() *show_n())* Sy;
            std::complex<double>wSy(wSy_real , 0);
            double kw = k / (show_st() * show_n());
            std::complex<double>m = M(a, show_xl(), kw, show_xl());
            Fh = ((j * dft_wt[k]) / wSy) * m - ((j * dft_wt[k]) / wSy);
            dft_h[k] = Fh;
        }
    }

    // 水头变化的傅里叶逆变换
    std::vector<double>h(show_n());

//    std::vector<double>h_fftw(show_n());
//    // 输入和输出数组
//    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * show_n());
//    for(int k=0;k<show_n();k++){
//        in[k][0] = dft_h[k].real();
//        in[k][1] = dft_h[k].imag();
//    }
//    double* out = (double*) fftw_malloc(sizeof(double) * show_n());
//    // 创建逆变换计划
//    fftw_plan plan = fftw_plan_dft_c2r_1d(show_n(), in, out, FFTW_ESTIMATE);
//    // 执行逆变换
//    fftw_execute(plan);
//    for(int n = 0; n < show_n(); n++){
//        h_fftw[n] = out[n] / show_n();
//    }
//    // 释放资源
//    fftw_destroy_plan(plan);
//    fftw_free(in);
//    fftw_free(out);

    for(int n = 0; n < show_n(); n++){
        std::complex<double> h_(0,0);
        for (int k = 0; k < show_n(); k++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (k * n)/ n_;
            std::complex<double>eh(0, imag);
            h_ = h_ + dft_h[k] * exp(eh);
        }
        h[n] = h_.real()* (1/(show_st() *show_n())) / 3.1415926;
    }

    // 误差补偿
    double h_wucha = 0 - h[0];
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        if(k == 0){ // 初始条件为0
            std::complex<double>Fh0(h_wucha * show_n(), 0);
            dft_h[k] = Fh0;
        }
        else{
            std::complex<double>Fh(0, 0);
            double wSy_real = k / (show_st() *show_n())* Sy;
            std::complex<double>wSy(wSy_real , 0);
            double kw = k / (show_st() * show_n());
            std::complex<double>m = M(a, show_xl(), kw, show_xl());
            Fh = ((j * dft_wt[k]) / wSy) * m - ((j * dft_wt[k]) / wSy);
            dft_h[k] = Fh;
        }
    }

    for(int n = 0; n < show_n(); n++){
        std::complex<double> h_(0,0);
        for (int k = 0; k < show_n(); k++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (k * n)/ n_;
            std::complex<double>eh(0, imag);
            h_ = h_ + dft_h[k] * exp(eh);
        }
        h[n] = h_.real()* (1/(show_st() *show_n())) / 3.1415926;
    }

    for(int n = 0; n < show_n(); n++){
        H_ALL(n, (show_m() - 1)) = h[n] + show_ic();
    }
    return H_ALL;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve_an_h_l_t()
{
    a_(); // 重新计算压力扩散系数a
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个离散位置的水头值
    H_ALL.setZero();
    std::vector<double> ht(show_n()); // 离散化左边界储存数组
    double t = 0.0;
    // 离散化左边界赋值
    for(int i = 0; i < show_n(); i++){
        ht[i] = list_h_l_cal(t) - h_l[1];
        t += show_st();
    }
    // 把左边界做离散傅里叶变换
    std::vector<std::complex<double>> dft_hlt(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 * 2 * 3.1415926 * k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = ht[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_hlt[k] = Fw;
    }

//    std::vector<std::complex<double>> dft_wt(show_n());
//        // 创建输出数组，用于存储FFT结果
//    fftw_complex *out_w = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * show_n());
//    // 创建FFT计划
//    fftw_plan plan_w = fftw_plan_dft_r2c_1d(show_n(), ht.data(), out_w, FFTW_ESTIMATE);
//    // 执行FFT计算
//    fftw_execute(plan_w);
//    for(int k = 0; k < show_n(); k++){
//        dft_wt[k] = std::complex<double>(out_w[k][0], out_w[k][1]);
//    }
//    // 释放资源
//    fftw_destroy_plan(plan_w);
//    fftw_free(out_w);

    // 求水头变化量离散傅里叶变换, 就右边界一条线
    std::vector<std::complex<double>> dft_h(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        if(k == 0){ // 初始条件为0
            std::complex<double>Fh0(0, 0);
            dft_h[k] = Fh0;
        }
        else{
            std::complex<double>Fh(0, 0);
            double kw = 6.28 * k * (1 /(show_st() *show_n()));
            std::complex<double>m = M(a, show_xl(), kw, show_xl());
            Fh = dft_hlt[k] * m;
            dft_h[k] = Fh;
        }
    }

//    std::vector<double>h_fftw(show_n());
//    // 输入和输出数组
//    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * show_n());
//    for(int k=0;k<show_n();k++){
//        in[k][0] = dft_h[k].real();
//        in[k][1] = dft_h[k].imag();
//    }
//    double* out = (double*) fftw_malloc(sizeof(double) * show_n());
//    // 创建逆变换计划
//    fftw_plan plan = fftw_plan_dft_c2r_1d(show_n(), in, out, FFTW_ESTIMATE);
//    // 执行逆变换
//    fftw_execute(plan);
//    for(int n = 0; n < show_n(); n++){
//        h_fftw[n] = out[n] / show_n();
//    }
//    // 释放资源
//    fftw_destroy_plan(plan);
//    fftw_free(in);
//    fftw_free(out);

    // 水头变化的傅里叶逆变换
    std::vector<double>h(show_n());
    for(int n = 0; n < show_n(); n++){
        std::complex<double> h_(0,0);
        //double h__ = 0;
        for (int k = 0; k < show_n(); k++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = 2 * 3.1415926 * (k * n)/ n_;
            std::complex<double>eh(0, imag);
            h_ = h_ + dft_h[k] * exp(eh);
        //h__ = dft_h[k].real() * std::cos(2.0 * M_PI * k * n / n_) - dft_h[k].imag() * std::sin(2.0 * M_PI * k * n / n_);
        }
        h[n] = h_.real() * (1/(show_st() *show_n())) * 2;// 3.1415926;
        //h[n] = h__ / show_n();
    }
    for(int n = 0; n < show_n(); n++){
        H_ALL(n, (show_m() - 1)) = h[n] + show_ic();
    }
    return H_ALL;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve_an_wt_h_l_t()
{  // 源汇项随时间变化与左边界随时间变化的解析解求解计算
    a_(); // 重新计算压力扩散系数a
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个离散位置的水头值
    H_ALL.setZero();

    std::vector<double> ht(show_n()); // 离散化左边界储存数组
    double t = 0.0;
    // 离散化左边界赋值
    for(int i = 0; i < show_n(); i++){
        ht[i] = list_h_l_cal(t) - h_l[1];
        t += show_st();
    }
    // 把左边界做离散傅里叶变换
    std::vector<std::complex<double>> dft_hlt(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 * 2 * 3.1415926 * k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = ht[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_hlt[k] = Fw;
    }

    std::vector<double> wt(show_n()); // 离散化源汇项储存数组
    double t_ = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < show_n(); i++){
        wt[i] = source_sink_term(t_);
        t_ += show_st();
    }
    // 把源汇项做离散傅里叶变换
    std::vector<std::complex<double>> dft_wt(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 * 2 * 3.1415926 * k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = wt[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_wt[k] = Fw;
    }

    std::complex<double>j(0, 1); // 此式子代表虚数j
    // 求水头变化量离散傅里叶变换, 就右边界一条线
    std::vector<std::complex<double>> dft_h(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        if(k == 0){ // 初始条件为0
            std::complex<double>Fh0(0, 0);
            dft_h[k] = Fh0;
        }
        else{
            std::complex<double>Fh(0, 0);
            double wSy_real = k * (1 /(show_st() *show_n()))* Sy * 6.2831852;
            std::complex<double>wSy(wSy_real , 0);
            double kw = 6.2831852 * k * (1 /(show_st() *show_n()));
            std::complex<double>m = M(a, show_xl(), kw, show_xl());
            Fh = (((j * dft_wt[k]) / wSy) + dft_hlt[k])* m - ((j * dft_wt[k]) / wSy);
            dft_h[k] = Fh;
        }
    }

    // 水头变化的傅里叶逆变换
    std::vector<double>h(show_n());
    for(int n = 0; n < show_n(); n++){
        std::complex<double> h_(0,0);
        for (int k = 0; k < show_n(); k++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (6.2831852 * k * n)/ n_;
            std::complex<double>eh(0, imag);
            h_ = h_ + dft_h[k] * exp(eh);
        }
        h[n] = h_.real()* (1/(show_st() *show_n())) *2; // 3.1415926;
    }
    for(int n = 0; n < show_n(); n++){
        H_ALL(n, (show_m() - 1)) = h[n] + show_ic();
    }
    return H_ALL;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve_an_wt_h_l_t_half()
{
    // 源汇项随时间变化与左边界随时间变化的解析解求解计算(x = l/2时)
    a_(); // 重新计算压力扩散系数a
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个离散位置的水头值
    H_ALL.setZero();

    std::vector<double> ht(show_n()); // 离散化左边界储存数组
    double t = 0.0;
    // 离散化左边界赋值
    for(int i = 0; i < show_n(); i++){
        ht[i] = list_h_l_cal(t) - h_l[1];
        t += show_st();
    }
    // 把左边界做离散傅里叶变换
    std::vector<std::complex<double>> dft_hlt(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 /* 2 * 3.1415926 */* k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = ht[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_hlt[k] = Fw;
    }

    std::vector<double> wt(show_n()); // 离散化源汇项储存数组
    double t_ = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < show_n(); i++){
        wt[i] = source_sink_term(t_);
        t_ += show_st();
    }
    // 把源汇项做离散傅里叶变换
    std::vector<std::complex<double>> dft_wt(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        std::complex<double>Fw(0, 0);
        for (int n = 0; n < show_n(); n++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (-1 /* 2 * 3.1415926*/ * k * n )/ n_;
            std::complex<double>eF(0, imag);
            std::complex<double>F_ = wt[n] * exp(eF);
            Fw = Fw + F_;
        }
        Fw = Fw * show_st();
        dft_wt[k] = Fw;
    }

    std::complex<double>j(0, 1); // 此式子代表虚数j
    // 求水头变化量离散傅里叶变换, 就右边界一条线
    std::vector<std::complex<double>> dft_h(show_n());
    // 按频率索引开始赋值
    for(int k = 0; k < show_n(); k++){
        if(k == 0){ // 初始条件为0
            std::complex<double>Fh0(0, 0);
            dft_h[k] = Fh0;
        }
        else{
            std::complex<double>Fh(0, 0);
            double wSy_real = k * (1 /(show_st() *show_n()))* Sy;
            std::complex<double>wSy(wSy_real , 0);
            double kw = k * (1 /(show_st() *show_n()));
            std::complex<double>m = M(a, show_xl() / 2, kw, show_xl());
            Fh = (((j * dft_wt[k]) / wSy) + dft_hlt[k])* m - ((j * dft_wt[k]) / wSy);
            dft_h[k] = Fh;
        }
    }

    // 水头变化的傅里叶逆变换
    std::vector<double>h(show_n());
    for(int n = 0; n < show_n(); n++){
        std::complex<double> h_(0,0);
        for (int k = 0; k < show_n(); k++) {
            double n_ = show_n(); // 我也不知道为啥非得先这样。。。。。。
            double imag = (k * n)/ n_;
            std::complex<double>eh(0, imag);
            h_ = h_ + dft_h[k] * exp(eh);
        }
        h[n] = h_.real()* (1/(show_st() *show_n())) / 3.1415926;
    }
    for(int n = 0; n < show_n(); n++){
        H_ALL(n, (show_m() - 1)) = h[n] + show_ic();
    }
    return H_ALL;
}

Eigen::VectorXd Random_one_dimension_boussinesq::fast_fourier_transfrom(Eigen::MatrixXd solution, int n)
{
    //mutex.lock(); // 互斥量保护一下线程，以便于后续多线程计算
    fftw_plan_with_nthreads(1);
    set_n_m();
    // 创建输出数组，用于存储FFT结果
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    // 创建FFT计划
    fftw_plan plan = fftw_plan_dft_r2c_1d(n, solution.data(), out, FFTW_ESTIMATE);
    // 执行FFT计算
    fftw_execute(plan);

    // 归一化和取一半处理
    int N = std::round(n / 2);
    Eigen::VectorXd Amplitude(N); // 用于存储频谱的振幅
    for (int i = 0; i < N; i++) {
        Amplitude(i) = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]); // 取复数的绝对值得到振幅
        if(i == 0){
            Amplitude(i) =  Amplitude(i) / n; // 转换为真实振幅0;
        }
        else{
            Amplitude(i) = Amplitude(i) / (n / 2); // 转换为真实振幅
        }
    }
    // 销毁计划
    fftw_destroy_plan(plan);
    fftw_free(out);
    //mutex.unlock();
    return Amplitude;
}

Eigen::VectorXd Random_one_dimension_boussinesq::power_spectral_density(Eigen::MatrixXd solution, int n)
{
    set_n_m();
    // 创建输出数组，用于存储FFT结果
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    // 创建FFT计划
    fftw_plan plan = fftw_plan_dft_r2c_1d(n, solution.data(), out, FFTW_ESTIMATE);
    // 执行FFT计算
    fftw_execute(plan);

    // 归一化和取一半处理
    int N = std::round(n / 2);
    Eigen::VectorXd Amplitude(N); // 用于功率谱密度的振幅
    for (int i = 0; i < N; i++) {
        Amplitude(i) = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]); // 取复数的绝对值得到振幅
        if(i == 0){
            Amplitude(i) =  (Amplitude(i) / n) * (Amplitude(i) / n) / show_tl(); // 转换为真实振幅0;Sx[k]=(1/N)∣X[k]∣^2
        }
        else{
            Amplitude(i) = (Amplitude(i) / (n / 2)) * (Amplitude(i) / (n / 2)) / show_tl(); // 转换为真实振幅的模的平方
        }
    }
    // 销毁计划
    fftw_destroy_plan(plan);
    fftw_free(out);
    return Amplitude;
}

Eigen::VectorXd Random_one_dimension_boussinesq::amplitude_complete_fdm(Eigen::MatrixXd solution, int l) // l代表位置
{
    Eigen::VectorXd solution_h = solution.col(l);
    Eigen::VectorXd Amplitude_h(show_n());
    for(int i = 0; i < show_n(); i++){
        Amplitude_h[i] = solution_h[i] - show_ic();
    }
    Eigen::VectorXd Amplitude_h_fft = fast_fourier_transfrom(Amplitude_h, show_n());// flow.fast_fourier_transfrom(solve_fdm.col(a), flow.show_n())
    Eigen::VectorXd Amplitude_w(show_n());
    double t = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < show_n(); i++){
        Amplitude_w[i] = source_sink_term(t) - we;
        t += show_st();
    }
    Eigen::VectorXd Amplitude_w_fft = fast_fourier_transfrom(Amplitude_w, show_n());
    Eigen::VectorXd Amplitude_fdm(show_n() / 2);
    //std::cout<<"sbwycwds"<<std::endl;
    for(int j = 0; j < (show_n() / 2); j++){
        Amplitude_fdm[j] = fabs(Amplitude_h_fft[j] / Amplitude_w_fft[j]);
    }
    return Amplitude_fdm;
}

Eigen::VectorXd Random_one_dimension_boussinesq::amplitude_complete_fdm_hl(Eigen::MatrixXd solution, int l) // l代表位置
{
    Eigen::VectorXd Amplitude_h_fft = fast_fourier_transfrom(solution.col(l), show_n());// flow.fast_fourier_transfrom(solve_fdm.col(a), flow.show_n())
    Eigen::VectorXd Amplitude_hl(show_n());
    double t = 0.0;
    // 离散化左边界赋值
    for(int i = 0; i < show_n(); i++){
        Amplitude_hl[i] = list_h_l_cal(t) - h_l[1];
        t += show_st();
    }
    Eigen::VectorXd Amplitude_hl_fft = fast_fourier_transfrom(Amplitude_hl, show_n());
    Eigen::VectorXd Amplitude_fdm(show_n() / 2);
    for(int j = 0; j < (show_n() / 2); j++){
        Amplitude_fdm[j] = fabs(Amplitude_h_fft[j] / Amplitude_hl_fft[j]);
    }
    return Amplitude_fdm;
}

Eigen::VectorXd Random_one_dimension_boussinesq::amplitude_complete_analyze()
{
    a_();
    Eigen::VectorXd Amplitude_analyze(show_n() / 2);
    double e = 0.0;
    for(int i = 0; i < (show_n() / 2); i++){
        e = A(i * (1 / show_tl()) * 6.28);
        Amplitude_analyze[i] = e;
    }
    Amplitude_analyze[0] = Amplitude_analyze[1]; // 把不存在的0的位置给替换掉
    return Amplitude_analyze;
}

Eigen::VectorXd Random_one_dimension_boussinesq::amplitude_complete_analyze_hl()
{
    a_();
    Eigen::VectorXd Amplitude_analyze(show_n() / 2); //取一半有效振幅
    double e = 0.0;
    for(int i = 0; i < (show_n() / 2); i++){
        e = A_hl(i * (1 / show_tl()) * 6.28); // 1/tl = i*st/n=i*st/tl/st 此处角频率记得转化
        Amplitude_analyze[i] = e;
    }
    Amplitude_analyze[0] = 1; // 在频率为0时，理论振幅比应该为1
    return Amplitude_analyze;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve_zhuiganfa(Eigen::MatrixXd a,Eigen::MatrixXd b)
{
    Eigen::MatrixXd x(show_m(), 1); // 解值矩阵
    double *BB,*GG;
    int i;
    //为系数数组分配内存空间
    BB=new double[show_m()];
    GG=new double[show_m()];
    //系数追踪
    BB[0]=a(0,0);
    GG[0]=b(0)/BB[0];
    for(i=1;i<show_m();i++)
    {
        BB[i]=a(i,i)-a(i, i-1)*a(i-1, i)/BB[i-1];
        GG[i]=(b(i)-a(i, i-1)*GG[i-1])/BB[i];
    }
    //求解未知数
    for(i=(show_m()-1);i>=0;i--)
    {
        if(i==(show_m()-1))x(i)=GG[i];
        else  x(i)=GG[i]-a(i, i+1)*x(i+1)/BB[i];
    }
    delete[] BB;
    delete[] GG;
    return x;
}

double Random_one_dimension_boussinesq::A(double w){
    //w = w;

    /*double z_imag = w * show_xl() * show_xl() / a;
    std::complex<double> z(0, z_imag);
    std::complex<double> sqrt_z = std::sqrt(z);
    // 计算双曲余弦的反函数
    std::complex<double> acosh_z = std::acosh(sqrt_z);
    double awl_real = a / (w * show_xl() * show_xl());
    std::complex<double> awl(awl_real, 0);
    double result = (show_xl() * show_xl() / (a * Sy)) * abs(awl * acosh_z);*/


    std::complex<double>j(0, 1); // 此式子代表虚数j
    std::complex<double>one(1,0); // 代表实数1
    std::complex<double>Fh(0, 0);
    double wSy_real = w * Sy;
    std::complex<double>wSy(wSy_real , 0);
    std::complex<double>m = M(a, show_xl(), w, show_xl());
    std::complex<double>Fjwsy = j / wSy;
    Fh = Fjwsy * (m-one);
    double result = fabs(Fh);
    return result;
}

double Random_one_dimension_boussinesq::A_hl(double w){
//    //w = w;
//    double z_imag = w * show_xl() * show_xl() / a;
//    std::complex<double> z(0, z_imag);
//    std::complex<double> sqrt_z = std::sqrt(z);
//    // 计算双曲余弦的反函数
//    std::complex<double> acosh_z = std::acosh(sqrt_z);
//    double awl_real = a / (w * show_xl() * show_xl());
//    std::complex<double> awl(awl_real, 0);
//    double result =  abs(awl * acosh_z);
    double imag_r = w / a;
    std::complex<double>r(0,  imag_r);
    std::complex<double>r1 = sqrt(r) ;
    std::complex<double>r2(-r1.real(), -r1.imag());
    std::complex<double>r_(2, 0);
    std::complex<double>z1 = r_ / (exp(r1 * show_xl()) + exp(r2 * show_xl()));
    double result = abs(z1);
    return result;
}

void Random_one_dimension_boussinesq::set_angle(double a)
{
    angle = a;
}

void Random_one_dimension_boussinesq::create_plate()
{
    plate.clear();
    double m = show_m();
    int sl = show_sl();
    for(int i = 0;i < m;i++)
    {
        plate.emplace_back(tan(angle * 3.1415926 / 180) * i * sl);
    }
}

double Random_one_dimension_boussinesq::plate_elevation(double i)
{
    return plate[i];
}

void Random_one_dimension_boussinesq::set_white_noise_time(int s) // 设置是否启用白噪声模式
{
    if(s==0) use_white_noise_time = false;
    else
    {
        set_n_m();
        numberical_value_w.clear();
        // 创建一个随机数引擎
        double t0 = time(0); // 使用当前系统时间戳为随机数种子
        std::mt19937 gen(t0); // 使用Mersenne Twister引擎
        // 定义两个随机数分布器
        std::uniform_real_distribution<> dis(0, 2 * we); // 生成-1 * we ~ 2 * we均匀分布的随机小数
        use_white_noise_time = true;
        for(int i=0;i<show_n();i++)
        {
            numberical_value_w.emplace_back(dis(gen));// 白噪声数值填充
        }
    }
}

double Random_one_dimension_boussinesq::actual_expectations_white_noise_time()
{
    double expectations = 0;
    for(double e : numberical_value_w)
    {
        expectations += e;
    }
    return expectations / show_n();
}

double Random_one_dimension_boussinesq::fangcha_white_noise_time()
{
    double expectation =0;
    double expectations = 0;
    double d = 0;
    for(double e : numberical_value_w)
    {
        expectation += e;
    }
    expectations = expectation / show_n();
    for(double i : numberical_value_w)
    {
        d += (i - expectations) * (i - expectations);
    }
    return d / show_n();
}

void Random_one_dimension_boussinesq::set_white_noise_h_l(int s) // 设置是否启用白噪声模式
{
    if(s==0) use_white_noise_h_l = false;
    else
    {
        set_n_m();
        h_l_value.clear();
        // 创建一个随机数引擎
        double t0 = time(0); // 使用当前系统时间戳为随机数种子
        std::mt19937 gen(t0); // 使用Mersenne Twister引擎
        // 定义两个随机数分布器
        std::uniform_real_distribution<> dis(0.95 * h_l[1], 1.05 * h_l[1]); // 生成0.95 * h_l[1]~ 1.05 * h_l[1]均匀分布的随机小数
        use_white_noise_h_l = true;
        for(int i=0;i<show_n();i++)
        {
            h_l_value.emplace_back(dis(gen));// 白噪声数值填充
        }
    }
}
