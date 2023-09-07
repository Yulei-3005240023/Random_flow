#include "Random_flow.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "fftw3.h"

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

double Random_flow_lib::show_ic(){
    return ic;
}

Random_one_dimension_boussinesq::Random_one_dimension_boussinesq() : Random_flow_lib()
{
    we = 0.05;
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
}

Random_one_dimension_boussinesq::~Random_one_dimension_boussinesq()
{

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

void Random_one_dimension_boussinesq::random_source_sink_term()
{
    // 创建一个随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd()); // 使用Mersenne Twister引擎
    // 定义一个随机数分布器
    std::uniform_int_distribution<> dis(1, 51); // 生成1到50之间的随机整数
    // 随机振幅
    double amplitude = we / dis(gen);
    // 随机周期
    double cycle = 0.0;
    while(true){
        cycle = show_tl() / dis(gen); // 依据香农采样定理采样频率必须大于信号频率的两倍
        if (cycle >= 3 * show_st())  // 所以信号周期的随机生成必须大于采样周期的两倍，本程序取三倍
            break;
    }
    set_list_source_sink_term(0, amplitude, cycle);
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
    double w_ = we;
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

std::vector<std::vector<double>> Random_one_dimension_boussinesq::share_list_source_sink_term()
{
    return list_source_sink_term;
}

Eigen::MatrixXd Random_one_dimension_boussinesq::solve(int how_to_solve)
{
    // 如果未设定压力扩散系数
    if(a==-1){
        a = K / Sy;
    }
    set_n_m(); // 设置X,T轴差分点数
    Eigen::MatrixXd H_ALL(show_n(), show_m());//创建一个矩阵，用于存放各个差分位置的水头值
    H_ALL.setZero();
    for(int k = 0; k < show_n(); k++) //对行(时间轴)进行扫描
    {
        int iteration_times = 0;  // 迭代运算次数计数
        Eigen::MatrixXd H_previous_iteration(1, show_m());
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
                    H_b(l_a, 0) = h_l[1];
                }
                else if((i - 1) < 0 && h_l[0] == 2){
                    // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) - Sy / (K * show_st()) * H_ALL(k - 1, i) -
                                  2 * show_sl() * h_l[1] * (H_previous_iteration(0, i) + h_l[1] * 0.5 * show_sl()) / (show_sl() * show_sl());
                    // 给位置为(i, k)处的水头赋上系数值
                    H_a(l_a, l_a) = -(H_previous_iteration(0, i + 1) + H_previous_iteration(0, i)) / (2 * show_sl() * show_sl()) -
                                    (H_previous_iteration(0,i) + h_l[1] * 0.5 * show_sl()) / (show_sl() * show_sl()) - Sy / (K * show_st());
                    // 给位置为(i+1, k)处的水头赋上系数值
                    H_a(l_a, l_a + 1) = (H_previous_iteration(0, i + 1) + H_previous_iteration(0, i)) / (
                                            2 * show_sl() * show_sl()) + (H_previous_iteration(0, i) + h_l[1] * 0.5 * show_sl()) / (show_sl() * show_sl());
                }
                // 右边界赋值
                else if((i + 1) == show_m() && h_r[0] == 1){
                    H_a(l_a, l_a) = 1;
                    H_b(l_a, 0) = h_r[1];
                }
                else if((i + 1) == show_m() && h_r[0] == 2){
                    // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) - Sy / (K * show_st()) * H_ALL(k - 1, i) +
                                  2 * show_sl() * h_r[1] * (H_previous_iteration(0, i) + h_r[1] * 0.5 * show_sl()) / (show_sl() * show_sl());
                    // 给位置为(i, k)处的水头赋上系数值
                    H_a(l_a, l_a) = - (H_previous_iteration(0, i) + h_r[1] * 0.5 * show_sl()) / (show_sl() * show_sl()) -
                                    (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (2 * show_sl() * show_sl()) - Sy / (K *show_st());
                    // 给位置为(i-1, k)处的水头赋上系数值
                    H_a(l_a, l_a - 1) = (H_previous_iteration(0, i) + H_previous_iteration(0, i - 1)) / (
                                            2 * show_sl() * show_sl()) + (H_previous_iteration(0, i) + h_r[1] * 0.5 * show_sl()) / (show_sl() * show_sl());
                }
                // 非边界部分赋值
                else{
                        // 源汇项赋值
                    H_b(l_a, 0) = - source_sink_term(k * show_st()) - Sy / (K * show_st()) * H_ALL(k - 1, i);
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
            else if (how_to_solve == 3){ //使用Cholesky分解法进行求解
                H = H_a.ldlt().solve(H_b);
            }
            if(k == 0){  // 第零时刻不参与迭代计算
                break;
            }
            // 判断是否满足精度需求
            int precision = 0;
            for(int u = 0; u < show_m(); u++){
                if (abs(H_previous_iteration(0, u) - H(u, 0)) > 0.001){
                    precision = 1;
                }
            }
            if (precision != 1){
                break;
            }
            else{
                iteration_times += 1;
                H_previous_iteration.row(0) = H.col(0);
            }
            if(iteration_times > 100)
            break;
        }
        //std::cout<<iteration_times<<std::endl;
        for(int o = 0; o < show_m(); o++)  // 对空间进行扫描，整合成所有适合的计算水头
            H_ALL(k, o) = H(o, 0);

    }
    return H_ALL;
}

Eigen::VectorXd Random_one_dimension_boussinesq::fast_fourier_transfrom(Eigen::MatrixXd solution, int n)
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
    Eigen::VectorXd Amplitude(N); // 用于存储频谱的振幅
    for (int i = 0; i < N; i++) {
        Amplitude(i) = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]); // 取复数的绝对值得到振幅
        if(i == 0){
            Amplitude(i) =  0;
        }
        else{
            Amplitude(i) = Amplitude(i) / (n / 2); // 转换为真实振幅
        }
    }
    // 销毁计划
    fftw_destroy_plan(plan);
    fftw_free(out);
    return Amplitude;
}
//int main() {
//    Random_one_dimension_boussinesq flow;
//    flow.set_list_source_sink_term(1,0.06,5);
//    Eigen::MatrixXd h = flow.solve(1);
//    std::cout<<h;

//    return 0;
//}
