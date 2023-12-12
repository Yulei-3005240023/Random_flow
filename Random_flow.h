#ifndef RANDOM_FLOW_H
#define RANDOM_FLOW_H

#include<string>
#include<vector>
#include<cmath>
#include<QtCore>
#include<Eigen/Dense>

class Random_flow_lib
{
private:
    std::string name_chinese = "潜水含水层随机流";
    double ic;  // 初始条件
    double tl; // 时间轴轴长 原则上单位为天
    double st; // 时间差分步长
    int n; // 时间轴差分点数目
    double xl; // X轴轴长 原则上单位为米
    double sl; // X轴差分步长
    int m; // X轴差分点数目

    double B=1; // 默认一维流的宽度为1

public:
    Random_flow_lib();
    ~Random_flow_lib();
    void step_length(double step_length);
    void step_time(double step_time);
    double show_sl();
    double show_st();
    void x_length(double length);
    void t_length(double length);
    double show_tl();
    double show_xl();
    void set_n_m();
    int show_n();
    int show_m();
    void initial_condition(double initial_condition);
    double show_ic();
    void width(double width);
};

class Random_one_dimension_boussinesq : public Random_flow_lib
{
private:
    double we; // 源汇项期望
    double we_x;
    double Sy; // 给水度
    double K; // 渗透系数
    double a; // 压力扩散系数
    double a_as; // 解析解压力扩散系数
    double ha; // 参考厚度
    double h_r[2]; // 左右边界
    double h_l[2];
    int series_number; // 表示随机源汇项的傅里叶级数的项数
    std::vector<std::vector<double>> list_source_sink_term; // 制定一个不定长二维数组来实现对源汇项傅里叶级数的储存
    int series_number_x; // 表示随机源汇项(x轴)的傅里叶级数的项数
    std::vector<std::vector<double>> list_source_sink_term_x; // 制定一个不定长二维数组来实现对源汇项(x轴)傅里叶级数的储存
    std::string name_chinese = "潜水含水层随机非稳定一维流";
    double angle; // 底板倾斜角度
    std::vector<double> plate; // 底板高程数组
    std::vector<double> numberical_value_w; // 源汇项数值组

protected:
    bool use_white_noise_time;
    bool use_white_noise_length;
    qint64 seed; // 随机数种子
    QRandomGenerator *rand;


public:
    Random_one_dimension_boussinesq();
    ~Random_one_dimension_boussinesq();

    void a_(); // 计算压力扩散系数
    void reference_thickness(double ha_);
    void pressure_diffusion_coefficient(double a_);
    void source_sink_expectation(double we_);
    void source_sink_expectation_x(double we_);
    void hydraulic_conductivity(double K_);
    void specific_yield(double Sy_);
    void r_boundary(double h, bool Dirichlet=false, bool Neumann=false, bool Robin=false);
    void l_boundary(double h, bool Dirichlet=false, bool Neumann=false, bool Robin=false);

    void random_source_sink_term();
    void set_list_source_sink_term(double type_function, double amplitude, double cycle);
    void del_last_list_source_sink_term();
    void clear_list_source_sink_term();

    void random_source_sink_term_x();
    void set_list_source_sink_term_x(double type_function, double amplitude, double cycle);
    void del_last_list_source_sink_term_x();
    void clear_list_source_sink_term_x();
    void show();

    double source_sink_term(double x);
    double source_sink_term_x(double x);

    std::vector<std::vector<double>> share_list_source_sink_term();
    std::vector<std::vector<double>> share_list_source_sink_term_x();
    Eigen::MatrixXd solve(int how_to_solve);
    Eigen::VectorXd fast_fourier_transfrom(Eigen::MatrixXd solution, int n); // 快速傅里叶变换代码
    Eigen::VectorXd power_spectral_density(Eigen::MatrixXd solution, int n); // 功率谱密度变换代码
    Eigen::VectorXd amplitude_complete_fdm(Eigen::MatrixXd solution, int l); // 功率谱振幅比值代码-数值解
    Eigen::VectorXd amplitude_complete_analyze(); // 功率谱振幅比-解析解

    Eigen::MatrixXd solve_zhuiganfa(Eigen::MatrixXd a,Eigen::MatrixXd b); // 追赶法算法
    std::complex<double> M(double x, double w, double l); // 解析解中的M
    double A(double w); // 解析解，返回给定频率的功率谱振幅比值

    void set_angle(double a); // 设置底板倾斜角度
    void create_plate(); // 创建底板高程数组
    double plate_elevation(double i); // 底板高程数据输出

    void set_white_noise_time(int s);
    double actual_expectations_white_noise_time(); // 白噪声实际期望值计算



};
#endif // RANDOM_FLOW_H
