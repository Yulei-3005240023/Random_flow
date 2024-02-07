#include "random_flow_window.h"
#include "ui_random_flow_window.h"
#include "set_fdm.h"
#include "set_hydrogeological_parameter.h"
#include "set_new_wave.h"
#include "Random_flow.h"
#include <iostream>
#include <chrono>


Random_flow_Window::Random_flow_Window(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::Random_flow_Window)
{
    ui->setupUi(this);
    step_length = 1;
    step_time = 1;
    how_to_solve = 1;
    hydraulic_conductivity = 10;
    specific_yield = 0.1;
    chart_head = new QChart();
    series_head = new QLineSeries();
    series_plate = new QLineSeries();
    series_analyze = new QLineSeries();
    axis_head = new QValueAxis();
    axis_x = new QValueAxis();
    logAxisX = new QLogValueAxis();
    create_chart_head();
    chart_W = new QChart();
    series_W = new QLineSeries();
    axis_W = new QValueAxis();
    axis_w = new QValueAxis();

    // 判断图表交互操作的状态
    move_the_chart = false;
    select_the_chart = false;

    create_chart_W();
    // 记录当前系统时间
    auto now = std::chrono::system_clock::now();
    // 输出当前时间
    std::time_t time_now = std::chrono::system_clock::to_time_t(now);
    QString str = std::ctime(&time_now);
    ui->textBrowser->append("当前程程序时间：" + str);
    lab_chartXY = new QLabel("Chart X=, Y=", this);  // 状态栏显示坐标
    ui->statusbar->addWidget(lab_chartXY);
    connect(ui->graphicsView, SIGNAL(mouseMovePoint(QPoint)), this, SLOT(do_mouseMovePoint(QPoint)));

    // 蒙特卡洛多线程部分定义
    //MC_amplitude_complete_fdm = new Eigen::VectorXd;
    MCThread_uniform_wt_amp1 = new MCThread_uniform_wt_amp(this);
    connect(MCThread_uniform_wt_amp1, &MCThread_uniform_wt_amp::which_time, this, &Random_flow_Window::get_MC_times);
    connect(MCThread_uniform_wt_amp1, &MCThread_uniform_wt_amp::MC_amp_fdm, this, &Random_flow_Window::get_MC_amplitude_complete_fdm);
    connect(MCThread_uniform_wt_amp1, &MCThread_uniform_wt_amp::finished, this, &Random_flow_Window::get_MC_finished);

    MCThread_uniform_wt_amp2 = new MCThread_uniform_wt_amp(this);
    connect(MCThread_uniform_wt_amp2, &MCThread_uniform_wt_amp::which_time, this, &Random_flow_Window::get_MC_times_1);
    connect(MCThread_uniform_wt_amp2, &MCThread_uniform_wt_amp::MC_amp_fdm, this, &Random_flow_Window::get_MC_amplitude_complete_fdm_1);
    connect(MCThread_uniform_wt_amp2, &MCThread_uniform_wt_amp::finished, this, &Random_flow_Window::get_MC_finished_1);
}

Random_flow_Window::~Random_flow_Window()
{
    delete ui;
    delete chart_head;
    delete series_head;
    delete series_plate;
    delete series_analyze;
    delete axis_head;
    delete axis_x;
    delete chart_W;
    delete series_W;
    delete axis_W;
    delete axis_w;
    delete logAxisX;
    //delete MC_amplitude_complete_fdm;
    std::cout<<"game is over"<<std::endl;
}

void Random_flow_Window::create_chart_head() // 创建初始图表
{
    ui->graphicsView->setChart(chart_head);
    ui->graphicsView->setRenderHint(QPainter::Antialiasing);  // 抗锯齿
    chart_head->setTitle(tr("水头绘制曲线:示例"));
    QMargins W;
    W.setBottom(0);
    W.setTop(0);
    W.setLeft(0);
    W.setRight(0);
    chart_head->setMargins(W);
    axis_x->setRange(0, 10);
    axis_x->setLabelFormat("%.1f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("X轴(m)");

    axis_head->setLabelFormat("%.1f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("水头(m)");

    series_head->setName("水头曲线: 本图表为示例");
    series_head->append(0, 1);
    series_head->append(2, 7);
    series_head->append(6, 2);
    series_head->append(8, 9);

    chart_head->addSeries(series_head);
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);

}

void Random_flow_Window::create_chart_W() // 创建初始图表
{
    ui->graphicsView_W->setChart(chart_W);
    ui->graphicsView_W->setRenderHint(QPainter::Antialiasing);  // 抗锯齿
    chart_W->setTitle(tr("源汇项绘制曲线:示例"));
    QMargins W;
    W.setBottom(0);
    W.setTop(0);
    W.setLeft(0);
    W.setRight(0);
    chart_W->setMargins(W);
    axis_w->setRange(0, 10);
    axis_w->setLabelFormat("%.1f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("T轴(d)");

    axis_W->setLabelFormat("%.1f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("水头(m)");

    for(int i = 0; i < 11; i++){
        series_W->append(i, std::sin(i));
    }

    chart_W->addSeries(series_W);
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::left_boundary()
{
    QString str = QString::number(ui->doubleSpinBox_left_boundary->value());

    if(ui->comboBox_left_boundary->currentText() == "一类边界(给定水头)"){
        flow.l_boundary(ui->doubleSpinBox_left_boundary->value(), true, false, false);
        ui->textBrowser->append("左边界（给定水头边界）：\n水头值为：" + str);
    }
    else if(ui->comboBox_left_boundary->currentText() == "二类边界(给定通量)"){
        flow.l_boundary(ui->doubleSpinBox_left_boundary->value(), false, true, false);
        ui->textBrowser->append("左边界（给定通量边界）：\n通量值为：" + str);
    }
}

void Random_flow_Window::right_boundary()
{
    QString str = QString::number(ui->doubleSpinBox_right_boundary->value());

    if(ui->comboBox_right_boundary->currentText() == "一类边界(给定水头)"){
        flow.r_boundary(ui->doubleSpinBox_right_boundary->value(), true, false, false);
        ui->textBrowser->append("右边界（给定水头边界）：\n水头值为：" + str);
    }
    else if(ui->comboBox_right_boundary->currentText() == "二类边界(给定通量)"){
        flow.r_boundary(ui->doubleSpinBox_right_boundary->value(), false, true, false);
        ui->textBrowser->append("右边界（给定通量边界）：\n通量值为：" + str);
    }
}

void Random_flow_Window::set_time_choose_box() // 设置时刻选择条
{
    QString str_tl = QString::number(flow.show_tl());
    QString str_st = QString::number(flow.show_st());
    QString str_n = QString::number(flow.show_n());
    ui->spinBox_time->setMaximum(flow.show_n() - 1);
    ui->textBrowser_time->setPlainText("计算时长为：" + str_tl + "天，时间分割步长为：" + str_st + "天，共计有时刻" + str_n + "个。");
}

void Random_flow_Window::set_location_choose_box()  // 设置位置选择条
{
    QString str_xl = QString::number(flow.show_xl());
    QString str_sl = QString::number(flow.show_sl());
    QString str_m = QString::number(flow.show_m());
    ui->spinBox_X->setMaximum(flow.show_m() - 1);
    ui->textBrowser_X->setPlainText("计算轴长为：" + str_xl + "米，空间分割步长为：" + str_sl + "米，共计有位置离散点" + str_m + "个。");
}

void Random_flow_Window::solve()
{
    flow.initial_condition(ui->doubleSpinBox_ic->value());
    left_boundary();
    right_boundary();
    flow.x_length(ui->spinBox_x_length->value());
    flow.t_length(ui->spinBox_t_length->value());
    solve_fdm =  flow.solve(how_to_solve);
    ui->textBrowser->append("计算完毕");
    set_time_choose_box();
    set_location_choose_box();
}

void Random_flow_Window::clear_chart_head()
{
    chart_head->removeSeries(series_head); // 清除原有图表
    chart_head->removeSeries(series_plate);
    chart_head->removeAxis(axis_x);
    chart_head->removeAxis(axis_head);
    chart_head->removeAxis(logAxisX);
    chart_head->removeSeries(series_analyze);
    series_analyze->clear();
    series_head->clear();
    series_plate->clear();
}

void Random_flow_Window::save_the_data(QString filename, int model){
    QFile afile(filename);
    if(!afile.open(QIODevice::WriteOnly | QIODevice::Text)){ // 以只写模式打开
        qDebug() << "无法打开文件";
    }
    QTextStream astream(&afile);
    if(model == 0){ // 写入数值解数据
        for(double h:solve_fdm.col(flow.show_m() - 1)){
            QString str = QString::number(h, 'f', 8);
            astream<<str<<"\n";
        };
    }
    else if (model == 1){ // 写入蒙特卡洛数据
        for(double h:MC_amplitude_complete_fdm){
            QString str = QString::number(h, 'f', 8);
            astream<<str<<"\n";
        }
    }
    else if (model == 2){ // 写入解析解数据
        Eigen::MatrixXd solve_as1 = flow.solve_an_wt();
        for(double h:solve_as1.col(flow.show_m() - 1)){
            QString str = QString::number(h, 'f', 8);
            astream<<str<<"\n";
        }
    }
    else if(model == 3){ // 写入解析解数据(振幅比)
        Eigen::VectorXd amplitude_complete_analyze = flow.amplitude_complete_analyze();
        for(double h:amplitude_complete_analyze){
            QString str = QString::number(h, 'f', 8);
            astream<<str<<"\n";
        }
    }
    afile.close();
}

void Random_flow_Window::get_wave_info(double cycle, double amplitue) // 主窗口获得波动信息的槽函数
{
    flow.set_list_source_sink_term(0, amplitue, cycle);
    std::vector<std::vector<double>> list_source_sink_term = flow.share_list_source_sink_term();
    QString str = QString::number(ui->doubleSpinBox_rain->value()/ (1000 * 365));
    for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
        if(list_source_sink_term[i][0] == 0){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "sin((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
        }
        else if(list_source_sink_term[i][0] == 1){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "cos((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
        }
    };
    ui->textBrowser_rain_function->clear();
    ui->textBrowser_rain_function->append(str);
}

void Random_flow_Window::get_MC_times(int time) // 主窗口获得蒙特卡洛模拟进行次数的槽函数
{
    mutex.lock();
    MC_act_time = time ;
    int t = floor((MC_act_time + MC_act_time_1) * 100 / MC_times) ;
    ui->progressBar_MC_uwt->setValue(t);
    mutex.unlock();
}

void Random_flow_Window::get_MC_amplitude_complete_fdm(Eigen::VectorXd vector1) // 主窗口获得蒙特卡洛模拟结果的槽函数
{
    mutex.lock();
    if(MC_act_time == 1) MC_amplitude_complete_fdm_1 = vector1;
    else
    {
        for(int i = 0; i < (flow.show_n() / 2);i++)
        {
            MC_amplitude_complete_fdm_1[i] = (MC_amplitude_complete_fdm_1[i] * (MC_act_time - 1) + vector1[i]) / MC_act_time; // 取加权平均数，即为所以结果总共的平均数

        }
    }
    mutex.unlock();
}

void Random_flow_Window::get_MC_finished() //获得进程结束的槽函数
{
    mutex.lock();
    MCThread_uniform_wt_amp1_work = false;
    ui->textBrowser->append("线程1已结束");
    if(MCThread_uniform_wt_amp2_work == false)
    {
        ui->progressBar_MC_uwt->setValue(100);
        ui->MC_un_wt_amp_start->setEnabled(true);
        ui->textBrowser->append("所有线程均已结束");
        for(int i = 0; i < (flow.show_n() / 2);i++)
        {
            MC_amplitude_complete_fdm[i] = (MC_amplitude_complete_fdm_1[i] + MC_amplitude_complete_fdm_2[i]) / 2; // 两个线程的结果取平均数以得到最后的结果
        }
    }
    else
    {
        MC_amplitude_complete_fdm = MC_amplitude_complete_fdm_1;
    }
    mutex.unlock();
}

void Random_flow_Window::get_MC_times_1(int time) // 主窗口获得蒙特卡洛模拟进行次数的槽函数
{
    mutex.lock();
    MC_act_time_1 = time ;
    int t = floor((MC_act_time + MC_act_time_1) * 100 / MC_times) ;
    ui->progressBar_MC_uwt->setValue(t);
    mutex.unlock();
}

void Random_flow_Window::get_MC_amplitude_complete_fdm_1(Eigen::VectorXd vector1) // 主窗口获得蒙特卡洛模拟结果的槽函数
{
    mutex.lock();
    if(MC_act_time_1 == 1) MC_amplitude_complete_fdm_2 = vector1;
    else
    {
        for(int i = 0; i < (flow.show_n() / 2);i++)
        {
            MC_amplitude_complete_fdm_2[i] = (MC_amplitude_complete_fdm_2[i] * (MC_act_time_1 - 1) + vector1[i]) / MC_act_time_1; // 取加权平均数，即为所以结果总共的平均数
        }
    }
    mutex.unlock();
}

void Random_flow_Window::get_MC_finished_1() //获得进程结束的槽函数
{
    mutex.lock();
    MCThread_uniform_wt_amp2_work = false;
    ui->textBrowser->append("线程2已结束");
    if(MCThread_uniform_wt_amp1_work == false)
    {
        ui->progressBar_MC_uwt->setValue(100);
        ui->MC_un_wt_amp_start->setEnabled(true);
        ui->textBrowser->append("所有线程均已结束");
        for(int i = 0; i < (flow.show_n() / 2);i++)
        {
            MC_amplitude_complete_fdm[i] = (MC_amplitude_complete_fdm_1[i] + MC_amplitude_complete_fdm_2[i]) / 2; // 两个线程的结果取平均数以得到最后的结果
        }
    }
    else
    {
        MC_amplitude_complete_fdm = MC_amplitude_complete_fdm_2;
    }
    mutex.unlock();
}

void Random_flow_Window::do_mouseMovePoint(QPoint point)
{
    QPointF pt = chart_head->mapToValue(point);
    QString str = QString::asprintf("Chart X=%.2f,Y=%.2f", pt.x(), pt.y());
    lab_chartXY->setText(str);
}

void Random_flow_Window::on_actionSet_fdm_triggered() // 此函数用于打开设置数值解求解设置的子窗口并设置数值
{
    Set_FDM *Set_FDM_window = new Set_FDM(this);
    Set_FDM_window->set_window(step_length, step_time, how_to_solve);
    int ret = Set_FDM_window->exec();
    if(ret == QDialog::Accepted){
        step_length = Set_FDM_window->step_length(); // 获取差分步长
        step_time =  Set_FDM_window->step_time(); // 获取差分时长
        Set_FDM_window->how_to_solve(); // 判断求解方式
        how_to_solve = Set_FDM_window->solve(); // 获取求解方式
        //std::cout<<step_length<<std::endl<<step_time<<std::endl<<how_to_solve<<std::endl;
        flow.step_length(step_length);
        flow.step_time(step_time);
    }
    delete Set_FDM_window;
}

void Random_flow_Window::on_actionSet_hydrogeological_parameter_triggered() // 此函数用于打开设置水文地质学参数的子窗口并设置数值
{
    Set_hydrogeological_parameter *Set_hydrogeological_parameter_window = new Set_hydrogeological_parameter(this);
    Set_hydrogeological_parameter_window->set_window(hydraulic_conductivity, specific_yield);
    int ret = Set_hydrogeological_parameter_window->exec();
    if(ret == QDialog::Accepted){
        hydraulic_conductivity = Set_hydrogeological_parameter_window->hydraulic_conductivity(); // 获取渗透系数
        specific_yield=  Set_hydrogeological_parameter_window->specific_yield(); // 获取给水度
        flow.hydraulic_conductivity(hydraulic_conductivity);
        flow.specific_yield(specific_yield);
    }
    delete Set_hydrogeological_parameter_window;
}

void Random_flow_Window::on_solve_FDM_clicked()
{
    // 开始时间
    auto start = std::chrono::high_resolution_clock::now();

    solve();

    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();

    // 计算时间差
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // 输出执行时间
    ui->textBrowser->append("用时" + QString::number(duration) + "毫秒\n");

}

void Random_flow_Window::on_doubleSpinBox_rain_valueChanged(double arg1)
{
    QString str = QString::number(arg1/ (1000 * 365));
    flow.source_sink_expectation(arg1 / (1000 * 365));  // 一年的降雨量期望转化为一天的
    ui->textBrowser_rain_function->clear();  // 重置展示降雨量函数的框
    flow.clear_list_source_sink_term();  // 重置降雨量波动列表
    ui->textBrowser_rain_function->append(str);
}

void Random_flow_Window::on_doubleSpinBox_rain_x_valueChanged(double arg1)
{
    QString str = QString::number(arg1/ (1000 * 365));
    flow.source_sink_expectation_x(arg1 / (1000 * 365));  // 一年的降雨量期望转化为一天的
    ui->textBrowser_rain_function_x->clear();  // 重置展示降雨量函数的框
    flow.clear_list_source_sink_term_x();  // 重置降雨量波动列表
    ui->textBrowser_rain_function_x->append(str);
}

void Random_flow_Window::on_spinBox_x_length_valueChanged(int arg1)
{
    flow.x_length(arg1);
}

void Random_flow_Window::on_spinBox_t_length_valueChanged(int arg1)
{
    flow.t_length(arg1);
}

void Random_flow_Window::on_random_new_wave_clicked()
{
    if(ui->spinBox_t_length->value() != 0) flow.t_length(ui->spinBox_t_length->value());
    flow.random_source_sink_term();
    std::vector<std::vector<double>> list_source_sink_term = flow.share_list_source_sink_term();
    QString str = QString::number(ui->doubleSpinBox_rain->value()/ (1000 * 365));
    for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
        if(list_source_sink_term[i][0] == 0){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "sin((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
        }
        else if(list_source_sink_term[i][0] == 1){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "cos((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
        }
    };
    ui->textBrowser_rain_function->clear();
    ui->textBrowser_rain_function->append(str);
}

void Random_flow_Window::on_new_wave_clicked() // 此函数用于手动新建一个波动的子窗口并设置数值
{
    if(ui->spinBox_t_length->value() != 0) flow.t_length(ui->spinBox_t_length->value());
    set_new_wave *set_new_wave_window = new set_new_wave;
    set_new_wave_window->setAttribute(Qt::WA_DeleteOnClose); // 对话框关闭时自动删除
    set_new_wave_window->set_text(ui->doubleSpinBox_rain->value());
    connect(set_new_wave_window, &set_new_wave::wave_info, this, &Random_flow_Window::get_wave_info); // 绑定信号和槽
    set_new_wave_window->setModal(false);
    set_new_wave_window->show();
}

void Random_flow_Window::on_delete_wave_clicked()
{
    if(flow.share_list_source_sink_term().size() == 0) return;
    else{
        flow.del_last_list_source_sink_term();
        std::vector<std::vector<double>> list_source_sink_term = flow.share_list_source_sink_term();
        QString str = QString::number(ui->doubleSpinBox_rain->value()/ (1000 * 365));
        for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
            if(list_source_sink_term[i][0] == 0){
                str += " + " + QString::number(list_source_sink_term[i][1]) + "sin((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
            }
            else if(list_source_sink_term[i][0] == 1){
                str += " + " + QString::number(list_source_sink_term[i][1]) + "cos((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * t)";
            }
        };
        ui->textBrowser_rain_function->clear();
        ui->textBrowser_rain_function->append(str);
    }
}

void Random_flow_Window::on_draw_solve_line_clicked()
{
    clear_chart_head();

    QString title = "数值解，空间差分步长为" +QString::number(flow.show_sl()) + "时间差分步长为" + QString::number(flow.show_st()) + "，绘图时刻为第" + QString::number(ui->spinBox_time->value()) + "时刻。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    series_head->setName("水头曲线: 数值解");
    series_plate->setName("底板");
    int a = ui->spinBox_time->value();
    double x = 0.0;
    double min_h = 100.0;
    double max_h = 0.0;
    for (double h : solve_fdm.row(a)) {
        series_head->append(x, h);
        x += flow.show_sl();
        if(h < min_h) min_h = h;
        if(h > max_h) max_h = h;
    }
    double x_1 = 0.0;
    for (int z = 0;z < flow.show_m();z++)
    {
        series_plate->append(x_1,flow.plate_elevation(z));
        x_1 += flow.show_sl();
    }

    axis_x->setRange(0, flow.show_xl());
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("X轴(m)");

    axis_head->setRange(0, max_h);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("水头(m)");

    chart_head->addSeries(series_head); // 更新图表
    chart_head->addSeries(series_plate);
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);
    series_plate->attachAxis(axis_x);
    series_plate->attachAxis(axis_head);

}

void Random_flow_Window::on_draw_solve_line_location_clicked()
{
    clear_chart_head();

    QString title = "数值解，空间差分步长为" +QString::number(flow.show_sl()) + "时间差分步长为" + QString::number(flow.show_st()) + "，绘图位置为第" + QString::number(ui->spinBox_X->value()) + "位置。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    series_head->setName("水头曲线: 数值解");
    int a = ui->spinBox_X->value();
    double x = 0.0;
    double min_h = 800.0;
    double max_h = 0.0;
    for (double h : solve_fdm.col(a)) {
        series_head->append(x, h);
        x += flow.show_st();
        if(h < min_h) min_h = h;
        if(h > max_h) max_h = h;
    }

    axis_x->setRange(0, flow.show_tl());
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("时间轴(d)");

    axis_head->setRange(min_h, max_h);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("水头(m)");

    chart_head->addSeries(series_head); // 更新图表
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);
}

void Random_flow_Window::on_time_field_figure_clicked()
{
    chart_W->removeSeries(series_W); // 清除原有图表
    chart_W->removeAxis(axis_W);
    chart_W->removeAxis(axis_w);
    series_W->clear();

    QString title = "源汇项时域图像";
    ui->graphicsView_W->setChart(chart_W);
    chart_W->setTitle(title);
    series_W->setName("源汇项波动曲线");
    int n = floor(flow.show_tl() / flow.show_st()) + 1;
    double t = 0.0;
    double min_h = 100.0;
    double max_h = 0.0;
    for(int i = 0; i < n; i++){
        double h =flow.source_sink_term(t);
        series_W->append(t, h);
        t += flow.show_st();
        if(h < min_h) min_h = h;
        if(h > max_h) max_h = h;
    }

    axis_w->setRange(0, flow.show_tl());
    axis_w->setGridLineVisible(false);
    axis_w->setLabelFormat("%.1f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("T轴(d)");

    axis_W->setRange(min_h, max_h);
    axis_W->setGridLineVisible(false);
    axis_W->setLabelFormat("%.5f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("水头(m)");

    chart_W->addSeries(series_W); // 更新图表
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::on_frequency_field_figure_clicked()
{
    chart_W->removeSeries(series_W); // 清除原有图表
    chart_W->removeAxis(axis_W);
    chart_W->removeAxis(axis_w);
    series_W->clear();

    QString title = "源汇项频域图像";
    ui->graphicsView_W->setChart(chart_W);
    chart_W->setTitle(title);
    series_W->setName("源汇项频谱曲线");

    flow.set_n_m();
    Eigen::VectorXd Amplitude(flow.show_n());
    double t = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < flow.show_n(); i++){
        Amplitude[i] = flow.source_sink_term(t);
        t += flow.show_st();
    }

    // 快速傅里叶变换振幅
    Eigen::VectorXd Amplitude_fft = flow.fast_fourier_transfrom(Amplitude, flow.show_n());

    // 找到最大振幅值，以便设置绘图坐标轴
    double max_A = Amplitude_fft.maxCoeff();
    ui->graphicsView_W->setChart(chart_W); // 这段必须放在前面
    ui->textBrowser->append("最大振幅为:" + QString::number(max_A));

    double w = 0.0;
    for (double W : Amplitude_fft) {
        series_W->append(w, W);
        w += (1/ flow.show_tl());
    }

    // 绘图
    axis_w->setRange(0, 1 / (2*flow.show_st())); // 转换为取一半的频率坐标
    axis_w->setLabelFormat("%.2f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("频率(次/天)");

    axis_W->setRange(0, max_A);
    axis_W->setLabelFormat("%.5f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("振幅(m)");

    chart_W->addSeries(series_W); // 更新图表
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::on_draw_solve_line_location_fft_clicked()
{
    clear_chart_head();

    QString title = "数值解，空间差分步长为" +QString::number(flow.show_sl()) + "时间差分步长为" + QString::number(flow.show_st()) + "，绘图位置为第" + QString::number(ui->spinBox_X->value()) + "位置。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    series_head->setName("水头曲线: 数值解(频域)");

    int a = ui->spinBox_X->value();
    // 快速傅里叶变换振幅
    Eigen::VectorXd Amplitude_fft = flow.fast_fourier_transfrom(solve_fdm.col(a), flow.show_n());
    // 找到最大振幅值，以便设置绘图坐标轴
    double max_A = Amplitude_fft.maxCoeff();
    double w = 0.0;
    for (double W : Amplitude_fft) {
        series_head->append(w, W);
        w += (1/ flow.show_tl());
    }

    axis_x->setRange(0,  1 / (2*flow.show_st())); // 转换为取一半的频率坐标
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("频率(次/天)");

    axis_head->setRange(0, max_A);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("振幅(m)");

    chart_head->addSeries(series_head); // 更新图表
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);

}

void Random_flow_Window::on_random_new_wave_x_clicked()
{
    if(ui->spinBox_x_length->value() != 0) flow.x_length(ui->spinBox_x_length->value());
    flow.random_source_sink_term_x();
    std::vector<std::vector<double>> list_source_sink_term = flow.share_list_source_sink_term_x();
    QString str = QString::number(ui->doubleSpinBox_rain_x->value()/ (1000 * 365));
    for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
        if(list_source_sink_term[i][0] == 0){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "sin((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * x)";
        }
        else if(list_source_sink_term[i][0] == 1){
            str += " + " + QString::number(list_source_sink_term[i][1]) + "cos((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * x)";
        }
    };
    ui->textBrowser_rain_function_x->clear();
    ui->textBrowser_rain_function_x->append(str);
}

void Random_flow_Window::on_delete_wave_x_clicked()
{
    if(flow.share_list_source_sink_term_x().size() == 0) return;
    else{
        flow.del_last_list_source_sink_term_x();
        std::vector<std::vector<double>> list_source_sink_term = flow.share_list_source_sink_term_x();
        QString str = QString::number(ui->doubleSpinBox_rain_x->value()/ (1000 * 365));
        for (unsigned long long i = 0; i < list_source_sink_term.size(); i++) {
            if(list_source_sink_term[i][0] == 0){
                str += " + " + QString::number(list_source_sink_term[i][1]) + "sin((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * x)";
            }
            else if(list_source_sink_term[i][0] == 1){
                str += " + " + QString::number(list_source_sink_term[i][1]) + "cos((2pi / " + QString::number(list_source_sink_term[i][2]) + ") * x)";
            }
        };
        ui->textBrowser_rain_function_x->clear();
        ui->textBrowser_rain_function_x->append(str);
    }
}

void Random_flow_Window::on_draw_solve_line_fft_clicked()
{
    clear_chart_head();

    QString title = "数值解，空间差分步长为" +QString::number(flow.show_sl()) + "时间差分步长为" + QString::number(flow.show_st()) + "，绘图时刻为第" + QString::number(ui->spinBox_time->value()) + "时刻。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    series_head->setName("水头曲线: 数值解(频域)");

    int a = ui->spinBox_time->value();
    // 快速傅里叶变换振幅
    Eigen::VectorXd Amplitude_fft = flow.fast_fourier_transfrom(solve_fdm.row(a), flow.show_m());
    // 找到最大振幅值，以便设置绘图坐标轴
    double max_A = Amplitude_fft.maxCoeff();
    double w = 0.0;
    for (double W : Amplitude_fft) {
        series_head->append(w, W);
        w += (1/ flow.show_xl());
    }

    axis_x->setRange(0,  1 / (2*flow.show_sl())); // 转换为取一半的频率坐标
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("频率(次/天)");

    axis_head->setRange(0, max_A);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("振幅(m)");

    chart_head->addSeries(series_head); // 更新图表
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);
}

void Random_flow_Window::on_time_field_figure_x_clicked()
{
    chart_W->removeSeries(series_W); // 清除原有图表
    chart_W->removeAxis(axis_W);
    chart_W->removeAxis(axis_w);
    series_W->clear();

    QString title = "源汇项(空间随机)时域图像";
    ui->graphicsView_W->setChart(chart_W);
    chart_W->setTitle(title);
    series_W->setName("源汇项波动曲线");
    int n = floor(flow.show_xl() / flow.show_sl()) + 1;
    double x = 0.0;
    double min_h = 100.0;
    double max_h = 0.0;
    for(int i = 0; i < n; i++){
        double h =flow.source_sink_term_x(x);
        series_W->append(x, h);
        x += flow.show_sl();
        if(h < min_h) min_h = h;
        if(h > max_h) max_h = h;
    }

    axis_w->setRange(0, flow.show_xl());
    axis_w->setGridLineVisible(false);
    axis_w->setLabelFormat("%.1f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("X轴(m)");

    axis_W->setRange(min_h, max_h);
    axis_W->setGridLineVisible(false);
    axis_W->setLabelFormat("%.5f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("水头(m)");

    chart_W->addSeries(series_W); // 更新图表
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::on_frequency_field_figure_x_clicked()
{
    chart_W->removeSeries(series_W); // 清除原有图表
    chart_W->removeAxis(axis_W);
    chart_W->removeAxis(axis_w);
    series_W->clear();

    QString title = "源汇项(空间随机)频域图像";
    ui->graphicsView_W->setChart(chart_W);
    chart_W->setTitle(title);
    series_W->setName("源汇项频谱曲线");

    flow.set_n_m();
    Eigen::VectorXd Amplitude(flow.show_m());
    double x = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < flow.show_m(); i++){
        Amplitude[i] = flow.source_sink_term_x(x);
        x += flow.show_sl();
    }

    // 快速傅里叶变换振幅
    Eigen::VectorXd Amplitude_fft = flow.fast_fourier_transfrom(Amplitude, flow.show_m());

    // 找到最大振幅值，以便设置绘图坐标轴
    double max_A = Amplitude_fft.maxCoeff();
    ui->graphicsView_W->setChart(chart_W); // 这段必须放在前面
    ui->textBrowser->append("最大振幅为:" + QString::number(max_A));

    double w = 0.0;
    for (double W : Amplitude_fft) {
        series_W->append(w, W);
        w += (1/ flow.show_xl());
    }

    // 绘图
    axis_w->setRange(0, 1 / (2*flow.show_sl())); // 转换为取一半的频率坐标
    axis_w->setLabelFormat("%.2f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("频率(次/天)");

    axis_W->setRange(0, max_A);
    axis_W->setLabelFormat("%.5f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("振幅(m)");

    chart_W->addSeries(series_W); // 更新图表
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::on_move_clicked()
{
    QFlags<QChartView::RubberBand> flags = QChartView::ClickThroughRubberBand; // 必须有ClickThroughRubberBand才能把信号传递给序列
    if (ui->move->isChecked())
    {
        move_the_chart = true;
        ui->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);

        if(ui->select->isChecked())  // 两个按钮的互斥代码
        {
            ui->select->setChecked(false);
            select_the_chart = false;
            ui->graphicsView->viewport()->unsetCursor();
            ui->graphicsView->setCursor(Qt::ArrowCursor);
            flags |= QChartView::NoRubberBand;
            ui->graphicsView->setRubberBand(flags);
        }

    }
    else
    {
        move_the_chart = false;
        ui->graphicsView->setDragMode(QGraphicsView::NoDrag);
    }
}

void Random_flow_Window::on_frame_destroyed()
{

}

void Random_flow_Window::on_graphicsView_destroyed()
{

}

void Random_flow_Window::on_select_clicked()
{
    //原因是QGraphicsView继承自QAbstractScrollArea类，只能通过viewport()获取view Widget。
    //所以在QGraphicsView子类中使用viewport()->setCursor()而不是直接setCursor()，这样才能真正改变视觉上的鼠标形状。
    ui->graphicsView->viewport()->unsetCursor();  // 必须加上viewport!!才能更新视图中的鼠标形状！！！坑死我了
    ui->graphicsView->setCursor(Qt::CrossCursor);
    QFlags<QChartView::RubberBand> flags = QChartView::ClickThroughRubberBand; // 必须有ClickThroughRubberBand才能把信号传递给序列
    if(ui->select->isChecked())
    {
        flags |= QChartView::RectangleRubberBand;
        select_the_chart = true;

        if(ui->move->isChecked()) // 两个按钮的互斥代码
        {
            ui->move->setChecked(false);
            move_the_chart = false;
            ui->graphicsView->setDragMode(QGraphicsView::NoDrag);
        }
    }
    else
    {
        ui->graphicsView->viewport()->unsetCursor();
        ui->graphicsView->setCursor(Qt::ArrowCursor);
        flags |= QChartView::NoRubberBand;
        select_the_chart = false;
    }
    ui->graphicsView->setRubberBand(flags);
}

void Random_flow_Window::on_reset_1_clicked()
{
    chart_head->zoomReset();
}

void Random_flow_Window::on_spinBox_angle_valueChanged(double arg1)
{
    flow.set_angle(arg1);
}

void Random_flow_Window::on_doubleSpinBox_angle_valueChanged(double arg1)
{
    flow.set_angle(arg1);
}

void Random_flow_Window::on_use_white_noise_checkBox_time_clicked()
{
    if(ui->use_white_noise_checkBox_time->isChecked() == true)
    {
        flow.set_white_noise_time(1);
    }
    else
    {
        flow.set_white_noise_time(0);
    }
}

void Random_flow_Window::on_actual_expectations_clicked()
{
    QString str = QString::number(flow.actual_expectations_white_noise_time());
    ui->textBrowser->append("在本次随机数生成中，源汇项实际期望值为：" + str);
    QString str1 = QString::number(flow.fangcha_white_noise_time());
    ui->textBrowser->append("在本次随机数生成中，源汇项实际方差为：" + str1);
}

void Random_flow_Window::on_power_spectral_density_figure_clicked()
{
    chart_W->removeSeries(series_W); // 清除原有图表
    chart_W->removeAxis(axis_W);
    chart_W->removeAxis(axis_w);
    series_W->clear();

    QString title = "源汇项功率谱密度图像";
    ui->graphicsView_W->setChart(chart_W);
    chart_W->setTitle(title);
    series_W->setName("源汇项功率谱密度曲线");

    flow.set_n_m();
    Eigen::VectorXd Amplitude(flow.show_n());
    double t = 0.0;
    // 离散化源汇项赋值
    for(int i = 0; i < flow.show_n(); i++){
        Amplitude[i] = flow.source_sink_term(t);
        t += flow.show_st();
    }

    // 快速傅里叶变换振幅
    Eigen::VectorXd Amplitude_psd = flow.power_spectral_density(Amplitude, flow.show_n());

    // 找到最大振幅值，以便设置绘图坐标轴
    double max_A = Amplitude_psd.maxCoeff();
    ui->graphicsView_W->setChart(chart_W); // 这段必须放在前面
    ui->textBrowser->append("最大振幅为:" + QString::number(max_A));

    double w = 0.0;
    for (double W : Amplitude_psd) {
        series_W->append(w, W);
        w += (1/ flow.show_tl());
    }

    // 绘图
    axis_w->setRange(0, 1 / (2*flow.show_st())); // 转换为取一半的频率坐标
    axis_w->setLabelFormat("%.2f"); // 标签格式
    axis_w->setTickCount(6);
    axis_w->setMinorTickCount(1);
    axis_w->setTitleText("频率(次/天)");

    axis_W->setRange(0, max_A);
    axis_W->setLabelFormat("%.5f"); // 标签格式
    axis_W->setTickCount(6);
    axis_W->setMinorTickCount(1);
    axis_W->setTitleText("瓦特每赫兹（W/Hz）");

    chart_W->addSeries(series_W); // 更新图表
    chart_W->addAxis(axis_w, Qt::AlignBottom);
    chart_W->addAxis(axis_W, Qt::AlignLeft);
    series_W->attachAxis(axis_w);
    series_W->attachAxis(axis_W);
}

void Random_flow_Window::on_amplitude_complete_figure_clicked()
{
    clear_chart_head();

    QString title = "功率谱振幅比值，空间差分步长为" +QString::number(flow.show_sl()) + "时间差分步长为" + QString::number(flow.show_st()) + "，绘图位置为第" + QString::number(ui->spinBox_X->value()) + "位置。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    series_head->setName("功率谱振幅比值: 数值解");
    series_analyze->setName("功率谱振幅比值: 解析解");
    int a = ui->spinBox_X->value();
    double x = 0.0;
    double w = 0.0;
    Eigen::VectorXd amplitude_complete_fdm = flow.amplitude_complete_fdm(solve_fdm, a);
    Eigen::VectorXd amplitude_complete_analyze = flow.amplitude_complete_analyze();

    // 找到最大振幅值，以便设置绘图坐标轴
    //double max_A = amplitude_complete_fdm.maxCoeff();
    double max_B = amplitude_complete_analyze.maxCoeff();

    // 使用蒙特卡洛绘图
    if (use_MC_to_draw){
        for (double A_fdm : MC_amplitude_complete_fdm) {
            series_head->append(x,A_fdm);
            x += (1/ flow.show_tl());
        }
    }
    else{
        for (double A_fdm : amplitude_complete_fdm) {
            series_head->append(x,A_fdm);
            x += (1/ flow.show_tl());
        }
    }
    for (double A_analyze : amplitude_complete_analyze) {
        series_analyze->append(w, A_analyze);
        w += (1/ flow.show_tl());
    }

    axis_x->setRange(0, 1 / (2*flow.show_st()));
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("频率(次/天)");

    //if(max_A>max_B)axis_head->setRange(0, max_A);
    //else
    axis_head->setRange(0, max_B);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("振幅比");

    chart_head->addSeries(series_head); // 更新图表
    chart_head->addSeries(series_analyze);
    chart_head->addAxis(axis_x, Qt::AlignBottom);
    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_head->attachAxis(axis_x);
    series_head->attachAxis(axis_head);

}

void Random_flow_Window::on_MC_un_wt_amp_start_clicked() // 关于源汇项随时间成均匀分布的蒙特卡洛线程启动函数
{
    //MC_amplitude_complete_fdm.setZero();
    MCThread_uniform_wt_amp1->set_flow(flow);
    MCThread_uniform_wt_amp1->set_times(floor(MC_times/2));
    MCThread_uniform_wt_amp1->start();
    MCThread_uniform_wt_amp1_work = true;

    MCThread_uniform_wt_amp2->set_flow(flow);
    MCThread_uniform_wt_amp2->set_times(floor(MC_times/2));
    MCThread_uniform_wt_amp2->start();
    MCThread_uniform_wt_amp2_work = true;

    ui->MC_un_wt_amp_start->setEnabled(false);
    ui->progressBar_MC_uwt->setValue(0); // 进度条设置为0

    MC_act_time = 0;
    MC_act_time_1 = 0;
}

void Random_flow_Window::on_spinBox_MC_un_wt_amp_valueChanged(int arg1)
{
    MC_times = arg1;
}

void Random_flow_Window::on_checkBox_use_MC_clicked()
{
    if(ui->checkBox_use_MC->isChecked() == true) use_MC_to_draw = true;
    else use_MC_to_draw = false;
}


void Random_flow_Window::on_draw_solve_line_as_clicked()
{
    clear_chart_head();

    QString title = "解析解与数值解对比图，绘图位置为右边界位置。";
    ui->graphicsView->setChart(chart_head);
    chart_head->setTitle(title);
    //series_head->setName("水头曲线: 数值解");
    series_analyze->setName("水头曲线: 解析解");
    //series_head->setColor(QColorConstants::Red);
    //std::vector<double>solve_fdm_l(flow.show_n());
//    int a = ui->spinBox_X->value();
//    int x_fdm = 0;
//    for (double h : solve_fdm.col(a)) {
//        solve_fdm_l[x_fdm] = h;
//        x_fdm += 1;
//    }
//    for (int n = 0; n < flow.show_n(); n++) {
//        series_head->append(n, solve_fdm_l[n]);
//    }

    double x = 0.0;
    double min_h = 800.0;
    double max_h = 0.0;
    Eigen::MatrixXd solve_an = flow.solve_an_wt();
    for (double h : solve_an.col((flow.show_m()-1))) {
        series_analyze->append(x, h);
        x += flow.show_st();
        if(h < min_h) min_h = h;
        if(h > max_h) max_h = h;
    }

    axis_x->setRange(0, flow.show_tl());
    axis_x->setLabelFormat("%.2f"); // 标签格式
    axis_x->setTickCount(11);
    axis_x->setMinorTickCount(1);
    axis_x->setTitleText("时间轴(d)");

    chart_head->addAxis(axis_x, Qt::AlignBottom);
    series_analyze->attachAxis(axis_x);

    axis_head->setRange(min_h, max_h);
    axis_head->setLabelFormat("%.4f"); // 标签格式
    axis_head->setTickCount(11);
    axis_head->setMinorTickCount(1);
    axis_head->setTitleText("水头(m)");

    chart_head->addAxis(axis_head, Qt::AlignLeft);
    series_analyze->attachAxis(axis_head);

    // 更新图表
    chart_head->addSeries(series_analyze);
//    series_head->attachAxis(axis_x);
//    series_head->attachAxis(axis_head);
//    chart_head->addSeries(series_head);
}

void Random_flow_Window::on_actionsave_fdm_triggered() // 存储当前数值解数据
{
    QString curPath = QDir::currentPath();  // 获取应用程序当前目录
    QString dlgTitle = "保存当前数值解数据";
    QString filter = "文本文件(*.txt)";
    QString aFileName = QFileDialog::getSaveFileName(this, dlgTitle, curPath, filter);
    ui->textBrowser->append(aFileName);
    save_the_data(aFileName, 0);
}

void Random_flow_Window::on_actionsave_MC_triggered()
{
    QString curPath = QDir::currentPath();  // 获取应用程序当前目录
    QString dlgTitle = "保存当前蒙特卡洛数据";
    QString filter = "文本文件(*.txt)";
    QString aFileName = QFileDialog::getSaveFileName(this, dlgTitle, curPath, filter);
    ui->textBrowser->append(aFileName);
    save_the_data(aFileName, 1);
}

void Random_flow_Window::on_actionsave_as_triggered()
{
    QString curPath = QDir::currentPath();  // 获取应用程序当前目录
    QString dlgTitle = "保存当前解析解数据";
    QString filter = "文本文件(*.txt)";
    QString aFileName = QFileDialog::getSaveFileName(this, dlgTitle, curPath, filter);
    ui->textBrowser->append(aFileName);
    save_the_data(aFileName, 2);
}

void Random_flow_Window::on_actionsave_as_complete_triggered()
{
    QString curPath = QDir::currentPath();  // 获取应用程序当前目录
    QString dlgTitle = "保存当前解析解数据(振幅比)";
    QString filter = "文本文件(*.txt)";
    QString aFileName = QFileDialog::getSaveFileName(this, dlgTitle, curPath, filter);
    ui->textBrowser->append(aFileName);
    save_the_data(aFileName, 3);
}

void Random_flow_Window::on_doubleSpinBox_left_boundary_valueChanged(double arg1)
{
    QString str = QString::number(arg1);
    ui->textBrowser_function_h_l->clear();  // 重置展示边界条件函数的框
    flow.clear_list_h_l();  // 重置左边界波动列表
    ui->textBrowser_function_h_l->append(str);
}

void Random_flow_Window::on_doubleSpinBox_right_boundary_valueChanged(double arg1)
{
    QString str = QString::number(arg1);
    ui->textBrowser_function_h_r->clear();  // 重置展示边界条件函数的框
    flow.clear_list_h_r();  // 重置右边界波动列表
    ui->textBrowser_function_h_r->append(str);
}

void Random_flow_Window::on_random_new_wave_h_l_clicked()
{
    if(ui->spinBox_t_length->value() != 0) flow.t_length(ui->spinBox_t_length->value());
    flow.random_h_l();
    std::vector<std::vector<double>> list_h_l = flow.share_list_h_l();
    QString str = QString::number(ui->doubleSpinBox_left_boundary->value());
    for (unsigned long long i = 0; i < list_h_l.size(); i++) {
        if(list_h_l[i][0] == 0){
            str += " + " + QString::number(list_h_l[i][1]) + "sin((2pi / " + QString::number(list_h_l[i][2]) + ") * t)";
        }
        else if(list_h_l[i][0] == 1){
            str += " + " + QString::number(list_h_l[i][1]) + "cos((2pi / " + QString::number(list_h_l[i][2]) + ") * t)";
        }
    };
    ui->textBrowser_function_h_l->clear();
    ui->textBrowser_function_h_l->append(str);
}

void Random_flow_Window::on_delete_wave_h_l_clicked()
{
    if(flow.share_list_h_l().size() == 0) return;
    else{
        flow.del_last_list_h_l();
        std::vector<std::vector<double>> list_h_l = flow.share_list_h_l();
        QString str = QString::number(ui->doubleSpinBox_left_boundary->value());
        for (unsigned long long i = 0; i < list_h_l.size(); i++) {
            if(list_h_l[i][0] == 0){
                str += " + " + QString::number(list_h_l[i][1]) + "sin((2pi / " + QString::number(list_h_l[i][2]) + ") * t)";
            }
            else if(list_h_l[i][0] == 1){
                str += " + " + QString::number(list_h_l[i][1]) + "cos((2pi / " + QString::number(list_h_l[i][2]) + ") * t)";
            }
        };
        ui->textBrowser_function_h_l->clear();
        ui->textBrowser_function_h_l->append(str);
    }
}

void Random_flow_Window::on_random_new_wave_h_r_clicked()
{
    if(ui->spinBox_t_length->value() != 0) flow.t_length(ui->spinBox_t_length->value());
    flow.random_h_r();
    std::vector<std::vector<double>> list_h_r = flow.share_list_h_r();
    QString str = QString::number(ui->doubleSpinBox_right_boundary->value());
    for (unsigned long long i = 0; i < list_h_r.size(); i++) {
        if(list_h_r[i][0] == 0){
            str += " + " + QString::number(list_h_r[i][1]) + "sin((2pi / " + QString::number(list_h_r[i][2]) + ") * t)";
        }
        else if(list_h_r[i][0] == 1){
            str += " + " + QString::number(list_h_r[i][1]) + "cos((2pi / " + QString::number(list_h_r[i][2]) + ") * t)";
        }
    };
    ui->textBrowser_function_h_r->clear();
    ui->textBrowser_function_h_r->append(str);
}

void Random_flow_Window::on_delete_wave_h_r_clicked()
{
    if(flow.share_list_h_r().size() == 0) return;
    else{
        flow.del_last_list_h_r();
        std::vector<std::vector<double>> list_h_r = flow.share_list_h_r();
        QString str = QString::number(ui->doubleSpinBox_left_boundary->value());
        for (unsigned long long i = 0; i < list_h_r.size(); i++) {
            if(list_h_r[i][0] == 0){
                str += " + " + QString::number(list_h_r[i][1]) + "sin((2pi / " + QString::number(list_h_r[i][2]) + ") * t)";
            }
            else if(list_h_r[i][0] == 1){
                str += " + " + QString::number(list_h_r[i][1]) + "cos((2pi / " + QString::number(list_h_r[i][2]) + ") * t)";
            }
        };
        ui->textBrowser_function_h_r->clear();
        ui->textBrowser_function_h_r->append(str);
    }
}

