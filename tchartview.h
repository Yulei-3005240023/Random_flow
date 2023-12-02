#ifndef TCHARTVIEW_H
#define TCHARTVIEW_H

#include <QtCharts>
class TChartView : public QChartView
{
    Q_OBJECT
private:
    QPoint beginPoint; // 起点
    QPoint endPoint;  // 终点
   // bool m_custom_Zoom = false; // 是否使用自定义矩形放大模式
protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
public:
    TChartView(QWidget *parent = nullptr);
signals:
    void mouseMovePoint(QPoint point);
};

#endif // TCHARTVIEW_H
