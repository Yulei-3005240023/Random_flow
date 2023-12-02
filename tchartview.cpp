#include "tchartview.h"

TChartView::TChartView(QWidget *parent):QChartView(parent)
{
    this->setMouseTracking(true);
}

void TChartView::mousePressEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton && this->dragMode() == QGraphicsView::ScrollHandDrag)
    {
        beginPoint = event->pos();
    }
    QChartView::mousePressEvent(event);
}

void TChartView::mouseReleaseEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton && this->dragMode() == QGraphicsView::ScrollHandDrag)
    {
        endPoint = event->pos();
        chart()->scroll(beginPoint.x() - endPoint.x(), endPoint.y() - beginPoint.y());
    }
    QChartView::mouseReleaseEvent(event);
}

void TChartView::mouseMoveEvent(QMouseEvent *event)
{// 鼠标移动事件
    QPoint point = event->pos();
    emit mouseMovePoint(point);
    QChartView::mouseMoveEvent(event); // 父类继续处理这件事，必须这么写
}
