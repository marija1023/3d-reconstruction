#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLabel>
#include <QPushButton>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QMouseEvent>
#include <vector>
#include "Eigen/Dense"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public slots:
    void swap();
    void check();

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void mousePressEvent(QMouseEvent *evt);

    std::vector<std::tuple<double, double, double>> getLeftPoints() const
    {
        return leftPoints;
    }

    std::vector<std::tuple<double, double, double>> getRightPoints() const
    {
        return rightPoints;
    }

private:
    Ui::MainWindow *ui;
    bool leftImg;
    std::vector<std::tuple<double, double, double>> leftPoints;
    std::vector<std::tuple<double, double, double>> rightPoints;

    Eigen::MatrixXd left;
    Eigen::MatrixXd right;
};

#endif // MAINWINDOW_H
