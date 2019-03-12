#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <Eigen/Geometry>

typedef Eigen::Matrix<double, 9, 1> Vector9d;

void MainWindow::swap()
{
    leftImg = !leftImg;
    if(leftImg) {
        // left
        ui->lblImage->clear();
        ui->lblImage->setPixmap(QPixmap(":/images/leva_orig.JPG"));
    } else {
        // right
        ui->lblImage->clear();
        ui->lblImage->setPixmap(QPixmap(":/images/desna_orig.JPG"));
    }
    ui->lblImage->setScaledContents(true);
    ui->lblImage->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);

}

Eigen::Vector3d af_coords(const Eigen::Vector3d& e)
{
    return 1.0/e(2) * e;
}

Eigen::Matrix3d points_normalization(const Eigen::MatrixXd& m)
{
    double Tx = 0;
    double Ty = 0;
    for(unsigned i = 0; i < m.cols(); i++) {
        Tx += m(0, i);
        Ty += m(1, i);
    }
    Tx = Tx / m.cols();
    Ty = Ty / m.cols();
    std::cout << "Tx" << std::endl;
    std::cout << Tx << std::endl;
    std::cout << "Ty" << std::endl;
    std::cout << Ty << std::endl;

    // translira se teziste u koordinatni pocetak
    Eigen::Matrix3d G;
    G << 1, 0, -Tx,
         0, 1, -Ty,
         0, 0, 1;
    std::cout << "G" << std::endl;
    std::cout << G << std::endl;

    // transliraju se tacke
    Eigen::MatrixXd m_trans;
    m_trans.resize(3, m.cols());
    for(unsigned i = 0; i < m.cols(); i++) {
        m_trans.col(i) = G * m.col(i);
    }
    std::cout << "m_trans" << std::endl;
    std::cout << m_trans << std::endl;

    // prosecno rasojanje od koordinatnog pocetka
    double avg = 0;
    for(unsigned i = 0; i < m_trans.cols(); i++) {
        avg += sqrt(m_trans(0, i)*m_trans(0, i) + m_trans(1, i)*m_trans(1, i));
    }
    avg = avg / m_trans.cols();
    std::cout << "avg" << std::endl;
    std::cout << avg << std::endl;

    // skaliraju se tacke tako da udaljenost tacke od koordinatnog pocetka budu sqrt(2)
    Eigen::Matrix3d S;
    S << sqrt(2)/avg, 0, 0,
         0, sqrt(2)/avg, 0,
         0, 0, 1;
    std::cout << "S" << std::endl;
    std::cout << S << std::endl;

    return S * G;
}

Eigen::Matrix3d normalized_xpoint_algorithm(Eigen::MatrixXd& x, Eigen::MatrixXd& xp)
{
    Eigen::Matrix3d T = points_normalization(x);
    Eigen::Matrix3d Tp = points_normalization(xp);

    Eigen::MatrixXd x_norm, xp_norm;
    x_norm.resize(3, x.cols());
    xp_norm.resize(3, xp.cols());
    for(unsigned i = 0; i < x.cols(); i++) {
        x_norm.col(i) = T * x.col(i);
        xp_norm.col(i) = Tp * xp.col(i);
    }

    std::cout << "x_norm" << std::endl;
    std::cout << x_norm << std::endl;
    std::cout << "xp_norm" << std::endl;
    std::cout << xp_norm << std::endl;

    Eigen::MatrixXd A(x_norm.cols(), 9);

    for(unsigned i = 0; i < x_norm.cols(); i++) {
        A.row(i) << xp_norm(0, i)*x_norm(0, i), xp_norm(0, i)*x_norm(1, i), xp_norm(0, i),
                xp_norm(1, i)*x_norm(0, i), xp_norm(1, i)*x_norm(1, i), xp_norm(1, i),
                x_norm(0, i), x_norm(1, i), 1;
    }
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;

    Eigen::Matrix3d F(3, 3);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto V = svd.matrixV();
    Vector9d f = V.col(V.cols() - 1);
    F << f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8];
    std::cout << "F" << std::endl;
    std::cout << F << std::endl;
    std::cout << "det(F): " << F.determinant() << std::endl;

    std::cout << "provera" << std::endl;
    for(unsigned i = 0; i < xp_norm.cols(); i++) {
        std::cout << xp_norm.col(i).transpose() * F * x_norm.col(i) << std::endl;
    }

    Eigen::JacobiSVD<Eigen::Matrix3d> svd_F(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto DD = svd_F.singularValues();
    Eigen::Matrix3d DDMatrix;
    DDMatrix << DD(0), 0, 0,
                0, DD(1), 0,
                0, 0, 0;

    Eigen::Matrix3d F1 = svd_F.matrixU() * DDMatrix * svd_F.matrixV().transpose();
    std::cout << "F1" << std::endl;
    std::cout << F1 << std::endl;

    std::cout << "provera" << std::endl;
    for(unsigned i = 0; i < xp_norm.cols(); i++) {
        std::cout << xp_norm.col(i).transpose() * F1 * x_norm.col(i) << std::endl;
    }
    std::cout << "det(F1): " << F1.determinant() << std::endl;

    return Tp.transpose() * F1 * T;
}

void camera_info(const Eigen::MatrixXd& cam)
{
    std::cout << "***INFO***" << std::endl;

    Eigen::Matrix3d T0;
    T0 << cam(0, 0), cam(0, 1), cam(0, 2),
          cam(1, 0), cam(1, 1), cam(1, 2),
          cam(2, 0), cam(2, 1), cam(2, 2);
    Eigen::Vector3d T0C = cam.col(3);
    Eigen::Vector3d C;
    C = T0.inverse() * T0C;
    std::cout << "C" << std::endl;
    std::cout << C << std::endl;

//    Eigen::Matrix3d P;
//    P << 0, 0, 1,
//         0, 1, 0,
//         1, 0, 0;
//    Eigen::Matrix3d T0_;
//    T0_ = P * T0;

//    Eigen::HouseholderQR<Eigen::Matrix3d> qr(T0_.transpose());
    Eigen::HouseholderQR<Eigen::Matrix3d> qr(T0.inverse());
    Eigen::Matrix3d Q = qr.householderQ();
    Eigen::Matrix3d R = qr.matrixQR().triangularView<Eigen::Upper>();
//    Eigen::Matrix3d A = P * Q.transpose();
//    Eigen::Matrix3d K = P * R.transpose() * P;
    Eigen::Matrix3d A = Q.transpose();
    Eigen::Matrix3d K = R.inverse();

//    std::cout << "Q" << std::endl;
//    std::cout << Q << std::endl;
//    std::cout << "R" << std::endl;
//    std::cout << R << std::endl;
//    std::cout << "Q*R" << std::endl;
//    std::cout << Q*R << std::endl;

    std::cout << "K" << std::endl;
    std::cout << K << std::endl;
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;

    std::cout << "***KRAJ INFO***" << std::endl;
}

Eigen::Vector3d intersection_point(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3)
{
    Eigen::Vector3d line01 = p0.cross(p1);
    Eigen::Vector3d line23 = p2.cross(p3);

    return line01.cross(line23);
}

void add_inivisible_points(Eigen::MatrixXd &mat)
{
    Eigen::Vector3d x_1 = mat.col(0);
    Eigen::Vector3d x_2 = mat.col(1);
    Eigen::Vector3d x_3 = mat.col(2);
    Eigen::Vector3d x_4 = mat.col(3);
    Eigen::Vector3d xp_1 = mat.col(4);
    Eigen::Vector3d xp_2 = mat.col(5);
    Eigen::Vector3d xp_3 = mat.col(6);
    Eigen::Vector3d x_5 = mat.col(7);
    Eigen::Vector3d x_6 = mat.col(8);
    Eigen::Vector3d x_7 = mat.col(9);
    Eigen::Vector3d x_8 = mat.col(10);
    Eigen::Vector3d xp_6 = mat.col(11);
    Eigen::Vector3d xp_7 = mat.col(12);

    // nepoznate:
    // xp_4
    // xp_5
    // xp_8

    // sibice - xp_4
    Eigen::Vector3d X1_inf, Y1_inf;
    // {X1_inf} = x_2*x_3 /\ xp_2*xp_3
    X1_inf = intersection_point(x_2, x_3, xp_2, xp_3);
    // {Y1_inf} = x_2*xp_2 /\ x_1*xp_1
    Y1_inf = intersection_point(x_2, xp_2, x_1, xp_1);

    // {xp_4} = xp_1*X1_inf /\ x_4*Y1_inf
    Eigen::Vector3d xp_4;
    xp_4 = intersection_point(xp_1, X1_inf, x_4, Y1_inf);

    // kutija - xp_5, xp_8
    Eigen::Vector3d X2_inf, Y2_inf, Z2_inf;
    // {X2_inf} = x_5*x_6 /\ x_7*x_8
    X2_inf = intersection_point(x_5, x_6, x_7, x_8);
    // {Y2_inf} = x_6*xp_6 /\ x_7*xp_7
    Y2_inf = intersection_point(x_6, xp_6, x_7, xp_7);
    // {Z2_inf} = x_6*x_7 /\ xp_6*xp_7
    Z2_inf = intersection_point(x_6, x_7, xp_6, xp_7);

    // {xp_8} = xp_7*X2_inf /\ x_8*Y2_inf
    Eigen::Vector3d xp_8;
    xp_8 = intersection_point(xp_7, X2_inf, x_8, Y2_inf);

    // {xp_5} = xp_6*X2_inf /\ x_5*Y2_inf
    Eigen::Vector3d xp_5;
    xp_5 = intersection_point(xp_6, X2_inf, x_5, Y2_inf);

    mat.col(13) = xp_4;
    mat.col(14) = xp_5;
    mat.col(15) = xp_8;

    mat.col(13) = mat.col(13) / mat(2, 13);
    mat.col(14) = mat.col(14) / mat(2, 14);
    mat.col(15) = mat.col(15) / mat(2, 15);
}

void MainWindow::check()
{
    if(leftPoints.size() == rightPoints.size() && leftPoints.size() >= 8) {
//        left.resize(3, leftPoints.size());
//        right.resize(3, rightPoints.size());
        left.resize(3, 13);
        right.resize(3, 13);

//        for(unsigned i = 0; i < leftPoints.size(); i++) {
//            qDebug() << i << std::get<0>(leftPoints[i]) << " " << std::get<1>(leftPoints[i]) << " " << std::get<2>(leftPoints[i]);
//            left.col(i) << std::get<0>(leftPoints[i]), std::get<1>(leftPoints[i]), std::get<2>(leftPoints[i]);
//        }
//        for(unsigned i = 0; i < rightPoints.size(); i++) {
//            qDebug() << i << std::get<0>(rightPoints[i]) << " " << std::get<1>(rightPoints[i]) << " " << std::get<2>(rightPoints[i]);
//            right.col(i) << std::get<0>(rightPoints[i]), std::get<1>(rightPoints[i]), std::get<2>(rightPoints[i]);
//        }

//        left.col(0) << 958, 38, 1;
//        left.col(1) << 1117, 111, 1;
//        left.col(2) << 874, 285, 1;
//        left.col(3) << 707, 218, 1;
//        left.col(4) << 292, 569, 1;
//        left.col(5) << 770, 969, 1;
//        left.col(6) << 770, 1465, 1;
//        left.col(7) << 317, 1057, 1;

//        right.col(0) << 933, 33, 1;
//        right.col(1) << 1027, 132, 1;
//        right.col(2) << 692, 223, 1;
//        right.col(3) << 595, 123, 1;
//        right.col(4) << 272, 360, 1;
//        right.col(5) << 432, 814, 1;
//        right.col(6) << 414, 1284, 1;
//        right.col(7) << 258, 818, 1;

        left.col(0) << 343, 234, 1;
        left.col(1) << 372, 245, 1;
        left.col(2) << 369, 317, 1;
        left.col(3) << 344, 307, 1;
        left.col(4) << 386, 204, 1;
        left.col(5) << 414, 215, 1;
        left.col(6) << 413, 284, 1;
        left.col(7) << 270, 293, 1;
        left.col(8) << 355, 362, 1;
        left.col(9) << 353, 444, 1;
        left.col(10) << 274, 376, 1;
        left.col(11) << 479, 297, 1;
        left.col(12) << 473, 379, 1;

        right.col(0) << 349, 232, 1;
        right.col(1) << 365, 247, 1;
        right.col(2) << 361, 316, 1;
        right.col(3) << 346, 298, 1;
        right.col(4) << 410, 217, 1;
        right.col(5) << 425, 232, 1;
        right.col(6) << 418, 297, 1;
        right.col(7) << 295, 271, 1;
        right.col(8) << 322, 347, 1;
        right.col(9) << 318, 426, 1;
        right.col(10) << 293, 349, 1;
        right.col(11) << 474, 327, 1;
        right.col(12) << 464, 405, 1;

        std::cout << "left" << std::endl;
        std::cout << left << std::endl;
        std::cout << "right" << std::endl;
        std::cout << right << std::endl;

        // nevidljive tacke
        left.conservativeResize(Eigen::NoChange, left.cols()+3);
        right.conservativeResize(Eigen::NoChange, right.cols()+3);
        add_inivisible_points(left);
        add_inivisible_points(right);

        std::cout << "left" << std::endl;
        std::cout << left << std::endl;
        std::cout << "right" << std::endl;
        std::cout << right << std::endl;

        Eigen::Matrix3d FF1 = normalized_xpoint_algorithm(left, right);
        std::cout << "FF1" << std::endl;
        std::cout << FF1 << std::endl;

        // epipoli
        Eigen::JacobiSVD<Eigen::Matrix3d> svd_e(FF1, Eigen::ComputeFullU | Eigen::ComputeFullV);
        auto VF = svd_e.matrixV();
        auto UF = svd_e.matrixU();

        Eigen::Vector3d e1 = VF.col(VF.cols() - 1).normalized();
        std::cout << "e1" << std::endl;
        std::cout << e1 << std::endl;
        Eigen::Vector3d e2 = UF.col(UF.cols() - 1).normalized();
        std::cout << "e2" << std::endl;
        std::cout << e2 << std::endl;

        // kao kod Rodrigeza, dobija se matrica vektorskog mnozenja jedinicnim vektorom e2
        Eigen::Matrix3d E2;
        E2 << 0, -e2(2), e2(1),
              e2(2), 0, -e2(0),
              -e2(1), e2(0), 0;
        std::cout << "E2" << std::endl;
        std::cout << E2 << std::endl;
        Eigen::Matrix3d M = E2 * FF1;
        std::cout << "M" << std::endl;
        std::cout << M << std::endl;


        // nekalibrisane matrice kamera
        Eigen::MatrixXd T1;
        T1.resize(3, 4);
        T1 << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0;
        std::cout << "T1" << std::endl;
        std::cout << T1 << std::endl;
        std::cout << "T1 info" << std::endl;
        camera_info(T1);

        Eigen::MatrixXd T2;
        T2.resize(3, 4);
        T2 << E2 * FF1, e2;
        std::cout << "T2" << std::endl;
        std::cout << T2 << std::endl;
        std::cout << "T2 info" << std::endl;
        camera_info(T2);

        // 3D koordinate
        Eigen::MatrixXd X_hom;
        X_hom.resize(4, left.cols());

        for(unsigned i = 0; i < left.cols(); i++) {
            Eigen::Matrix4d sys;
            sys << left(0, i) * T1.row(2) - T1.row(0),
                   -left(1, i) * T1.row(2) + T1.row(1),
                   right(0, i) * T2.row(2) - T2.row(0),
                   -right(1, i) * T2.row(2) + T2.row(1);
            Eigen::JacobiSVD<Eigen::Matrix4d> svd_sys(sys, Eigen::ComputeFullU | Eigen::ComputeFullV);
//            std::cout << "sys " << i << std::endl;
//            std::cout << sys << std::endl;
            auto VX = svd_sys.matrixV();
            X_hom.col(i) = VX.col(VX.cols() - 1);
        }

        std::cout << "X_hom" << std::endl;
        std::cout << X_hom << std::endl;

        Eigen::MatrixXd X;
        X.resize(3, left.cols());
        std::cout << "X" << std::endl;
        std::cout << "{ ";
        for(unsigned i = 0; i < left.cols(); i++) {
            X.col(i) << X_hom(0, i) / X_hom(3, i), X_hom(1, i) / X_hom(3, i), X_hom(2, i) / X_hom(3, i);
            std::cout << "{" << X(0, i) << ", " << X(1, i) << ", " << X(2, i) << "}, ";
        }
        std::cout << " }" << std::endl;

        Eigen::MatrixXd X_sort;
        X_sort.resize(3, left.cols());

        X_sort.col(0) = X.col(4);
        X_sort.col(1) = X.col(5);
        X_sort.col(2) = X.col(1);
        X_sort.col(3) = X.col(0);
        X_sort.col(4) = X.col(13);
        X_sort.col(5) = X.col(6);
        X_sort.col(6) = X.col(2);
        X_sort.col(7) = X.col(3);
        X_sort.col(8) = X.col(7);
        X_sort.col(9) = X.col(8);
        X_sort.col(10) = X.col(9);
        X_sort.col(11) = X.col(10);
        X_sort.col(12) = X.col(14);
        X_sort.col(13) = X.col(11);
        X_sort.col(14) = X.col(12);
        X_sort.col(15) = X.col(15);

        std::cout << "X sortirano" << std::endl;
        std::cout << "{ ";
        for(unsigned i = 0; i < left.cols(); i++) {
            std::cout << "{" << X_sort(0, i) << ", " << X_sort(1, i) << ", " << X_sort(2, i) << "}, ";
        }
        std::cout << " }" << std::endl;
    }
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    leftImg(true)
{
    ui->setupUi(this);

    ui->lblImage->setPixmap(QPixmap(":/images/leva_orig.JPG"));
    ui->lblImage->setScaledContents(true);
    ui->lblImage->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);

    connect(ui->btnSwap, SIGNAL(clicked(bool)),
            this, SLOT(swap()));

    connect(ui->btnDone, SIGNAL(clicked(bool)),
            this, SLOT(check()));
//    view->show();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::mousePressEvent(QMouseEvent *evt)
{
    if(evt->button() == Qt::LeftButton) {
//        qDebug() << evt->pos();
        QPointF e(mapToParent(evt->pos()));
        if(leftImg) {
            leftPoints.push_back(std::make_tuple<double, double, double>(e.x(), e.y(), 1));
        } else {
            rightPoints.push_back(std::make_tuple<double, double, double>(e.x(), e.y(), 1));
        }
        qDebug() << leftPoints.size() << " " << rightPoints.size();
    }
}
