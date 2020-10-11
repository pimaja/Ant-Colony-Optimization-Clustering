#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "1cpp.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    r=0; mi=0;
    bestfit="", podaci="", kmfit="";
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_lineEdit_textEdited(const QString &arg1)
{
    r=arg1.toInt();
}

void MainWindow::on_lineEdit_2_textEdited(const QString &arg1)
{
    mi=arg1.toInt();
}

/*void MainWindow::on_lineEdit_4_textEdited(const QString &arg1)
{
    podaci=arg1.toStdString();
}*/

void MainWindow::on_pushButton_released()
{

    klasa KLASA (r, mi, ui->comboBox->currentText().toStdString());
    KLASA.main2(bfit, kfit);
    bestfit=QString::number(bfit);
    kmfit=QString::number(kfit);
    ui->lineEdit_3->setText(bestfit);
    ui->lineEdit_4->setText(kmfit);
}

