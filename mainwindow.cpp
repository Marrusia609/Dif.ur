#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <vector>
using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // структура объявляется в заголовочном


    double  a = 0.0,    // начальное время
            b = 10.0;    // конечное время
    double  h = 0.02;   // шаг
    int   u0[] = {100, 50}; // количество в начале

    vector<Data> data;
    Data start;
    start.t = a;
    start.x = u0[0];
    start.y = u0[1];
    data.push_back(start);

    int stepCount = (b / h) + 1;
    double currentTime = a + h;

    double  x, y;
    double  x_prev = start.x,
            y_prev = start.y;
    double K[4];
    for (int j = 1; j <= stepCount; j++)
    {
        // Решаем X
        K[0] = (2.0-0.02*y_prev) * x_prev;
        K[1] = (2.0-0.02*(y_prev + h*K[0]/2.0))  * (x_prev + h*K[0]/2.0);
        K[2] = (2.0-0.02*(y_prev + h*K[1]/2.0))  * (x_prev + h*K[1]/2.0);
        K[3] = (2.0-0.02*(y_prev + h*K[2]))    * (x_prev + h*K[2]);
        x   = x_prev + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решаем Y
        K[0] = (-1.0*y_prev  + 0.01*x_prev*y_prev);
        K[1] = (-1.0*(y_prev + h*K[0]/2.0) + 0.01*(x_prev + h*K[0]/2.0)   * (y_prev + h*K[0]/2.0));
        K[2] = (-1.0*(y_prev + h*K[1]/2.0) + 0.01*(x_prev + h*K[1]/2.0)   * (y_prev + h*K[1]/2.0));
        K[3] = (-1.0*(y_prev + h*K[2]/2.0) + 0.01*(x_prev + h*K[2])     * (y_prev + h*K[2]));
        y   = y_prev + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        Data currentValues;
        currentValues.x = x;
        currentValues.y = y;
        currentValues.t = currentTime;
        data.push_back(currentValues);

        x_prev = x;
        y_prev = y;

        currentTime += h;
    }

    for (Data dt : data)
    {
        cout << dt.t << "\t" <<  dt.x << "\t" << dt.y << endl;
    }

    cout << "Готово! Кол-во шагов: " << stepCount << endl;
    cout << "Размер: " << sizeof(Data)*data.size() << " bytes" << endl;
    //return 0;

     // решение адаптированно, изначальная задача проверялась с Matlab-овским ode45 (Тоже Рунге-Кутты метод).
    //Значения немного отличаются. Под конец разница на 15(и X и Y)
    // однако эьто объясняется ошибкой в справочной задаче в формулах записи
    // при переработке для наших уравнений, разницы в результатах быть не должно
}


MainWindow::~MainWindow()
{
    delete ui;
}

