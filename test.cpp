#include "piccoulomb.h"

int main() {
    MatrixXf p = MatrixXf::Zero(10, 10);
    Matrix2d a = Matrix2d::Identity();
    Matrix2d b = a;
    // Vector2d v = b.block(0,0,2,1);
    // Vector2d v = b(seq(0,1),seq(0,0));
    Vector2d v = b.col(0);

    MatrixXf lol = MatrixXf::Identity(100,100);
    cout << "asd" << endl;
    VectorXf vlol = lol.col(0);
    cout << "asd asd" << endl;
    lol(0,0) = 6;

    a(0,0) = 5;
    b(1,0) = 7;
    

    cout << "p " << p << endl << endl;
    cout << a << endl << endl;
    cout << b << endl << endl;
    cout << v << endl << endl;
    cout << lol(0,0) << endl << endl;
    cout << "asdasd " << vlol[0] << endl << endl;

    return 0;
}