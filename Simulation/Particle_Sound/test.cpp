#include <iostream>
#include <fstream>

using namespace std;

int main(){

    ifstream file("test.txt");

    char c;
    double x;
    while (!file.eof()){
        if (file.get() == ':'){
            file >> x;
            cout << x << endl;
        }
    }

    return 0;
}