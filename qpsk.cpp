#include<bits/stdc++.h>

using namespace std;
using Complex= complex<double>; //creting complex alias
constexpr double INV = 0.707;

//QPSK modulator function
Complex qpskmodulator(int b1, int b2){
    if ((b1 !=0 && b1 != 1)||(b2 !=0 && b2 !=1)){
    throw invalid_argument("out of range");}
    double I=(b1==0 ? 1.0:-1.0)*INV;
    double Q=(b2==0 ? 1.0:-1.0)*INV;
    return Complex(I,Q);
}

vector<int> qpskDemodulator(const Complex& symbol){
   int b1 = (symbol.real() <0.0) ? 1:0;
   int b2 = (symbol.imag() < 0.0) ? 1:0;
   return {b1,b2};
}

int main(){
    vector<int> bits = {0,1,1,1,0,0,0,0};
    vector<Complex> symbols;
    for (size_t i=0;i<bits.size();i+=2){
        Complex symbol = qpskmodulator(bits[i],bits[i+1]);
        symbols.push_back(symbol);
    }
    cout << "QPSK Transmitter";
vector<int> bit;
for (int bit:bits){
    cout << bit<< ' ';
}
for (size_t i = 0; i < symbols.size(); ++i) {
    cout << "Symbol " << i << ": "
         << symbols[i] << '\n';
}

vector<int> recoveredBits;
for (const Complex& symbol : symbols) {
    vector<int> pair = qpskDemodulator(symbol);
    recoveredBits.push_back(pair[0]);
    recoveredBits.push_back(pair[1]);
}

for (int bit:recoveredBits){
    cout << bit << ' ';
}
}