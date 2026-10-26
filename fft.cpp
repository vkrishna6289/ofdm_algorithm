#include <bits/stdc++.h>
using namespace std;

using Complex = complex<double>;
constexpr double PI = 3.14159265358979323846;

// Inverse DFT
vector<Complex> idft(const vector<Complex>& input) {
int N = input.size();
vector<Complex> output(N);

for (int n = 0; n < N; n++) {
    Complex sum(0.0, 0.0);

    for (int k = 0; k < N; k++) {
        double angle = 2.0 * PI * k * n / N;
        Complex W(cos(angle), sin(angle));
        sum += input[k] * W;
    }

    output[n] = sum / double(N);
}

return output;

}

// Forward DFT
vector<Complex> dft(const vector<Complex>& input) {
int N = input.size();
vector<Complex> output(N);

for (int k = 0; k < N; k++) {
    Complex sum(0.0, 0.0);

    for (int n = 0; n < N; n++) {
        double angle = -2.0 * PI * k * n / N;
        Complex W(cos(angle), sin(angle));
        sum += input[n] * W;
    }

    output[k] = sum;
}
return output;
}

int main(){
    vector<Complex> input={{2,0},{3,4},{4,6},{5.6}};
    vector<Complex> frequency=dft(input);
    //dft();

    vector<Complex> reconstucted=idft(frequency);
    cout << "DFT OUTPUT:\n";

    for (const Complex& x:frequency){
        cout << x << '\n';
    }
    cout << "IDFT OUTPUT:\n";

    for (const Complex& x:reconstucted){
        cout << x << '\n';
    }
    return 0;

}