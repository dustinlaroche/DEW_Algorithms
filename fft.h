#include <complex>
#include <vector>
#include <cmath>


vector<Complex> fft(vector<double> a, int transformType) {
    int e = a.size();
    if (e == 1) {
      return vector<Complex>(1, a[0]);
    }

    Complex wn = exp(2 * M_PI * complex<double>(0, 1) / complex<double>(e, 0));
    vector<Complex> y(e);

    vector<double> A0, A1;
    for (int i = 0; i < e; i++) {
        if (i % 2 == 0) {
            A0.push_back(a[i]);
        } else {
            A1.push_back(a[i]);
        }
    }

    vector<Complex> y0 = fft(A0, transformType);
    vector<Complex> y1 = fft(A1, transformType);

    Complex w = 1;
    for (int k = 0; k < e / 2; k++) {
        y[k] = y0[k] + w * y1[k];
        y[k + e / 2] = y0[k] - w * y1[k];
        w *= wn;
    }

    if (transformType == -1) { // inverse FFT
        for (int i = 0; i < e; i++) {
            y[i] /= e;
        }
    }
    return y;
}
