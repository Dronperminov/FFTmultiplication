#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cmath>

typedef std::complex<double> base;

void fft(std::vector<base> &a, bool invert = false) {
	size_t n = a.size();
 
	for (size_t i = 1, j = 0; i < n; i++) {
		size_t bit = n >> 1;

		while (j >= bit) {
			j -= bit;
			bit >>= 1;
		}

		j += bit;

		if (i < j)
			swap(a[i], a[j]);
	}
 
	for (size_t len = 1; len < n; len <<= 1) {
		double ang = M_PI / len * (invert ? -1 : 1);
		base wlen (cos(ang), sin(ang));

		size_t len2 = len << 1;

		for (size_t i = 0; i < n; i += len2) {
			base w(1);

			for (size_t j = 0; j < len; j++) {
				base u = a[i + j];
				base v = a[i + j + len] * w;

				a[i + j] += v;
				a[i + j + len] = u - v;
				w *= wlen;
			}
		}
	}

	if (invert)
		for (size_t i = 0; i < n; i++)
			a[i] /= n;
}


void multiply (const std::vector<int> &a, const std::vector<int> &b, std::vector<int> &res) {
	std::vector<base> fa (a.begin(), a.end());
	std::vector<base> fb (b.begin(), b.end());

	size_t n = 4;
	size_t maxSize = a.size() + b.size() + 1;

	while (n < maxSize)
		n <<= 1;

	fa.resize(n);
	fb.resize(n);
 
	fft(fa);
	fft(fb);

	for (size_t i = 0; i < n; i++)
		fa[i] *= fb[i];

	fft(fa, true);
 
	res.resize(maxSize);
	size_t carry = 0;

	for (size_t i = 0; i < maxSize; i++) {
		res[i] = fa[i].real() + 0.5 + carry;
		carry = res[i] / 10;
		res[i] %= 10;
	}
}

int main() {
	std::string s;

	while(getline(std::cin, s)) {
		std::vector<int> a;
		std::vector<int> b;
		std::vector<int> c;

		size_t i = s.length() - 1;
		while (s[i] != ' ')
			a.push_back(s[i--] - '0');

		while (i > 0)
			b.push_back(s[--i] - '0');

		multiply(a, b, c);

		for (size_t i = 0; i < c.size(); i++)
			std::cout << c[c.size() - 1 - i];

		std::cout << std::endl;
	}
}