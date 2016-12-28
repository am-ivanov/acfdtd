#ifndef FD_COEFF_H
#define FD_COEFF_H

#include "common.h"
#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include <vector>

inline double fact(const unsigned int n) {
	double res = 1;
	for (unsigned int i = 2; i <= n; ++i) res *= i;
	return res;
}

inline double partFact(const unsigned int start, const unsigned int end) {
	double res = 1;
	if (start+1 > end) throw std::logic_error("start+1 > end !!!");
	for (unsigned int i = start+1; i <= end; ++i) res *= i;
	return res;
}

class FDCoeff {
public:
	/**
	 * \brief Calculate FD coefficients
	 * \param N half of order
	 */
	void calc(int_t N) {
		// first derivative
		m_c1.resize(N+1, 0);
		m_c1.at(0) = 0;
		for (int_t n = 1; n <= N; ++n) {
			double mul = 1.0;
			for (int_t n_mul = 1; n_mul <= n; ++n_mul) {
				mul *= static_cast<double>(N - n + n_mul) / static_cast<double>(N + n_mul);
			}
			m_c1.at(n) = (n % 2 ? 1.0 : -1.0) * mul / n;
		}
		// second derivative
		m_c2.resize(N+1, 0);
		for (int_t n = 1; n <= N; ++n) {
			m_c2.at(n) = m_c1.at(n) * 2.0 / n;
			m_c2.at(0) += m_c2.at(n);
		}
		m_c2.at(0) *= -2;
		// first derivative staggered
		m_sc1.resize(N, 0);
		if (N == 1) m_sc1.at(0) = 1;
		else {
			for (int_t n = 0; n < N; ++n) {
				m_sc1.at(n) = pow(2,5) * pow(-1, n) * (2*N-1);
				m_sc1.at(n) /= pow(2*n+1, 2) * pow(2, 4*N);
				m_sc1.at(n) *= partFact(N-1, 2*N-1) * partFact(N-2, 2*N-3);
				m_sc1.at(n) /= fact(N-n-1) * fact(N+n);
			}
		}
	}
	/**
	 * Get specified FD coefficient
	 * \param n number of coefficient
	 */
	real_t c1(int_t n) {
		if (m_c1.empty())
			throw std::logic_error("Call calc() first");
		if (n >= 0) {
			return m_c1.at(n);
		} else {
			return - m_c1.at(-n);
		}
	}
	real_t c2(int_t n) {
		if (m_c1.empty())
			throw std::logic_error("Call calc() first");
		return m_c2.at(labs(n));
	}
	real_t sc1(const int_t n) const {
		#ifdef DEBUG
			if (m_c1.empty())
				throw std::logic_error("Call calc() first");
			}
			if (!(n >= 1)) {
				throw std::logic_error("Can't get staggered zero coefficient");
		#endif
		return m_sc1[n-1];
	}
private:
	std::vector<real_t> m_c1; // first derivatice coefficients
	std::vector<real_t> m_c2; // second derivative coefficients
	std::vector<real_t> m_sc1; // first derivatice staggered coefficients
};

#endif
