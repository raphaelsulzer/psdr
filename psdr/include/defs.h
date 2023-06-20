#pragma once

template <typename T>
bool jin(const T & a, const T & x, const T & b) {
	return (a <= x && x <= b);
}

template <typename T>
bool jinr(const T & a, const T & x, const T & b) {
	return (a < x && x < b);
}

template <typename T>
const T & jmin(const T & a, const T & b) {
	return (a < b ? a : b);
}

template <typename T>
const T & jmax(const T & a, const T & b) {
	return (a > b ? a : b);
}

template <typename T>
const T & jclamp(const T & a, const T & x, const T & b) {
	return (x < a ? a : (x > b ? b : x));
}

