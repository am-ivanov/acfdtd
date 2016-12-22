#ifndef STRING_H
#define STRING_H

#include "kutils/types.h"
#include "kutils/errors.h"

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#define STR(x) kutils::toString(x)

namespace kutils {
	std::string replace(const std::string& srcp, const std::string& what, const std::string& with);
	std::string trim(const std::string& str);
	std::string trim(const std::string& str, std::string whitespaces);
	void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);

	bool toInt(const std::string& input, int& value);
	bool toUChar(const std::string& input, unsigned char& value);
	bool toFloat(const std::string& input, float& value);
	bool toDouble(const std::string& input, double& value);
	bool toBool(const std::string& input, bool& value);

	int toInt(const std::string &inp);
	float toFloat(const std::string &inp);
	double toDouble(const std::string &inp);
	bool toBool(const std::string &inp);

	float toFloatDef(const std::string& input, const float defValue);
	int toIntDef(const std::string& input, const int defValue);

	std::string fromBytes(const char* bytes, int len);

	bool startWith(std::string str, std::string startwith);

	template <typename T> bool stringTo(const std::string& from, T& to)
	{
		std::istringstream iss(from);
		iss >> to;
		return true;
	}

	bool stringTo(const std::string& from, bool& to);


	template <typename T> std::string toString(const T& x)
	{
		std::ostringstream oss;
		oss << x;
		return oss.str();
	}

	template <typename T>
	bool toArray(const std::string &inp, std::vector<T> &v, const std::string &delim = ",")
	{
		std::vector<std::string> tok;
		tokenize(inp, tok, delim);
		for (int i = 0; i < (int)tok.size(); i++) {
			T t;
			if (!stringTo<T>(trim(tok[i]), t)) {
				return false;
			}
			v.push_back(t);
		}
		return true;
	}

	template <typename T>
	T toT(const std::string &inp)
	{
		T v;
		if (!stringTo(inp, v)) {
			throw ParseError("kutils::toT: Can't parse value \"" + inp + "\".");
		}
		return v;
	}


	std::string tolower(const std::string& s);
	std::string toupper(const std::string& s);

	static std::string WHITESPACES = " \n\t\r";

	static const int NUM_OF_VARIANTS_FOR_BOOL = 4;
	static const std::string VALUES_TRUE[NUM_OF_VARIANTS_FOR_BOOL]  = {"true", "on", "yes", "1"};
	static const std::string VALUES_FALSE[NUM_OF_VARIANTS_FOR_BOOL] = {"false", "off", "no", "0"};

}

#endif // STRING_H
