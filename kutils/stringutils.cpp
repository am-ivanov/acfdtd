#include "kutils/stringutils.h"

#include <string.h>

#include "kutils/errors.h"

//std::string kutils::WHITESPACES = " \n\t\r";

//const std::string kutils::VALUES_TRUE[kutils::NUM_OF_VARIANTS_FOR_BOOL] = {"true", "on", "yes", "1"};
//const std::string kutils::VALUES_FALSE[kutils::NUM_OF_VARIANTS_FOR_BOOL] = {"false", "off", "no", "0"};

std::string kutils::replace(const std::string& srcp, const std::string& what, const std::string& with)
{
	std::string src = srcp;
	if (what != with) {
		std::string temp;
		std::string::size_type prev_pos = 0, pos = src.find(what, 0);
		while (std::string::npos != pos) {
			temp += std::string(src.begin() + prev_pos, src.begin() + pos) + with;
			prev_pos = pos + what.size();
			pos = src.find(what, prev_pos);
		}
		if (!temp.empty()) {
			src = temp + std::string(src.begin() + prev_pos, src.end());
			if (std::string::npos == with.find(what)) {
				kutils::replace(src, what, with);
			}
		}
	}
	return src;
}


std::string kutils::trim(const std::string & str)
{
	return kutils::trim(str, kutils::WHITESPACES);
}

std::string kutils::trim(const std::string & str, std::string whitespaces)
{
	std::string st = str;
	//удаляем что впереди
	std::string::size_type notwhite = st.find_first_not_of(whitespaces);
	st.erase(0, notwhite);

	//удаляем сзади
	notwhite = st.find_last_not_of(whitespaces);
	st.erase(notwhite + 1);

	return st;
}

void kutils::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
	//Пропускаем разделители в начале строки.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	//Первый не разделитель.
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos) {
		//Ищем вхождения и бобавляем их в список.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		//Пропускаем разделители.
		lastPos = str.find_first_not_of(delimiters, pos);
		//Ищем следующий символ, не разделитель.
		pos = str.find_first_of(delimiters, lastPos);
	}
}

bool kutils::toInt(const std::string& input, int& value)
{
	return kutils::stringTo<int>(input, value);
}


bool kutils::toFloat(const std::string& input, float& value)
{
	return kutils::stringTo<float>(input, value);
}

bool kutils::toDouble(const std::string& input, double& value)
{
	return kutils::stringTo<double>(input, value);
}


bool kutils::toUChar(const std::string& input, unsigned char& value)
{
	return kutils::stringTo<unsigned char>(input, value);
}

float kutils::toFloatDef(const std::string& input, const float defValue)
{
	float value;
	return kutils::toFloat(input, value) ? value : defValue;
}

int kutils::toIntDef(const std::string& input, const int defValue)
{
	int value;
	return kutils::toInt(input, value) ? value : defValue;
}

bool kutils::toBool(const std::string& input, bool& value)
{
	for (int i = 0; i < (int)kutils::NUM_OF_VARIANTS_FOR_BOOL; i++) {
		if (input == kutils::VALUES_TRUE[i]) {
			value = true;
			return true;
		}
		if (input == kutils::VALUES_FALSE[i]) {
			value = false;
			return true;
		}
	}
	return false;
}

bool kutils::stringTo(const std::string& from, bool& to)
{
	return kutils::toBool(from, to);
}

std::string kutils::fromBytes(const char* bytes, int len)
{
	char* str = new char[len + 1];
	memcpy(str, bytes, len);
	str[len] = '\0';
	return kutils::toString(str);
}


bool kutils::startWith(std::string str, std::string startwith)
{
	if (str.size() < startwith.size()) {
		return false;
	}
	int count = 0;
	for (int i = 0; i < (int)startwith.size(); i++) {
		if (str[i] == startwith[i]) {
			count++;
		}
	}
	return count == (int)startwith.size();
}

std::string kutils::tolower(const std::string& s)
{
	std::string d(s);
	std::transform(d.begin (), d.end (), d.begin (), (int(*)(int)) std::tolower);
	return d;
}

std::string kutils::toupper(const std::string& s)
{
	std::string d(s);
	std::transform(d.begin (), d.end (), d.begin (), (int(*)(int)) std::toupper);
	return d;
}

int kutils::toInt(const std::string &inp)
{
	int v;
	if (!kutils::toInt(inp, v)) {
		throw kutils::ParseError("kutils::toInt: Can't parse value \"" + inp + "\".");
	}
	return v;
}

float kutils::toFloat(const std::string &inp)
{
	float v;
	if (!kutils::toFloat(inp, v)) {
		throw kutils::ParseError("kutils::toFloat: Can't parse value \"" + inp + "\".");
	}
	return v;
}

double kutils::toDouble(const std::string &inp)
{
	double v;
	if (!kutils::toDouble(inp, v)) {
		throw kutils::ParseError("kutils::toDouble: Can't parse value \"" + inp + "\".");
	}
	return v;
}

bool kutils::toBool(const std::string &inp)
{
	bool v;
	if (!kutils::toBool(inp, v)) {
		throw kutils::ParseError("kutils::toBool: Can't parse value \"" + inp + "\".");
	}
	return v;
}
