
#include <fstream>
#include <sstream>

#include "kutils/config.h"
#include "kutils/stringutils.h"

using namespace kutils;

std::string ConfigEntry::DELIMETER = "/";

std::string Config::COMMENTS[Config::COMMENTS_NUM] = {"#", "//", ";"};
std::string Config::QUOTES[Config::QUOTES_NUM] = {"\"", "'"};
std::string Config::LEFT_BRACKET = "[";
std::string Config::RIGHT_BRACKET = "]";
std::string Config::SLASH = "/";
std::string Config::EQUALS = "=";

ConfigEntry::ConfigEntry()
{
	name = "";
	type = VALUE;
	value = "";
}


ConfigEntry::ConfigEntry(std::string name, int type, std::string value) throw(ConfigEntry::InvalidType)
{
	this->name = name;
	setType(type);
	if (type == VALUE) {
		this->value = value;
	}
}


ConfigEntry::~ConfigEntry()
{
	entries.clear();
}


void ConfigEntry::setType(int type) throw(ConfigEntry::InvalidType)
{
	if (type != SECTION &&
		type != VALUE ) {
		throw InvalidType("Unknown entry type: " + STR(type));
	}
	this->type = type;
}


std::string ConfigEntry::getValue() throw(ConfigEntry::InvalidEntry)
{
	if (type != VALUE) {
		throw InvalidEntry("This operation allowed only for VALUE type entries.");
	}
	return value;
}


void ConfigEntry::setValue(std::string value) throw(ConfigEntry::InvalidEntry)
{
	if (type != VALUE) {
		throw InvalidEntry("This operation allowed only for VALUE type entries.");
	}
	this->value = value;
}


std::vector<ConfigEntry> ConfigEntry::getEntries(std::vector<ConfigEntry>& entrs) throw(ConfigEntry::InvalidEntry)
{
	if (type != SECTION) {
		throw InvalidEntry("This operation allowed only for SECTION type entries.");
	}
	MapEntries::iterator it;
	for (it = entries.begin(); it != entries.end(); ++it) {
		entrs.push_back((*it).second);
	}
	return entrs;
}

ConfigEntry& ConfigEntry::find(std::string key) throw(ConfigEntry::InvalidEntry, NotFound)
{
	if (type != SECTION) {
		throw InvalidEntry("This operation allowed only for SECTION type entries.");
	}
	if (entries.count(key) < 1) {
		//STACK_TRACE;
		throw NotFound("Entry with key \"" + key + "\" not found.");
	}

	return (*entries.find(key)).second;
}

bool ConfigEntry::hasEntry(const std::string &path)
{
	bool ret = false;
	try {
		getEntry(path);
		ret = true;
	} catch(NotFound &e) {
	}
	return ret;
}

///todo: пока работает только над первыми вхождениями...
ConfigEntry& ConfigEntry::getEntry(std::string path) throw(InvalidEntry, NotFound)
{
	if (type != SECTION) {
		throw InvalidEntry("This operation allowed only for SECTION entries.");
	}

	//удаляем все с концов
	std::string newPath = trim(path, DELIMETER + WHITESPACES);

	//Теперь разбираем строку по вхождениям DELIMETER
	std::vector<std::string> tokens;
	tokenize(newPath, tokens, DELIMETER);
	//На всякий случай удаляем все пробельные символы во всех токенах (вообще их там быть не должно, но малоли кто будет использовать этот класс)
	for (int i = 0; i < (int)tokens.size(); i++) {
		tokens[i] = trim(tokens[i]);
	}

	if (tokens.size() == 0) {
		throw NotFound("Entry with path \"" + path + "\" not found.");
	} else if (tokens.size() == 1) {
		if (entries.count(tokens[0]) < 1) {
			throw NotFound("Entry with key \"" + path + "\" not found.");
		}
		return (*entries.find(tokens[0])).second;
	} else {
		if (entries.count(tokens[0]) < 1) {
			throw NotFound("Entry with key \"" + path + "\" not found.");
		}
		MapEntries::iterator it = entries.find(tokens[0]);
		if ((*it).second.getType() != SECTION) {
			throw InvalidEntry("This operation allowed only for SECTION entries.");
		}
		std::string nm = "";
		for (int j = 1; j < (int)tokens.size(); j++) {
			nm += tokens[j];
			nm += DELIMETER;
		}
		return (*it).second.getEntry(nm);
	}
}


void ConfigEntry::addEntry(ConfigEntry en) throw(ConfigEntry::InvalidEntry)
{
	if (type != SECTION) {
		throw InvalidEntry("This operation allowed only for SECTION entries.");
	}
	entries.insert(MapValue(en.getName(), en));
}



std::string kutils::ConfigEntry::toString(std::string space)
{
	std::ostringstream oss;
	if (type == VALUE) {
		oss << space << name << " = " << value;
	} else {
		oss << space << "[" << name << "]" << std::endl;
		MapEntries::iterator it;
		for (it = entries.begin(); it != entries.end(); ++it) {
			oss << (*it).second.toString("\t" + space) << std::endl;
		}
		oss << space << "[/" << name << "]";
	}
	return oss.str();
}

void ConfigEntry::getDoubleArray(std::vector<double> &a)
{
	if (!toArray<double>(getValue(), a)) {
		throw ParseError("ConfigEntry::getDoubleArray: Can't parse array: " + getValue());
	}
}

void ConfigEntry::getIntArray(std::vector<int> &a)
{
	if (!toArray<int>(getValue(), a)) {
		throw ParseError("ConfigEntry::getIntArray: Can't parse array: " + getValue());
	}
}

void ConfigEntry::getBoolArray(std::vector<bool> &a)
{
	if (!toArray<bool>(getValue(), a)) {
		throw ParseError("ConfigEntry::getBoolArray: Can't parse array: " + getValue());
	}
}

void ConfigEntry::getStringArray(std::vector<std::string> &a)
{
	if (!toArray<std::string>(getValue(), a)) {
		throw ParseError("ConfigEntry::getStringArray: Can't parse array: " + getValue());
	}
}




int kutils::Config::getLineType(std::string line)
{
	//если нулевой длинны - пропускаем
	if (line.size() == 0) {
		return LINE_COMMENT;
	}

	//проверяем на секции
	if (line.substr(0, 1) == LEFT_BRACKET && line.substr(line.length() - 1, 1) == RIGHT_BRACKET) {
		if (line.substr(1, 1) == SLASH) {
			return LINE_SECTION_END;
		} else {
			return LINE_SECTION_BEGIN;
		}
	}

	//ищем каменты
	for (int i = 0; i < COMMENTS_NUM; i++) {
		if (kutils::startWith(line, COMMENTS[i])) {
			return LINE_COMMENT;
		}
	}

	//значения
	if (line.find(EQUALS) > 0 && line.find(EQUALS) < line.length()) {
		return LINE_VALUE;
	}

	//не нашли - возвращаем не найдено
	return LINE_UNKNOWN;
}


std::string Config::getSectionName(std::string line, int type)
{
	std::string sec = line;
	sec.erase(0, 1);
	sec.erase(sec.length() - 1, 1);
	if (type == LINE_SECTION_END) {
		sec.erase(0, 1);
	}
	return sec;
}


std::pair<std::string, std::string> Config::getNameValue(std::string line)
{
	std::pair<std::string, std::string> ret;
	ret.first = trim(line.substr(0, line.find(EQUALS)));
	ret.second = trim(line.substr(line.find(EQUALS) + 1, line.length()));
	return ret;
}


void kutils::Config::parseSection(ConfigEntry& entry, std::istream& stream) throw(kutils::ParseError)
{
	if (!stream.good()) {
		throw kutils::ParseError("Bad stream.");
	}

	if (entry.getType() != ConfigEntry::SECTION) {
		throw ConfigEntry::InvalidEntry("This operation allowed only for SECTION entris.");
	}

	std::string line;

	while (!stream.eof()) {
		std::getline(stream, line);
		line = trim(line);
		//stream >> line;

		switch(getLineType(line)) {
			case LINE_COMMENT:
				//nothing to do
				break;
			case LINE_SECTION_BEGIN:
				{
					ConfigEntry ent;
					ent.setType(ConfigEntry::SECTION);
					ent.setName(getSectionName(line, LINE_SECTION_BEGIN));
					parseSection(ent, stream);
					entry.addEntry(ent);
				}
				break;
			case LINE_SECTION_END:
				if (getSectionName(line, LINE_SECTION_END) != entry.getName()) {
					throw ParseError("Invalid section end: \"" + line + "\", current section: \"" + entry.getName() + "\".");
				}
				return;
				break;
			case LINE_VALUE:
				entry.addEntry(ConfigEntry(getNameValue(line).first, ConfigEntry::VALUE, getNameValue(line).second));
				break;
			case LINE_UNKNOWN:
				throw ParseError("Some errors near line: \"" + line + "\"");
				break;
		}
	}
}


void kutils::Config::parseStream(std::istream& stream) throw(ParseError)
{
	entries.setName("");
	entries.setType(ConfigEntry::SECTION);

	parseSection(entries, stream);
}


void Config::parseFile(std::string path) throw(ParseError, OpenFailed)
{
	std::ifstream ifs;
	ifs.open(path.c_str());
	if (!ifs.good()) {
		throw OpenFailed("Error open file: \"" + path + "\".");
	}

	parseStream(ifs);

	ifs.close();
}


void Config::parseString(std::string str) throw(ParseError)
{
	std::istringstream iss(str);
	parseStream(iss);
}


