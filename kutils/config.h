/*
 * Copyright (C) 2009 Nikolay Khokhlov <k_h@inbox.ru>
 */

#ifndef KUTILS_CONFIG_H
#define KUTILS_CONFIG_H

#include "kutils/errors.h"
#include "kutils/stringutils.h"

#include <map>
#include <vector>


namespace kutils {
	class ConfigEntry {
		public:
			class InvalidType : public std::logic_error {
				public:
					InvalidType(const std::string& whatArg) :
						logic_error(whatArg) { }
			};
			
			class InvalidEntry : public std::logic_error {
				public:
					InvalidEntry(const std::string& whatArg) :
						logic_error(whatArg) { }
			};
			
			enum TYPE {
				SECTION = 0,
				VALUE
			};
			
			ConfigEntry();
			ConfigEntry(std::string name, int type, std::string value) throw(InvalidType);
			~ConfigEntry();
			
			void setName(std::string name) { this->name = name; }
			std::string getName() { return name; }
			
			void setType(int type) throw(InvalidType);
			int getType() { return type; }
			
			void setValue(std::string value) throw(InvalidEntry);
			std::string getValue() throw(InvalidEntry);
			
			std::vector<ConfigEntry> getEntries(std::vector<ConfigEntry>& entrs) throw(InvalidEntry);
			
			ConfigEntry& getEntry(std::string path) throw(InvalidEntry, NotFound);
			
			bool hasKey(const std::string &key) { return entries.count(key) > 0; }
			bool hasEntry(const std::string &path);
			
			void addEntry(ConfigEntry en) throw(InvalidEntry);
			
			std::string toString(std::string space = "");
			
			template <typename T>
			T get() { return kutils::toT<T>(getValue()); }
			template <typename T>
			T get(const char *path) { return getEntry(path).get<T>(); }
			
			template <typename T>
			T get(const char *path, const T& def) {
				if (hasEntry(path)) {
					return get<T>(path);
				}
				return def;
			}
			
			
			void getDoubleArray(std::vector<double> &a);
			void getIntArray(std::vector<int> &a);
			void getStringArray(std::vector<std::string> &a);
			void getBoolArray(std::vector<bool> &a);
			
			
		
		private:
			typedef std::pair<std::string, ConfigEntry> MapValue;
			typedef std::multimap<std::string, ConfigEntry> MapEntries;
			
			ConfigEntry& find(std::string key) throw(InvalidEntry, NotFound);
			
			std::string name;
			std::string value;
			int type;
			MapEntries entries;
			
			static std::string DELIMETER;
	};

	
	class Config {
		public:
			Config() {};
			~Config() {};
			void parseFile(std::string path) throw(OpenFailed, ParseError);
			void parseString(std::string str) throw(ParseError);
			void parseStream(std::istream& stream) throw(ParseError);
			ConfigEntry& getRootEntry() { return entries; }
			ConfigEntry& getEntry(std::string path) throw(ConfigEntry::InvalidEntry, NotFound) { return entries.getEntry(path); }
			std::string getValue(std::string path) throw(ConfigEntry::InvalidEntry, NotFound) { return entries.getEntry(path).getValue(); }
			template <typename T>
			T get(const std::string &path) { return entries.getEntry(path).get<T>(); }
			bool hasKey(const std::string &key) { return entries.hasKey(key); }
		private:
			int getLineType(std::string line);
			void parseSection(ConfigEntry& ent, std::istream& stream) throw(ParseError);
			std::string getSectionName(std::string line, int type);
			std::pair<std::string, std::string> getNameValue(std::string line);
			
			enum {
				LINE_COMMENT = 0,
				LINE_VALUE,
				LINE_SECTION_BEGIN,
				LINE_SECTION_END,
				LINE_UNKNOWN
			};
			
			static const int COMMENTS_NUM = 3;
			static std::string COMMENTS[COMMENTS_NUM];
			static const int QUOTES_NUM = 2;
			static std::string QUOTES[QUOTES_NUM];
			static std::string LEFT_BRACKET;
			static std::string RIGHT_BRACKET;
			static std::string SLASH;
			static std::string EQUALS;
			
			
			ConfigEntry entries;
	};
};

#endif // CONFIG_H
