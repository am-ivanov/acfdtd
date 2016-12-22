#ifndef ERRORS_H
#define ERRORS_H

#include "kutils/types.h"

#include <string>
#include <stdexcept>


namespace kutils {
	class OpenFailed : public std::runtime_error {
		public:
			OpenFailed(const std::string& whatArg) :
				runtime_error(whatArg) {}
	};
			
	class ParseError : public std::runtime_error {
		public:
			ParseError(const std::string& whatArg) :
				runtime_error(whatArg) {}
	};
	
	class NotFound : public std::runtime_error {
		public:
			NotFound(const std::string& whatArg) :
				runtime_error(whatArg) {}
	};
};

#endif // ERRORS_H
