#ifndef CFD24_MACROS_HPP
#define CFD24_MACROS_HPP

#include <stdexcept>
#include <sstream>

#define _THROW_NOT_IMP_ \
{\
	std::ostringstream oss; \
	oss << "NOT IMPLEMENTED method called." << std::endl; \
	oss << __PRETTY_FUNCTION__ << std::endl; \
	oss << "At" << std::endl; \
	oss << __FILE__ << ": " << __LINE__ << std::endl; \
	throw std::runtime_error(oss.str()); \
}

#endif
