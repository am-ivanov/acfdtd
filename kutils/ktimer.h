#ifndef KTIMER_H
#define KTIMER_H

#include <sys/time.h>

#include <string>
#include <list>

#define PROFILE(NAME) kutils::TimeProfiler timeprofiler(NAME);
#define PROFILE_REPORT() kutils::TimeProfiler::report();

namespace kutils {
	inline double timer()
	{
		timeval theStartTime;
		gettimeofday(&theStartTime, NULL);
		return theStartTime.tv_sec + 1e-6 * theStartTime.tv_usec;
	}
	
	class TimeProfiler {
		public:
			TimeProfiler(const std::string &name);
			~TimeProfiler();
			static void report();
		private:
			struct ProfileInfo {
				std::string name;
				double time;
			};
			typedef std::list<ProfileInfo> ListInfo;
			double time;
			std::string name;
			
			static ListInfo& data();
	};
}

#endif // KTIMER_H
