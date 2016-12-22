#include "kutils/ktimer.h"

#include <map>
#include <iostream>

using namespace std;
using namespace kutils;

TimeProfiler::TimeProfiler(const std::string &name)
{
	time = timer();
	this->name = name;
}

TimeProfiler::~TimeProfiler()
{
	ProfileInfo pi;
	pi.name = name;
	pi.time = timer() - time;
	data().push_back(pi);
}

void TimeProfiler::report()
{
	typedef map<string, pair<double, int> > MapType;
	MapType m;
	for (ListInfo::iterator it = data().begin(); it != data().end(); it++) {
		if (m.count((*it).name) > 0) {
			m[(*it).name].first += (*it).time;
			m[(*it).name].second++;
		} else {
			m[(*it).name] = make_pair((*it).time, 1);
		}
	}
	cout << "Name\t\tTrials\tTime, s\t\tAvg time, s" << endl;
	for (MapType::iterator it = m.begin(); it != m.end(); it++) {
		cout << (*it).first << "\t";
		cout << (*it).second.second << "\t";
		cout << (*it).second.first << "\t";
		cout << (*it).second.first / (*it).second.second << endl;
	}
}

TimeProfiler::ListInfo& TimeProfiler::data()
{
	static ListInfo data_;
	return data_;
}
