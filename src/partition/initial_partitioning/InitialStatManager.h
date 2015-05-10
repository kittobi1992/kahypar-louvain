/*
 * InitialStatManager.h
 *
 *  Created on: May 7, 2015
 *      Author: theuer
 */

#ifndef SRC_PARTITION_INITIAL_PARTITIONING_INITIALSTATMANAGER_H_
#define SRC_PARTITION_INITIAL_PARTITIONING_INITIALSTATMANAGER_H_

#define RESET   "\033[0m"
#define BOLD   "\033[1m"
#define BLUE    "\033[34m"

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

namespace ip = boost::interprocess;

namespace partition {

class InitialStatManager {

public:
	static InitialStatManager& getInstance() {
		static InitialStatManager stat;
		return stat;
	}

	void addStat(std::string group, std::string name, double stat,
			std::string style = "", bool result = false) {
		stats[group].push_back(std::make_pair(name, stat));
		styles[group].push_back(std::make_pair(name, style));
		results[group].push_back(std::make_pair(name, result));
	}

	double getStat(std::string group, std::string name) {
		for (auto v : stats[group]) {
			if (v.first.compare(name) == 0) {
				return v.second;
			}
		}
		return 0.0;
	}

	void updateStat(std::string group, std::string name, double value) {
		for (int i = 0; i < stats[group].size(); i++) {
			if (stats[group][i].first.compare(name) == 0) {
				stats[group][i] = std::make_pair(name, value);
				return;
			}
		}
		addStat(group, name, value);
	}

	void printStats(std::string filename, bool enable_style) {
		std::ostringstream oss;
		for (auto s : stats) {

			if (enable_style) {
				oss << BOLD << BLUE;
			}
			for (int i = 0; i < s.first.size() + 18; i++)
				oss << "=";
			oss << "\n======== " << s.first << " ========\n";
			for (int i = 0; i < s.first.size() + 18; i++)
				oss << "=";
			if (enable_style) {
				oss << RESET;
			}
			oss << "\n";

			for (int i = 0; i < s.second.size(); i++) {
				if (enable_style) {
					oss << styles[s.first][i].second << s.second[i].first
							<< " = " << s.second[i].second << RESET << "\n";
				}
				else {
					oss <<  s.second[i].first
							<< " = " << s.second[i].second << "\n";
				}
			}
		}
		if (!filename.empty()) {
			std::ofstream out_stream(filename.c_str(), std::ofstream::app);
			ip::file_lock f_lock(filename.c_str());
			{
				ip::scoped_lock<ip::file_lock> s_lock(f_lock);
				out_stream << oss.str();
				out_stream.flush();
			}
		} else {
			std::cout << oss.str() << std::endl;
		}
	}

	void setGraphName(std::string graph) {
		graphname = graph;
	}
	void setMode(std::string m) {
		mode = m;
	}

	void printResultLine(std::string filename) {
		std::ostringstream oss;
		oss << "RESULT graph=" << graphname << " mode=" << mode;
		for (auto s : stats) {

			for (int i = 0; i < s.second.size(); i++) {
				if (results[s.first][i].second) {
					oss << " " << s.second[i].first << "="
							<< s.second[i].second;
				}
			}
		}
		oss << "\n";
		if (!filename.empty()) {
			std::ofstream out_stream(filename.c_str(), std::ofstream::app);
			ip::file_lock f_lock(filename.c_str());
			{
				ip::scoped_lock<ip::file_lock> s_lock(f_lock);
				out_stream << oss.str();
				out_stream.flush();
			}
		} else {
			std::cout << oss.str() << std::endl;
		}
	}

	std::string graphname;
	std::string mode;

	std::map<std::string, std::vector<std::pair<std::string, double>>>stats;
	std::map<std::string, std::vector<std::pair<std::string, std::string>>> styles;
	std::map<std::string, std::vector<std::pair<std::string, bool>>> results;

private:
	InitialStatManager() {};

	~InitialStatManager() {};

};

}

#endif /* SRC_PARTITION_INITIAL_PARTITIONING_INITIALSTATMANAGER_H_ */