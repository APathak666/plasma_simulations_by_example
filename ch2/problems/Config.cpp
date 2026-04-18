#include "Config.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <set>

static std::string trim(const std::string &s) {
    auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

Config Config::load(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open config file: " + path);

    Config cfg{};
    std::set<std::string> seen;
    std::string line;
    int lineno = 0;

    while (std::getline(in, line)) {
        lineno++;
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        line = trim(line);
        if (line.empty()) continue;

        auto eq = line.find('=');
        if (eq == std::string::npos)
            throw std::runtime_error("Config line " + std::to_string(lineno) + ": missing '='");

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));

        if      (key == "xcount")   cfg.xcount   = std::stoi(val);
        else if (key == "ycount")   cfg.ycount   = std::stoi(val);
        else if (key == "zcount")   cfg.zcount   = std::stoi(val);
        else if (key == "dx")       cfg.dx       = std::stod(val);
        else if (key == "dy")       cfg.dy       = std::stod(val);
        else if (key == "dz")       cfg.dz       = std::stod(val);
        else if (key == "origin_x") cfg.origin_x = std::stod(val);
        else if (key == "origin_y") cfg.origin_y = std::stod(val);
        else if (key == "origin_z") cfg.origin_z = std::stod(val);
        else if (key == "ni")       cfg.ni       = std::stod(val);
        else if (key == "ne")       cfg.ne       = std::stod(val);
        else if (key == "dt")       cfg.dt       = std::stod(val);
        else if (key == "num_ts")   cfg.num_ts   = std::stoi(val);
        else throw std::runtime_error("Config line " + std::to_string(lineno) + ": unknown key '" + key + "'");

        seen.insert(key);
    }

    const std::set<std::string> required = {
        "xcount", "ycount", "zcount",
        "dx", "dy", "dz",
        "origin_x", "origin_y", "origin_z",
        "ni", "ne",
        "dt", "num_ts"
    };
    for (const auto &k : required)
        if (seen.find(k) == seen.end())
            throw std::runtime_error("Config: missing required key '" + k + "'");

    return cfg;
}
