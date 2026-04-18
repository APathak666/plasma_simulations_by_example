#ifndef _CONFIG_H
#define _CONFIG_H

#include <string>

struct Config {
    int xcount, ycount, zcount;
    double dx, dy, dz;
    double origin_x, origin_y, origin_z;
    double ni;
    double ne;
    double dt;
    int num_ts;

    static Config load(const std::string &path);
};

#endif
