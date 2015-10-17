#ifndef POINT_H
#define POINT_H

#include "math.h"

struct Point
{
    Point() : x(0), y(0) {}
    Point(int x, int y) : x(x), y(y) {}
    bool operator== (const Point& a) const;
    bool operator!= (const Point& a) const;
    Point& operator= (const Point& src);
    Point operator* (const int src);
    Point operator* (const double src);
    Point operator- (const Point& src);
    Point operator+ (const Point& src);
    double getDistance(Point a)
    {
        double dist;
        dist = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
        dist = sqrt(dist);
        return dist;

    }

    static double getDistance(Point a, Point b)
    {
        double dist;
        dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
        dist = sqrt(dist);
        return dist;
    }

    int x;
    int y;
};

#endif // POINT_H
