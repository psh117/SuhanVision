#include "suhanpoint.h"


bool Point::operator== (const Point& a) const
{
    if(a.x == x && a.y == y) return true;
    return false;
}
bool Point::operator!= (const Point& a) const
{
    if(a.x == x && a.y == y) return false;
    return true;
}
Point& Point::operator= (const Point& src)
{
    x = src.x;
    y = src.y;
    return *this;
}

Point Point::operator- (const Point& src)
{
    return Point(x - src.x, y - src.y);
}

Point Point::operator+ (const Point& src)
{
    return Point(x + src.x, y + src.y);
}

Point Point::operator* (const int src)
{
    return Point(x * src, y * src);
}

Point Point::operator* (const double src)
{
    return Point(x * src, y * src);
}
