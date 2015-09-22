#ifndef SUHANMATH_H
#define SUHANMATH_H

#include <cstdarg>


namespace SuhanMath{


/**
 * @brief max
 * @param count Max값을 구할 가변 파라미터의 길이
 * @return 최대값을 반환
 * @bug int, long, double, char만 사용 가능 unsigned char, float 등 사용시 오류 발생(Qt5.2 기준)
 */
template <class T>
T max(int count, ...)
{
    T max = 0;
    int i;

    va_list ap;
    va_start(ap, count);

    for(i=0; i<count; i++)
    {
        T argValue = va_arg(ap, T);
        if(argValue > max)
        {
            max = argValue;
        }
    }
    va_end(ap);
    return max;
}

/**
 * @brief min
 * @param count Max값을 구할 가변 파라미터의 길이
 * @return 최대값을 반환
 * @bug int, long, double, char만 사용 가능 unsigned char, float 등 사용시 오류 발생(Qt5.2 기준)
 */
template <class T>
T min(int count, ...)
{
    T min;
    int i;

    va_list ap;
    va_start(ap, count);

    min = va_arg(ap, T);
    for(i=1; i<count; i++)
    {
        T argValue = va_arg(ap, T);
        if(argValue < min)
        {
            min = argValue;
        }
    }
    va_end(ap);
    return min;
}

}



#endif // SUHANMATH_H
