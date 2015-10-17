#ifndef SUHANMATH_H
#define SUHANMATH_H

#include <cstdarg>
#include <cstdlib>

namespace SuhanMath{

const double PI = 3.1415926535897932384626433832795;

/**
 * @brief max
 * @param count Max값을 구할 가변 파라미터의 길이
 * @return 최대값을 반환
 * @bug int, long, double, char만 사용 가능 unsigned char, float 등 사용시 오류 발생(Qt5.2 기준)
 */
template <class T>
T Max(int count, ...)
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
T Min(int count, ...)
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


template <class T>
T Square(T value)
{
    return (value * value);
}



template <class T>
class SHArray
{
public:
    SHArray() : pData(NULL), nDataLength(0) {}
    SHArray(const SHArray& ref) : pData(NULL)
    { *this = ref; }

    SHArray(int nSize) : pData(NULL), nDataLength(nSize)
    { Create(nSize); }

    virtual ~SHArray()
    {
        if(pData != NULL)
            delete pData;
    }
    void Create(int nSize)
    {
        if(pData != NULL)
            delete pData;
        nDataLength = nSize;
        pData = new T[nSize];
        for(int i=0; i<nSize;i++)
        {
            pData[i] = 0;
        }
    }
    void Reset() { if(pData != NULL) delete pData; nDataLength = 0; }

    // getter
    int Length() const { return nDataLength; }

    // opeartor overloading
    SHArray&    operator  =(const SHArray& ref)
    {
        Create(ref.Length());
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] = ref[i];
        }
    }


    SHArray&    operator +=(const SHArray& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] += ref[i];
        }
        return *this;
    }

    SHArray&    operator -=(const SHArray& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] -= ref[i];
        }
        return *this;
    }
    SHArray&    operator *=(const SHArray& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] *= ref[i];
        }
        return *this;
    }

    SHArray&    operator /=(const SHArray& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] /= ref[i];
        }
        return *this;
    }

    SHArray&    operator +=(const T& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] += ref;
        }
        return *this;
    }
    SHArray&    operator -=(const T& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] -= ref;
        }
        return *this;
    }
    SHArray&    operator *=(const T& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] *= ref;
        }
        return *this;
    }
    SHArray&    operator /=(const T& ref)
    {
        for(int i=0; i<nDataLength; i++)
        {
            pData[i] /= ref;
        }
        return *this;
    }

    SHArray     operator  +(const SHArray& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] += ref[i];
        }
        return dst;
    }
    SHArray     operator  +(const T& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] += ref;
        }
        return dst;
    }

    SHArray     operator  -(const SHArray& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] -= ref[i];
        }
        return dst;
    }
    SHArray     operator  -(const T& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] -= ref;
        }
        return dst;
    }
    SHArray     operator  *(const T& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] *= ref;
        }
        return dst;
    }

    SHArray     operator  /(const T& ref) const
    {
        SHArray dst(*this);
        for(int i=0; i<nDataLength; i++)
        {
            dst[i] /= ref;
        }
        return dst;
    }
    bool operator ==(const SHArray& ref) const
    {
        for(int i=0; i<nDataLength; i++)
        {
            if(pData[i] != ref[i]) return false;
        }
        return true;
    }

    T& operator [] (const int& i) { return pData[i]; }
    T& operator [] (const int& i) const { return pData[i]; }

protected:
    T* pData;
    int nDataLength;
};

template <class T>
class SHMatrix : public SHArray<T>
{
public:
    SHMatrix(int row, int col) : SHArray<T>(row*col), nRow(row), nCol(col) {}
    T& at(int row, int col) { return this->pData[row * nCol + col]; }
    virtual ~SHMatrix() {}

private:
    int nRow;
    int nCol;
};

}

#endif // SUHANMATH_H
