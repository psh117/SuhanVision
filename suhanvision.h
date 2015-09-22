/**
  @file suhanvision.h
  @author 박수한 ParkSuhan (psh117@gmail.com)
  @brief 수한 영상처리 라이브러리 Suhan Vision Process Library

  @bug 아직 없음 not now..
  @todo 만드는 중 making...
  */


#ifndef SUHANVISION_H
#define SUHANVISION_H



#include "kfc.h"
#include "suhanmath.h"

typedef unsigned char byte;

/**
 * @brief The IMG_TYPE enum
 * 이미지 형식 지정
 */
enum IMG_TYPE {BINARY,GRAY,HSV,RGB,BGR};


using SuhanMath::max;
using SuhanMath::min;

/**
 * @brief The SHImage class
 * 수한 이미지 클래스 Suhan Image Class
 */
class SHImage
{
public:
    SHImage();
    SHImage(int col, int row, int channel = 1, IMG_TYPE type = GRAY);
    SHImage(KImageColor& refImage);
    virtual ~SHImage();

    void Create(int col, int row, int channel = 1, IMG_TYPE type = GRAY);
    ///< 이미지 데이터 생성 함수 (pData의 메모리를 할당해 줌) This function is used for allocating "pData"
    void Reset();   ///< 현재 가지고 있는 데이터 초기화 함수 This function is used for clearing all data of this class.

    void ToKImageColor(KImageColor& dstImage);

    byte& operator[] (const int& i) { return pData[i]; }
    byte& pixel(int col, int row, int channel = 0) { return pData[ ( col + (row * nCol) ) * nChannel + channel]; }

    void FillData(int channel, byte data);

    void ConvertToHSV(SHImage& dstImage);
    void ConvertToRGB(SHImage& dstImage);
    void ConvertToBGR(SHImage& dstImage);

private:
    IMG_TYPE nType; ///< 현재 이미지 타입 Current image type
    int nCol;       ///< 칸 수 (너비) a number of cols of image
    int nRow;       ///< 줄 수 (높이) a number of rows of image
    int nChannel;   ///< 채널 개수 a number of channel (1ch, 2ch, 3ch)
    byte *pData;    ///< 이미지 데이터 저장 포인터 Array pointer of image data
};

#endif // SUHANVISION_H
