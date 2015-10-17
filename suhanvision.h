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
#include "suhanpoint.h"
#include <vector>

typedef unsigned char byte;

/**
 * @brief The IMG_TYPE enum
 * 이미지 형식 지정
 */
enum IMG_TYPE {BINARY,GRAY,HSV,RGB,BGR};


using SuhanMath::Max;
using SuhanMath::Min;
using SuhanMath::Square;
using SuhanMath::PI;
using SuhanMath::SHMatrix;
using SuhanMath::SHArray;

using namespace std;

/**
 * @brief The SHImage class
 * 수한 이미지 클래스 Suhan Image Class
 */
class SHImage : SHArray<int>
{
public:
    // 생 소멸자
    SHImage();
    SHImage(int col, int row, int channel = 1, IMG_TYPE type = GRAY);
    SHImage(const SHImage& ref);
    SHImage(KImageColor& refImage);
    SHImage(KImageGray& refImage);
    virtual ~SHImage();

    // Management data 데이터 관리
    void Create(int col, int row, int channel = 1, IMG_TYPE type = GRAY);
    ///< 이미지 데이터 생성 함수 (pData의 메모리를 할당해 줌) This function is used for allocating "pData"
    void Reset();   ///< 현재 가지고 있는 데이터 초기화 함수 This function is used for clearing all data of this class.

    // kfc library converter kfc 형변환
    void ToKImageColor(KImageColor& dstImage);
    void ToKImageGray(KImageGray& dstImage);

    // 연산자 오버로딩
    int& operator[] (const int& i) { return pData[i]; }
    int& operator[] (const int& i) const { return pData[i]; }
    SHImage& operator =(const SHImage& src);

    // Getter 데이터 게터
    int Row() const { return nRow; }
    int Col() const { return nCol; }
    int Threshold() const { return nThreshold; }
    IMG_TYPE Type() const { return nType; }
    int Channel() const { return nChannel; }

    // 특수 게터
    int safeGetPixel(int col, int row, int channel = 0)
    {
        if(col < 0) return 0;
        if(row < 0) return 0;
        if(col >= nCol) return 0;
        if(row >= nRow) return 0;
        return pData[ ( col + (row * nCol) ) * nChannel + channel];
    }

    int& pixel(int col, int row, int channel = 0)
    {  return pData[ ( col + (row * nCol) ) * nChannel + channel]; }

    // 데이터 부분관리
    void FillData(int channel, int data);

    // 히스토그램 관련 함수
    void MakeHistogram();
    void MakeHistogram(SHImage& dstHisto);

    // 이미지 색상공간 변환
    void ConvertToHSV(SHImage& dstImage);
    void ConvertToRGB(SHImage& dstImage);
    void ConvertToBGR(SHImage& dstImage);

    // 응용 함수
    void SpeiaToneTransform(double dHue, double dSat);
    void OtsuThresholding();
    void Binarization();
    void BinaryErode(int nMaskSize);
    void BinaryDilate(int nMaskSize);
    void AddGaussianNoise(double dSigma);
    void AddSlatPepperNoise(double dSigma);
    void GaussianSmoothing(double dSigma);
    void MedianFiltering(int nMask);

protected:
    IMG_TYPE nType; ///< 현재 이미지 타입 Current image type
    int nCol;       ///< 칸 수 (너비) a number of cols of image
    int nRow;       ///< 줄 수 (높이) a number of rows of image
    int nChannel;   ///< 채널 개수 a number of channel (1ch, 2ch, 3ch)

    int *pHistogramData; ///< Histogram 데이터가 저장되는 곳
    bool bHistogram; ///< Histogram 데이터가 생성되었는지 나타내 줌
    int nThreshold; ///< Otsu's Threshold 사용 시 Threshold 값
};



class SHLabel : public SHImage
{
public:
    SHLabel() : SHImage() {}
    ~SHLabel() {}

    void SequentialLabeling(SHImage &refImage);
    int FindTopEquLabel(int index);
    void SyncEquList();
    void UpdateVectorList(vector<int> &refVector, int index);
    void PrintBlobColor(SHImage &colorImage, int nMinSize);
private:
    //vector<int> vEquList;
    vector< vector<int> > vvEquList;
    vector< vector<Point> > vvLabelList;
};

#endif // SUHANVISION_H
