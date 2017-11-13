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
#include <stack>
#include "QDebug"

typedef unsigned char byte;

/**
 * @brief The IMG_TYPE enum \n
 * 이미지 형식 지정
 */
enum IMG_TYPE {BINARY,GRAY,HSV,RGB,BGR};


using SuhanMath::Max;
using SuhanMath::Min;
using SuhanMath::Square;
using SuhanMath::PI;
using SuhanMath::SHMatrix;
using SuhanMath::SHArray;
using SuhanMath::SHRect;

using namespace std;


struct SHRGB
{
    SHRGB(int r, int g, int b) : R(r), G(g), B(b) {}
    SHRGB(double h, double s, double v) {


        // 0~1, 0~180으로 스케일된 값
        double h_scale = h;
        double s_scale = s;
        double v_scale = v;

        // Chroma의 값
        double C = v_scale * s_scale;

        // H'의 값
        double H_ = h_scale / 60;

        // H'의 범위를 Case로 구분짓는 identity
        int idt = H_;

        // H' mod 2를 연산할 때 floating point에 %연산자 사용이 불가하므로
        // int로 나머지를 구한 뒤 소숫점을 더해주기 위해 소숫점 아랫부분 값을 먼저 구함
        double ff = H_ - idt;

        // X = C(1- |H' mod 2 - 1|) 공식으로 Intermediate 값으로 사용. 두번째 큰 값
        double X = C * (1 - fabs( idt % 2 + ff - 1 ) );

        // R,G,B 모두 더해주는 값
        double m = v_scale - C;

        // 0~1범위의 R' G' B'
        double R_ = 0.0, G_ = 0.0, B_ = 0.0;


        switch(idt)
        {
        case 0: R_ = C; G_ = X; break;
        case 1: R_ = X; G_ = C; break;
        case 2: G_ = C; B_ = X; break;
        case 3: G_ = X; B_ = C; break;
        case 4: R_ = X; B_ = C; break;
        case 5: R_ = C; B_ = X; break;
        default:
            break;
        }

        R = (R_+m) * 255;
        G = (G_+m) * 255;
        B = (B_+m) * 255;

    }

    int R;
    int G;
    int B;
};

/**
 * @brief The SHImage class \n
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

    // 사이즈 변환
    void ToHalfSize(SHImage& dstImage);

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

    void DrawLine(double rho, double theta);
    void DrawCircle(int x, int y, int r);
    void DrawX(int x, int y);
    void DrawLineRGB(int x1, int y1, int x2, int y2, SHRGB rgb);
    void DrawArrowRGB(int x1, int y1, int x2, int y2, double length, SHRGB rgb);

    // 특수 게터
    int safeGetPixel(int col, int row, int channel = 0)
    {
        if(col < 0) return 0;
        if(row < 0) return 0;
        if(col >= nCol) return 0;
        if(row >= nRow) return 0;
        return pData[ ( col + (row * nCol) ) * nChannel + channel];
    }
    void safeSetPixel(int data, int col, int row, int channel = 0)
    {
        if(col < 0) return ;
        if(row < 0) return ;
        if(col >= nCol) return ;
        if(row >= nRow) return ;
        if(channel < 0) return ;
        if(channel >= nChannel) return ;
        pixel(col,row,channel) = data;
    }

    int& pixel(int col, int row, int channel = 0)
    {  return pData[ ( col + (row * nCol) ) * nChannel + channel]; }
    int& pixel(int col, int row, int channel = 0) const
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


/**
 * @brief The SHLabel class \n
 * Sequntial Labeling을 구현하는 SHLabel 클래스 \n
 * SHImage를 상속받아 내부 이미지 데이터를 Labeling하는 데에 사용됨.
 */
class SHLabel : public SHImage
{
public:
    SHLabel() : SHImage() {}
    virtual ~SHLabel() {}

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

/**
 * @brief The SHEdgeData struct \n
 * SHEdge클래스에 사용되는 EdgeData 구조체 \n
 * 좌표 및 방향, Magnitude가 포함되어 있음.
 */
struct SHEdgeData{
    int u;
    int v;
    int nAng;//0-360
    int nDir;//0,1,2,3
    double dBuf;
    double dMag;
};

/**
 * @brief The SHEdge class \n
 * Edge 검출용 SHEdge클래스. \n
 * SHImage를 상속받아 내부 이미지 데이터를 Edge결과물로써 사용.
 */
class SHEdge : public SHImage
{
public:
    SHEdge() : SHImage(), pEdgeData(NULL) {}
    virtual ~SHEdge() { if(pEdgeData) delete pEdgeData; }

    void Init(double dSigmaVal, int nMaskLength);
    void Canny(double dLow, double dHigh, SHImage imgRef);

    SHEdgeData& EdgePixel(int row, int col);
    vector<SHEdgeData> vEdgeData;

private:

    int nHalf;
    int nFull;

    double dSigma;
    SHRect rtArea;

    vector<SHEdgeData> vEdgeTemp;
    SHMatrix<double> mKernelX;
    SHMatrix<double> mKernelY;

    SHEdgeData *pEdgeData;
};

/**
 * @brief The HoughLine struct \n
 * Hough변환에 사용되는 Line Parameter 구조체
 */
struct HoughLine
{
    double rho;
    double theta;
};
/**
 * @brief The HoughCircle struct \n
 * Hough변환에 사용되는 Circle Parameter 구조체
 */
struct HoughCircle
{
    double x;
    double y;
    double rad;
};

/**
 * @brief The SHHough class \n
 * 하프 변환을 수행하는 클래스.
 */
class SHHough
{
public:
    SHHough() : pCircleSpace(NULL), pLineSpace(NULL), nRhoSize(0), nThetaSize(0) {}
    ~SHHough() { if(pLineSpace) delete pLineSpace; if(pCircleSpace) delete pCircleSpace; }
    void FindLine(SHEdge& edImage, double dRhoRes, double dThetaRes, int nThreshold);
    void FindCircle(SHEdge& edImage, double dRes, int nMinRad, int nMaxRad, int nThreshold);
    int &LineSpace(int rho, int theta) { return pLineSpace[rho * nThetaSize + theta]; }
    int &CircleSpace(int x, int y, int rad) { return pCircleSpace[ x * nYSize * nRadSize + y*nRadSize +rad]; }

    vector<HoughLine> vLineList; ///< Hough변환 Line 결과 저장 벡터
    vector<HoughCircle> vCircleList; ///< Hough변환 Circle 결과 저장 벡터

private:
    double dRhoResolution;
    double dThetaResolution;

    int nRhoSize;
    int nThetaSize;

    int nXSize;
    int nYSize;
    int nRadSize;

    int *pCircleSpace;
    int *pLineSpace;
};

struct SHCornerData
{
    int u;
    int v;
    double R;
};
struct SHCornerImage
{
    double dGradX;
    double dGradY;
    double dGradXX;
    double dGradXY;
    double dGradYY;
    double dR;
};

class SHCorner : public vector<SHCornerData>
{
public:
    SHCorner() : vector(), pCornerData(NULL) {}
    SHCorner(const double& dSigmaValue, const int& nBlockSize)
        :  vector(), pCornerData(NULL) { Init(dSigmaValue,nBlockSize); }

    SHCornerImage& CornerPixel(int row, int col);
    void Init(const double& dSigmaValue, const int& nBlockSize);
    void HarrisCorner(const double& dThresh, const SHImage& imgIn);

private:
    int nCol;
    int nRow;

    int nHalf;
    int nFull;
    int nHalfBlock;
    int nFullBlock;

    double dSigma;

    SHMatrix<double> mKernelX;
    SHMatrix<double> mKernelY;
    SHMatrix<double> mBlockWeight;
    vector<SHCornerData> vTemp;
    SHCornerImage *pCornerData;
};

class SHGaussianPiramid : public vector<SHImage>
{
private:
    double dSigma;
    int nOctave;
public:
    SHGaussianPiramid() : vector() {}
    SHGaussianPiramid(const SHImage& imgOrigin, double dSigma=1.4, int nOctave=0)
        : vector(), dSigma(dSigma), nOctave(nOctave)
    { Make(imgOrigin,dSigma,nOctave); }

    int Make(const SHImage& imgOrigin, double dSigma=1.4, int nOctave=0);
};

struct SHOpticalData
{
    SHOpticalData& operator = (int n);
    double dDx;
    double dDy;
    double dDxDx;
    double dDxDy;
    double dDyDy;

    double dDt;
    double dDtDx;
    double dDtDy;

    double u;
    double v;
};

class SHOpticalFlow : SHArray<SHOpticalData>
{
private:
    int nRow;
    int nCol;
    int nFull;
    int nHalf;

public:
    SHOpticalFlow() {}
    SHOpticalData& pixel(int col, int row)
    {  return pData[col + (row * nCol)]; }
    void findFlow(SHImage& imgT0, SHImage& imgT1, int nWindows = 3);
    void findFlowPyramid(SHImage& imgT0, SHImage& imgT1, int nWindows = 3);
    int Row() const { return nRow; }
    int Col() const { return nCol; }

};


#endif // SUHANVISION_H
