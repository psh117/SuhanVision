/**
  @file suhanvision.cpp
  @author ParkSuhan (psh117@gmail.com)
  @brief Suhan Vision Process Library

  */

#include "suhanvision.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <array>

// 생성자

/**
 * @brief SHImage::SHImage
 */
SHImage::SHImage() : pHistogramData(NULL) {}

/**
 * @brief SHImage::SHImage
 * @param col 생성할 이미지의 칸 수 (너비)
 * @param row 생성할 이미지의 줄 수 (높이)
 * @param channel 생성할 이미지의 채널 수 (채널)
 * @param type 생성할 이미지의 타입 (RGB, HSV 등)
 */
SHImage::SHImage(int col, int row, int channel, IMG_TYPE type)
    : SHArray()
{
    SHImage();
    pHistogramData = NULL;
    nType = type;
    Create(col,row,channel);
}
/**
 * @brief SHImage::SHImage
 * @param refImage 정문호교수님 kfc라이브러리의 KImageColor를 가져올 수 있는 파라미터
 */
SHImage::SHImage(KImageColor& refImage)
{
    SHImage();
    pHistogramData = NULL;

    int col,row;
    col = refImage.Col();
    row = refImage.Row();

    Create(col, row,3);

    int i,j;

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            int nPtr = i*col + j;
            pData[nPtr * 3 + 0] = refImage[i][j].r;
            pData[nPtr * 3 + 1] = refImage[i][j].g;
            pData[nPtr * 3 + 2] = refImage[i][j].b;
        }
    }
    nType = RGB;
}


SHImage::SHImage(const SHImage& ref)
{
    SHImage();
    pHistogramData = NULL;
    Create(ref.Col(),ref.Row(),ref.Channel(),ref.Type());
    nThreshold = ref.Threshold();

    SHArray::operator =(ref);
}

SHImage::SHImage(KImageGray& refImage)
{
    SHImage();
    pHistogramData = NULL;

    int col,row;
    col = refImage.Col();
    row = refImage.Row();

    Create(col, row,1);

    int i,j;

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            int nPtr = i*col + j;
            pData[nPtr] = refImage[i][j];
        }
    }
    nType = GRAY;
}

/**
 * @brief SHImage::~SHImage
 * 소멸자
 */
SHImage::~SHImage()
{
    if(pHistogramData != NULL)
        delete pHistogramData;
}

SHImage& SHImage::operator =(const SHImage& src)
{
    Create(src.Col(),src.Row(),src.Channel(),src.Type());
    nThreshold = src.Threshold();

    SHArray::operator =(src);
    return *this;
}

/**
 * @brief SHImage::Create
 * @param col 생성할 이미지의 칸 수 (너비)
 * @param row 생성할 이미지의 줄 수 (높이)
 * @param channel 생성할 이미지의 채널 수 (채널)
 * @param type 생성할 이미지의 타입 (RGB, HSV 등)
 */
void SHImage::Create(int col, int row, int channel, IMG_TYPE type)
{

    if(pHistogramData != NULL)
        delete pHistogramData;

    nCol = col;
    nRow = row;

    nChannel = channel;
    nType = type;
    nThreshold = -1;

    bHistogram = false;

    int length = nCol * nRow * nChannel;
    SHArray::Create(length);

}

/**
 * @brief SHImage::Reset
 */
void SHImage::Reset()
{
    SHArray::Reset();

    nCol = 0;
    nRow = 0;
    nChannel = 0;
    nType = GRAY;
}


/**
 * @brief SHImage::ToKImageColor
 * @param dstImage 정문호교수님 kfc라이브러리의 KImageColor로 변환할 수 있는 파라미터
 * @todo HSV, GRAY 등 각종 케이스 만들어야함
 */
void SHImage::ToKImageColor(KImageColor& dstImage)
{
    dstImage.Create(nRow,nCol);

    SHImage tmp;
    int i,j;
    switch (nType)
    {
    case RGB:
        if(nChannel != 3) return; // ERROR!

        for(i=0;i<nRow;i++)
        {
            for(j=0;j<nCol;j++)
            {
                int nPtr = i*nCol + j;
                dstImage[i][j].r = pData[nPtr * 3 + 0];
                dstImage[i][j].g = pData[nPtr * 3 + 1];
                dstImage[i][j].b = pData[nPtr * 3 + 2];
            }
        }
        break;
    case BGR:
        if(nChannel != 3) return; // ERROR!

        for(i=0;i<nRow;i++)
        {
            for(j=0;j<nCol;j++)
            {
                int nPtr = i*nCol + j;
                dstImage[i][j].b = pData[nPtr * 3 + 0];
                dstImage[i][j].g = pData[nPtr * 3 + 1];
                dstImage[i][j].r = pData[nPtr * 3 + 2];
            }
        }
        break;

    case HSV:
        if(nChannel != 3) return; // ERROR!
        ConvertToRGB(tmp);

        for(i=0;i<nRow;i++)
        {
            for(j=0;j<nCol;j++)
            {
                int nPtr = i*nCol + j;
                dstImage[i][j].r = (unsigned char)tmp[nPtr * 3 + 0];
                dstImage[i][j].g = (unsigned char)tmp[nPtr * 3 + 1];
                dstImage[i][j].b = (unsigned char)tmp[nPtr * 3 + 2];
            }
        }


        break;
    case GRAY:

        if(nChannel != 1) return; // ERROR!
        break;
    default:
        break;
    }

}
void SHImage::ToKImageGray(KImageGray& dstImage)
{
    if(nChannel != 1) return; // ERROR!

    dstImage.Create(nRow,nCol);

    int i,j;
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            int nPtr = i*nCol + j;
            dstImage[i][j] = (unsigned char)pData[nPtr];
            //printf("%d ",dstImage[i][j]);
        }
        //printf("\r\n");
    }
    //printf("\r\n");
}

void SHImage::ToHalfSize(SHImage& dstImage)
{
    int row = nRow >> 1;
    int col = nCol >> 1;
    dstImage.Create(col,row,nChannel,nType);

    int i,j,k;
    for(i=0;i<nChannel;i++)
    {
        for(j=0;j<row;j++)
        {
            for(k=0;k<col;k++)
            {
                dstImage.pixel(k,j,i) = pixel(k<<1,j<<1,i);
            }
        }
    }
}


/**
 * @brief SHImage::ConvertToRGB
 * @param dstImage RGB색공간으로 변환될 목표 이미지
 */
void SHImage::ConvertToRGB(SHImage& dstImage)
{
    dstImage.Create(nCol,nRow,3,RGB);

    int i,j;

    switch (nType)
    {
    case GRAY:
        for(i=0;i<nRow;i++) {
            for(j=0;j<nCol;j++) {
                int nPtr = i*nCol + j;
                if(nType == GRAY) {
                    dstImage[nPtr * 3 + 0] = pData[nPtr];
                    dstImage[nPtr * 3 + 1] = pData[nPtr];
                    dstImage[nPtr * 3 + 2] = pData[nPtr];
        }}}
        break;

    case HSV:
        for(i=0;i<nRow;i++) {
            for(j=0;j<nCol;j++) {
                int h,s,v;
                int r,g,b;

                int nPtr = i*nCol + j;
                h = pData[nPtr * 3 + 0];
                s = pData[nPtr * 3 + 1];
                v = pData[nPtr * 3 + 2];

                // 0~1, 0~180으로 스케일된 값
                double h_scale = h * 2.0;
                double s_scale = s / 255.0;
                double v_scale = v / 255.0;

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

                r = (R_+m) * 255;
                g = (G_+m) * 255;
                b = (B_+m) * 255;

                dstImage[nPtr * 3 + 0] = r;
                dstImage[nPtr * 3 + 1] = g;
                dstImage[nPtr * 3 + 2] = b;
            }
        }
        break;
    default:
        break;
    }

}

/**
 * @brief SHImage::ConvertToHSV
 * @param dstImage HSV색공간으로 변환될 목표 이미지
 */
void SHImage::ConvertToHSV(SHImage& dstImage)
{
    dstImage.Create(nCol,nRow,3,HSV);

    int i,j;

    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            int h,s,v;
            int r,g,b;

            int nPtr = i*nCol + j;

            switch (nType)
            {
            case RGB:
                r = pData[nPtr * 3 + 0];
                g = pData[nPtr * 3 + 1];
                b = pData[nPtr * 3 + 2];
                break;
            case BGR:
                b = pData[nPtr * 3 + 0];
                g = pData[nPtr * 3 + 1];
                r = pData[nPtr * 3 + 2];
                break;
            default:
                break;
            }

            // 최대, 최소, 차이 구하기
            int rgbMax = Max<int>(3,r,g,b);
            int rgbMin = Min<int>(3,r,g,b);
            int rgbDist = rgbMax - rgbMin;

            // Hue 를 구하기 위한 Raw Data
            double h_raw;
            if(rgbMax == r)
            {
                h_raw = 60 * ( (g-b) / (double)rgbDist ) + 0;
            }
            else if(rgbMax == g)
            {
                h_raw = 60 * ( (b-r) / (double)rgbDist ) + 120;
            }
            else
            {
                h_raw = 60 * ( (r-g) / (double)rgbDist ) + 240;
            }
            // 음수일 경우 360을 더해줌
            if(h_raw < 0) h_raw =+ 360;
            // 8비트로 하기 때문에 180으로 스케일
            h = int(h_raw / 2);
            if(rgbMax == 0)
                s = 0;
            else
                s = ((rgbMax - rgbMin) / (double)rgbMax) * 255;
            v = rgbMax;
            dstImage[nPtr * 3 + 0] = h;
            dstImage[nPtr * 3 + 1] = s;
            dstImage[nPtr * 3 + 2] = v;

        }
    }
}

/**
 * @brief SHImage::FillData
 * @param channel 데이터를 채울 채널
 * @param data 해당하는 채널에 쓸 데이터
 */
void SHImage::FillData(int channel, int data)
{
    int length = nCol * nRow;
    for(int i=0; i<length; i++)
    {
        pData[i*nChannel + channel] = data;
    }
}



void SHImage::MakeHistogram()
{
    if(pHistogramData != NULL)
        delete pHistogramData;
    bHistogram = true;

    int lLen = 256 * nChannel;
    pHistogramData = new int[lLen];

    // pHistogramData는 각 채널순서로 저장됨
    // 0~255 -> 0채널, 256~ 512 -> 1채널 등
    int i,j,k;
    for(i=0;i<lLen;i++)
    {
        pHistogramData[i] = 0;
    }
    for(i=0;i<nChannel;i++)
    {
        for(j=0;j<nRow;j++)
        {
            for(k=0;k<nCol;k++)
            {
                int len = j * nCol + k + i;
                pHistogramData[i*256 + pData[len]]++;
            }
        }
    }
}

void SHImage::MakeHistogram(SHImage& dstHisto)
{
    dstHisto.Create(256,256,1,GRAY);
    MakeHistogram();

    int i,j;
    int lMax = 0;
    for(i=0;i<256;i++)
    {
        if(lMax < pHistogramData[i]) lMax = pHistogramData[i];
    }
    for(i=0;i<256;i++) // col
    {
        if(i == nThreshold)
        {
            for(j=0;j<256;j++)
            {
                dstHisto[j*256+i] = 50;
            }
            continue;
        }
        int lRow = pHistogramData[i] * 256 / (double)lMax;
        for(j=0;j<256;j++) // row
        {
            if((256-j)<lRow)
                dstHisto[j*256 + i] = 0;
            else
                dstHisto[j*256 + i] = 255;
        }
    }
}

void SHImage::SpeiaToneTransform(double dHue, double dSat)
{
    int nHue = dHue / 2;
    int nSat = dSat * 255;

    SHImage imgHSV; // HSV 용

    switch (nType)
    {
    case RGB:
        ConvertToHSV(imgHSV); // 오리지날을 HSV로 컨버팅
        imgHSV.FillData(0,nHue); // HUE 채널 덮기
        imgHSV.FillData(1,nSat); // SAT 채널 덮기

        imgHSV.ConvertToRGB(*this);

        break;
    case BGR:
        break;
    case HSV:
        break;
    default:
        break;
    }

    if(nType == RGB)
    {
    }
}

void SHImage::OtsuThresholding()
{
    MakeHistogram();

    double dQ1[256], dMyu1[256], dMyu2[256];
    double dP[256];    // nomalized data
    int i,j;
    int lTotal = 0;
    double dSigmaB, dMaxSigma = 0.0;
    double dMyu = 0.0;
    for(i=0;i<256;i++)
    {
        dMyu1[i] = 0;
        dMyu2[i] = 1.0;
        lTotal += pHistogramData[i];        // sum
        //dP[i] = pHistogramData[i] / 256.0;    // nomalize
    }
    for(i=0;i<256;i++)
    {
        dP[i] = pHistogramData[i] / (double)lTotal;
        dMyu += dP[i] * i;
    }
    //double dMyu = lTotal / 256.0;

    dQ1[0] = dP[0];
    dMyu1[0] = 0.0;

    for(i=0; i<255;i++)
    {
        dQ1[i+1] = dQ1[i] + dP[i+1];
        if(dQ1[i+1] >= 1.0) // over percentage
            break;

        if(dQ1[i+1] == 0.0) // sum is still 0.0
        {
            dMyu1[i+1] = 0.0;   // mean is also 0.0
            continue;           // and check again
        }

        // Otsu Thresholding (Computer vision lecture #6 23p)
        // Recursion
        dMyu1[i+1] = (dQ1[i] * dMyu1[i] + (i+1) * dP[i+1]) / dQ1[i+1];
        dMyu2[i+1] = (dMyu - dQ1[i+1] * dMyu1[i+1]) / (1.0 - dQ1[i+1]);

        dSigmaB = dQ1[i+1] * (1.0 - dQ1[i+1]) * Square<double>(dMyu1[i+1] - dMyu2[i+1]);

        if(dSigmaB > dMaxSigma)
        {
            dMaxSigma = dSigmaB;
            nThreshold = i+1;
        }

    }

}

void SHImage::Binarization()
{
    int i,j;
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            int addr = i*nRow + j;
            if(pData[addr] >= nThreshold)   pData[addr] = 255;
            else                            pData[addr] = 0;
        }
    }
}
void SHImage::BinaryDilate(int nMaskSize)
{
    int i,j,x,y;
    int nHalfMask = nMaskSize / 2;
    SHImage dst;
    dst.Create(nCol,nRow);
    for(i=nHalfMask;i<nRow-nHalfMask;i++)
    {
        for(j=nHalfMask; j<nCol-nHalfMask; j++)
        {
            bool bAll = true;

            for(x=-nHalfMask; x<nHalfMask+1 && bAll; x++)
            {
                for(y=-nHalfMask; y<nHalfMask+1; y++)
                {
                    if(pixel(j+x,i+y) == 0)
                    {
                        bAll = false;
                        break;
                    }
                }
            }
            if(bAll == true)
            {
                dst.pixel(j,i) = 255;
            }
        }
    }
    // copy
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            pixel(j,i) = dst.pixel(j,i);
        }
    }
}

void SHImage::BinaryErode(int nMaskSize)
{
    int i,j,x,y;
    int nHalfMask = nMaskSize / 2;
    SHImage dst;
    dst.Create(nCol,nRow);
    for(i=nHalfMask;i<nRow-nHalfMask;i++)
    {
        for(j=nHalfMask; j<nCol-nHalfMask; j++)
        {
            bool bAll = true;

            for(x=-nHalfMask; x<nHalfMask+1 && bAll; x++)
            {
                for(y=-nHalfMask; y<nHalfMask+1; y++)
                {
                    if(pixel(j+x,i+y) == 255)
                    {
                        bAll = false;
                        break;
                    }
                }
            }
            if(bAll == false)
            {
                dst.pixel(j,i) = 255;
            }
        }
    }
    // copy
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            pixel(j,i) = dst.pixel(j,i);
        }
    }
}

void SHImage::AddGaussianNoise(double dSigma)
{
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,dSigma);

    int i,j,k;
    for(i=0;i<nChannel;i++)
    {
        for(j=0;j<nRow;j++)
        {
            for(k=0;k<nCol;k++)
            {
                pixel(k,j,i) += distribution(generator);
                if(pixel(k,j,i) > 255) pixel(k,j,i) = 255;
                if(pixel(k,j,i) < 0) pixel(k,j,i) = 0;
            }
        }
    }
}

void SHImage::AddSlatPepperNoise(double dSigma)
{
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,dSigma);

    int i,j,k;
    for(j=0;j<nRow;j++)
    {
        for(k=0;k<nCol;k++)
        {
            double dRandNum = distribution(generator);
            if(dRandNum > 10.0)
            {
                for(i=0;i<nChannel;i++)
                {
                    pixel(k,j,i) = 255;
                }
            }
            else if (dRandNum < -10.0)
            {
                for(i=0;i<nChannel;i++)
                {
                    pixel(k,j,i) = 0;
                }
            }
        }
    }
}



void SHImage::GaussianSmoothing(double dSigma)
{
    int i,j,k;
    int x, y;

    // 마스크 만들기
    int nHalf = (int)(dSigma*2); // 95%
    int nFull = 2 * nHalf + 1;
    SHMatrix<double> mMask(nFull,nFull);
    //double nMask[nFull][nFull];

    double dScale = 0.0;
    double dSigma2 = Square<double>(dSigma);

    double dConst = 1.0 / PI * 2 / dSigma2;

    y = -nHalf;
    for(i=0;i<nFull;i++)    // y
    {
        x = -nHalf;
        for(j=0;j<nFull;j++)    // x
        {
            mMask.at(i,j) = dConst*
                    exp((Square<int>(x) + Square<int>(y))
                        * -1.0 / 2.0 / dSigma2);
            dScale += mMask.at(i,j);
            x++;
        }
        y++;
    }
    mMask /= dScale;

    // 콘볼루쟌
    SHImage dst(*this);
    double dValue;
    for(i=nHalf;i<nRow-nHalf;i++)
    {
        for(j=nHalf;j<nCol-nHalf;j++)
        {
            for(k=0;k<nChannel;k++)
            {
                dValue = 0.0;
                for(y=-nHalf; y<=nHalf; y++)
                {
                    for(x=-nHalf; x<=nHalf; x++)
                    {
                        dValue += dst.pixel(j+x,i+y,k) * mMask.at(y+nHalf,x+nHalf);
                    }
                }
                pixel(j,i,k) = (int)dValue;
            }
        }
    }
}

void SHImage::MedianFiltering(int nMask)
{
    int i,j,k;
    int x,y;

    int nHalf = nMask / 2;
    // 에러 방지용
    int nFull = nHalf * 2 + 1;
    const int nSize = nFull * nFull;

    vector<int> vMedianVect(nSize);

    SHImage src(*this);
    for(i=nHalf; i<nRow-nHalf; i++)
    {
        for(j=nHalf; j<nCol-nHalf; j++)
        {
            for(k=0;k<nChannel;k++)
            {
                vMedianVect.clear();
                for(y=-nHalf;y<=nHalf;y++)
                {
                    for(x=-nHalf;x<=nHalf;x++)
                    {
                        vMedianVect.push_back(src.pixel(j+x,i+y,k));
                    }
                }
                sort(vMedianVect.begin(),vMedianVect.end());
                pixel(j,i,k) = vMedianVect[nSize / 2];
            }
        }
    }
}

// DRAW LINE - SHIMAGE Class

void SHImage::DrawLine(double rho, double theta)
{
    int x,y;
    for(x=0;x<nCol;x++)
    {
        y = (int)(( rho + x * sin(theta) ) / cos(theta));
        safeSetPixel(230,x,y,0);
        safeSetPixel(230,x+1,y,0);
        safeSetPixel(230,x,y+1,0);
        safeSetPixel(230,x+1,y+1,0);
        safeSetPixel(230,x,y-1,0);
        safeSetPixel(230,x-1,y,0);
        safeSetPixel(230,x-1,y-1,0);
        //safeSetPixel(60,x,y,1);
        //safeSetPixel(60,x,y,2);
    }
}



void SHImage::DrawCircle(int x, int y, int r)
{
    int deg;
    int x_p, y_p;
    for(deg=0;deg<360;deg++)
    {
        x_p = x + r*cos(deg * PI / 180);
        y_p = y + r*sin(deg * PI / 180);
        safeSetPixel(230,x_p,y_p,0);
        safeSetPixel(230,x_p+1,y_p,0);
        safeSetPixel(230,x_p,y_p+1,0);
        safeSetPixel(230,x_p+1,y_p+1,0);
        safeSetPixel(230,x_p,y_p-1,0);
        safeSetPixel(230,x_p-1,y_p,0);
        safeSetPixel(230,x_p-1,y_p-1,0);
    }
}
void SHImage::DrawX(int x, int y)
{
    safeSetPixel(230,x,y);
    for(int i =0;i<5;i++)
    {
        safeSetPixel(230,x+i,y+i);
        safeSetPixel(230,x-i,y+i);
        safeSetPixel(230,x+i,y-i);
        safeSetPixel(230,x-i,y-i);
    }
}
void SHImage::DrawLineRGB(int x1, int y1, int x2, int y2, SHRGB rgb)
{
    int dx = x2-x1;
    int dy = y2-y1;
    int i;
    int x,y;
    if(dx == 0 && dy == 0) return;

    //qDebug("dx == %d, dy == %d",dx,dy);
    if(abs(dx) < abs(dy))
    {
        for(i=y1; i!=y2; i+=(dy/abs(dy)) )
        {
            y = i;
            x = x1 + ((double)dx / dy) * (i-y1);
            safeSetPixel(rgb.R,x,y,0);
            safeSetPixel(rgb.G,x,y,1);
            safeSetPixel(rgb.B,x,y,2);
        }
        y = i;
        x = x1 + ((double)dx / dy) * (i-y1);
        // 마지막
        safeSetPixel(rgb.R,x,y,0);
        safeSetPixel(rgb.G,x,y,1);
        safeSetPixel(rgb.B,x,y,2);
    }
    else
    {

        for(i=x1; i!=x2; i+=(dx/abs(dx)) )
        {
            x = i;
            y = y1 + ((double)dy / dx) * (i-x1);
            safeSetPixel(rgb.R,x,y,0);
            safeSetPixel(rgb.G,x,y,1);
            safeSetPixel(rgb.B,x,y,2);
        }
        x = i;
        y = y1 + ((double)dy / dx) * (i-x1);
        // 마지막
        safeSetPixel(rgb.R,x,y,0);
        safeSetPixel(rgb.G,x,y,1);
        safeSetPixel(rgb.B,x,y,2);
    }
}

void SHImage::DrawArrowRGB(int x1, int y1, int x2, int y2, double length, SHRGB rgb)
{
    //qDebug("Draw Line");
    DrawLineRGB(x1,y1,x2,y2,rgb);

    double deg = atan2((y2 - y1), (x2 - x1));
    double deg_dst;
    //x2 -= cos(deg) * 12;
    //y2 -= sin(deg) * 12;

    //qDebug("Draw 1");
    deg_dst = deg + 30 * (3.141592) / 180 ;
    x1 = x2 - cos(deg_dst) * length;
    y1 = y2 - sin(deg_dst) * length;
    DrawLineRGB(x1,y1,x2,y2,rgb);

    //qDebug("Draw 2");
    deg_dst = deg - 30 * (3.141592) / 180;
    x1 = x2 - cos(deg_dst) * length;
    y1 = y2 - sin(deg_dst) * length;
    DrawLineRGB(x1,y1,x2,y2,rgb);

}

void SHLabel::UpdateVectorList(vector<int> &refVector, int index)
{
    unsigned int i,j;

    for(i=0;i<vvEquList[index].size();i++)
    {
        bool notExist = true;
        for(j=0;j<refVector.size();j++)
        {
            if(refVector[j] == vvEquList[index][i])
            {
                notExist = false;
            }

        }
        if(notExist == true)
        {
            refVector.push_back(vvEquList[index][i]);
            UpdateVectorList(refVector, vvEquList[index][i]);
            //isChanged = true;
        }
    }
}

void SHLabel::SyncEquList()
{
    unsigned int i,j;
    for(i=1;i<vvEquList.size();i++) // total member
    {
        //vector<int> total = vvEquList[i];
        for(j=0;j<vvEquList[i].size();j++)  // member's total member
        {
            UpdateVectorList(vvEquList[i],vvEquList[i][j]);
        }
        sort(vvEquList[i].begin(),vvEquList[i].end());
    }
}


int SHLabel::FindTopEquLabel(int index)
{
    int nPastIndex = index;
    int nFindIndex = index;

    while(true)
    {
        //nFindIndex = vEquList[nFindIndex];
        if(nFindIndex == nPastIndex)
            return nFindIndex;

        nPastIndex = nFindIndex;
    }
}



void SHLabel::SequentialLabeling(SHImage &refImage)
{
    Create(refImage.Col(),refImage.Row());

    int i,j,k;
    int nLabelCount = 1;
    vvEquList.push_back(vector<int>());  // Equ 0 = 0
    // Label은 1부터
    // 0은 Background
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            // Do nothing and proceed to next pixel
            if(refImage.pixel(j,i) == 0) continue;
            // If not background
            else if (safeGetPixel(j-1,i-1) != 0)
            {
                pixel(j,i) = safeGetPixel(j-1,i-1);
            }
            else if (safeGetPixel(j-1,i) == 0 && safeGetPixel(j,i-1) == 0)
            {
                pixel(j,i) = nLabelCount;
                vvEquList.push_back(vector<int>());
                // 본인
                vvEquList[nLabelCount].push_back(nLabelCount);
                nLabelCount++;
            }
            else if (safeGetPixel(j-1,i) != 0 && safeGetPixel(j,i-1) != 0)
            {

                pixel(j,i) = safeGetPixel(j-1,i);
                if(safeGetPixel(j-1,i) != safeGetPixel(j,i-1))
                {
                    vvEquList[safeGetPixel(j-1,i)].push_back(safeGetPixel(j,i-1));
                    vvEquList[safeGetPixel(j,i-1)].push_back(safeGetPixel(j-1,i));
                    //SyncEquList(safeGetPixel(j-1,i),safeGetPixel(j,i-1));
                }
            }
            else if (safeGetPixel(j-1,i) != 0)
            {
                pixel(j,i) = safeGetPixel(j-1,i);
            }
            else if (safeGetPixel(j,i-1) != 0)
            {
                pixel(j,i) = safeGetPixel(j,i-1);
            }
        }
    }
    SyncEquList();

    vector<int> vLabelList;
    vector<int>::iterator iter;

    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            if(pixel(j,i) == 0) continue;

            iter = find(vLabelList.begin(),vLabelList.end(),vvEquList[pixel(j,i)][0]);
            if(iter == vLabelList.end())
            {
                vLabelList.push_back(vvEquList[pixel(j,i)][0]);
                vvLabelList.push_back(vector<Point>());
            }

            pixel(j,i)=vvEquList[pixel(j,i)][0];
        }
    }
    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            if(pixel(j,i) == 0) continue;
            int labelNum = 0;
            for(k=0;k<vLabelList.size();k++)
            {
                if(pixel(j,i) == vLabelList[k])
                {
                    labelNum = k + 1;
                    break;
                }
            }
            pixel(j,i) = labelNum;
            vvLabelList[labelNum-1].push_back(Point(j,i));
        }
    }

}

void SHLabel::PrintBlobColor(SHImage &colorImage, int nMinSize)
{
    colorImage.Create(nCol,nRow,3,RGB);
    unsigned int i,j;
    for(i=0;i<vvLabelList.size();i++)
    {
        int r = rand() % 256;
        int g = rand() % 256;
        int b = rand() % 256;

        if(vvLabelList[i].size() >= nMinSize)
        {
            for(j=0;j<vvLabelList[i].size();j++)
            {
                colorImage.pixel(vvLabelList[i][j].x,vvLabelList[i][j].y,0) = r;
                colorImage.pixel(vvLabelList[i][j].x,vvLabelList[i][j].y,1) = g;
                colorImage.pixel(vvLabelList[i][j].x,vvLabelList[i][j].y,2) = b;
            }
        }
    }
}


/* SHEdge */

void SHEdge::Init(double dSigmaVal, int nMaskLength)
{
    // 내부 변수
    dSigma = dSigmaVal;

    // 마스크 만들기
    nHalf = nMaskLength / 2;
    nFull = nHalf * 2 + 1;

    //create the conv. mask
    mKernelX.Create(nFull,nFull);
    mKernelY.Create(nFull,nFull);

    //compute the mask
    int i,j;
    int x,y;
    double dTmp,dScale=0.0;
    double dSigma2 = Square<double>(dSigma);

    // FDG
    y=-nHalf;
    for(i=0;i<nFull;i++)
    {
        x=-nHalf;
        dTmp = -y*exp(-(y*y)) / (2*dSigma2);
        for(j=0;j<nFull;j++)
        {
            mKernelY.at(i,j) = dTmp*exp(-(x*x)/dSigma2*PI);
            mKernelX.at(j,i) = mKernelY.at(i,j);
            if(mKernelY.at(i,j) > 0)
                dScale += mKernelY.at(i,j);
            x++;
        }
        y++;
    }
    for(i=0; i<nMaskLength; i++)
    {
        for(j=0; j<nMaskLength; j++)
        {
            mKernelY.at(i,j) /= dScale;
            mKernelX.at(i,j) /= dScale;
        }
    }
    /*
    printf("\n");
    for(i=0;i<nFull;i++)
    {
        for(j=0;j<nFull;j++)
        {
            printf("%.4lf ", mKernelY.at(i,j));
        }
        printf("\n");
    }
    for(i=0;i<nFull;i++)
    {
        for(j=0;j<nFull;j++)
        {
            printf("%.4lf ", mKernelX.at(i,j));
        }
        printf("\n");
    }
    */
}
SHEdgeData& SHEdge::EdgePixel(int row, int col)
{
    return pEdgeData[row*nCol + col];
}

void SHEdge::Canny(double dLow, double dHigh, SHImage imgRef)
{
    vEdgeData.clear();
    vEdgeTemp.clear();

    SHImage::Create(imgRef.Col(),imgRef.Row());

    pEdgeData = new SHEdgeData[nCol*nRow];

    int i,j,r,c;
    int x,y;
    double dGradX, dGradY, dDeg;

    for(i=nHalf; i<nRow - nHalf; i++)
    {
        for(j=nHalf; j<nCol - nHalf; j++)
        {
            dGradX = 0.0;
            dGradY = 0.0;

            // convolution
            y=-nHalf;
            for(r=0;r<nFull;r++)
            {
                x=-nHalf;
                for(c=0;c<nFull;c++)
                {
                    dGradX += imgRef.pixel(j+x,i+y) * mKernelX.at(r,c);
                    dGradY += imgRef.pixel(j+x,i+y) * mKernelY.at(r,c);
                    x++;
                }
                y++;
            }
            //printf("(%d,%d) Mag = %.4lf, %.4lf ",i,j,dGradX,dGradY);
            // magnitude
            EdgePixel(i,j).dMag = abs(dGradX) + abs(dGradY);
            //printf("  %.4lf\n",EdgePixel(i,j).dMag);
            //EdgePixel(i,j).dMag = sqrt(Square<double>(dGradX) + Square<double>(dGradY));

            // direction
            if(EdgePixel(i,j).dMag > dLow)
            {
                dDeg = atan2(dGradY,dGradX);
                EdgePixel(i,j).nAng = (unsigned short)(dDeg * 180.0 / PI);
                EdgePixel(i,j).nDir = (unsigned char)((((int)(EdgePixel(i,j).nAng/22.5)+1)>>1) & 0x03);
            }
            else
            {
                EdgePixel(i,j).dMag = 0.0;
            }
        }
    }

    // non-maxima suppression

    int         nShiftX[4] = {1,1,0,-1};
    int         nShiftY[4] = {0,1,1,1};
    int         nH=nRow-nHalf-1, nW=nCol-nHalf-1;

    SHEdgeData oEdgeData;

    for(i=nHalf+1; i<nH; i++){
        for(j=nHalf+1; j<nW; j++){
            if(EdgePixel(i,j).dMag == 0.0) continue;


            if(EdgePixel(i,j).dMag > EdgePixel(i + nShiftY[EdgePixel(i,j).nDir],j + nShiftX[EdgePixel(i,j).nDir]).dMag &&
               EdgePixel(i,j).dMag > EdgePixel(i - nShiftY[EdgePixel(i,j).nDir],j - nShiftX[EdgePixel(i,j).nDir]).dMag)
            {
                if(EdgePixel(i,j).dMag > dHigh){
                    oEdgeData.u    = (unsigned short)j;
                    oEdgeData.v    = (unsigned short)i;
                    oEdgeData.nAng = EdgePixel(i,j).nAng;
                    oEdgeData.nDir = EdgePixel(i,j).nDir;
                    oEdgeData.dMag = EdgePixel(i,j).dMag;

                    vEdgeTemp.push_back(oEdgeData);
                }
                EdgePixel(i,j).dBuf = EdgePixel(i,j).dMag;
            }
            else if(EdgePixel(i,j).dMag == EdgePixel(i + nShiftY[EdgePixel(i,j).nDir],j + nShiftX[EdgePixel(i,j).nDir]).dMag) // 같은 값이 2개 연속일 때 예외처리
            {
                if(EdgePixel(i,j).dMag > dHigh){
                    oEdgeData.u    = (unsigned short)j;
                    oEdgeData.v    = (unsigned short)i;
                    oEdgeData.nAng = EdgePixel(i,j).nAng;
                    oEdgeData.nDir = EdgePixel(i,j).nDir;
                    oEdgeData.dMag = EdgePixel(i,j).dMag;

                    vEdgeTemp.push_back(oEdgeData);
                }
                EdgePixel(i,j).dBuf = EdgePixel(i,j).dMag;
            }
        }
    }

    // hysteresis thresholding
    int iy,jx;
    while(vEdgeTemp.empty() == false)
    {
        // store edge list
        oEdgeData = vEdgeTemp.back();
        vEdgeData.push_back(oEdgeData);
        vEdgeTemp.pop_back();

        // get edge coordinate
        jx = oEdgeData.u;
        iy = oEdgeData.v;

        // search neighbor edges
        for(i=-1; i<2; i++)
            for(j=-1; j<2; j++)
                if(EdgePixel(iy+i,jx+j).dBuf && EdgePixel(iy+i,jx+j).dBuf <= dHigh)
                {
                    oEdgeData.u    = (unsigned short)(jx+j);
                    oEdgeData.v    = (unsigned short)(iy+i);
                    oEdgeData.nAng = EdgePixel(i+iy,j+jx).nAng;
                    oEdgeData.nDir = EdgePixel(i+iy,j+jx).nDir;
                    oEdgeData.dMag = EdgePixel(i+iy,j+jx).dBuf;

                    vEdgeTemp.push_back(oEdgeData);

                    EdgePixel(i+iy,j+jx).dBuf = 0.0;
                }
    }

    for(i=0;i<(int)vEdgeData.size();i++)
    {
        pixel(vEdgeData[i].u,vEdgeData[i].v) = 255;
    }

}




void SHHough::FindLine(SHEdge &edImage, double dRhoRes, double dThetaRes, int nThreshold)
{
    if(pLineSpace)
        delete pLineSpace;

    int maxSize;
    maxSize = (edImage.Col() + edImage.Row()) * 2;

    nRhoSize = maxSize / dRhoRes;
    nThetaSize = PI / dThetaRes;

    // Space 생성
    pLineSpace = new int[nRhoSize * nThetaSize];

    int i,j;
    for(i=0;i<(int)edImage.vEdgeData.size();i++)
    {
        for(j=0;j<nThetaSize;j++)
        {
            int x = edImage.vEdgeData[i].u;
            int y = edImage.vEdgeData[i].v;
            double rho = (-x * sin(j * dThetaRes) + y * cos(j * dThetaRes) + maxSize / 2) * (1 / dRhoRes);
            LineSpace((int)(rho),j) += 1;
        }
    }

    vLineList.clear();
    for(i=0;i<nRhoSize;i++)
    {
        for(j=0;j<nThetaSize;j++)
        {
            if(LineSpace(i,j) > nThreshold)
            {
                HoughLine oHL;
                oHL.rho = (i * dRhoRes) - maxSize/2;
                oHL.theta = j * dThetaRes;
                vLineList.push_back(oHL);
            }
        }
    }
}

void SHHough::FindCircle(SHEdge& edImage, double dRes, int nMinRad, int nMaxRad, int nThreshold)
{
    if(pCircleSpace)
        delete pCircleSpace;

    nRadSize = nMaxRad - nMinRad;

    // Space 생성
    nXSize = edImage.Col() / dRes;
    nYSize = edImage.Row() / dRes;
    int nSize = ( nXSize * nYSize ) * nRadSize;
    pCircleSpace = new int[nSize];

    int i,j,k,l;
    for(i=0;i<(int)edImage.vEdgeData.size();i++)
    {
        for(j=1;j<nRadSize;j++)
        {
            int nRad = nMinRad + j;
            double dAng = edImage.vEdgeData[i].nAng * PI / 180;

            int x = (edImage.vEdgeData[i].u - nRad * cos(dAng)) / dRes;
            int y = (edImage.vEdgeData[i].v - nRad * sin(dAng)) / dRes;
            if(x<2) continue;
            if(y<2) continue;

            for(k=-2;k<3;k++)
            {
                for(l=-2;l<3;l++)
                {
                    CircleSpace(x+k,y+l,j-1) += 1;
                    CircleSpace(x+k,y+l,j) += 1;
                    CircleSpace(x+k,y+l,j+1) += 1;
                }
            }
            CircleSpace(x,y,j) += 1;

            x = (edImage.vEdgeData[i].u + nRad * cos(dAng)) / dRes;
            y = (edImage.vEdgeData[i].v + nRad * sin(dAng)) / dRes;

            for(k=-2;k<3;k++)
            {
                for(l=-2;l<3;l++)
                {
                    CircleSpace(x+k,y+l,j-1) += 1;
                    CircleSpace(x+k,y+l,j) += 1;
                    CircleSpace(x+k,y+l,j+1) += 1;
                }
            }
            CircleSpace(x,y,j) += 1;

        }
    }
    vCircleList.clear();
    for(i=0;i<nYSize;i++)
    {
        for(j=0;j<nXSize;j++)
        {
            for(k=0;k<nRadSize;k++)
            {
                if(CircleSpace(j,i,k) > nThreshold)
                {
                    int a,b,c;
                    bool fail = false;
                    for(a=-1;a<2;a++)
                    {
                        for(b=-1;b<2;b++)
                        {
                            for(c=-1;c<2;c++)
                            {
                                if(CircleSpace(j,i,k) <= CircleSpace(j+a,i+b,k+c))
                                {
                                    if(a==0 && b==0 && c==0) continue;
                                    fail = true;
                                }
                            }
                        }
                    }
                    for(a=-1;a<1;a++)
                    {
                        for(b=-1;b<1;b++)
                        {
                            for(c=-1;c<1;c++)
                            {
                                if(CircleSpace(j,i,k) == CircleSpace(j+a,i+b,k+c))
                                {
                                    if(a==0 && b==0 && c==0) continue;
                                    fail = false;
                                }
                            }
                        }
                    }
                    if(fail == false)
                    {

                        HoughCircle oHC;
                        oHC.rad = k+nMinRad;
                        oHC.x = j * dRes;
                        oHC.y = i * dRes;
                        vCircleList.push_back(oHC);
                    }
                }
            }
        }
    }
}



SHCornerImage& SHCorner::CornerPixel(int row, int col)
{
    return pCornerData[row*nCol + col];
}

void SHCorner::Init(const double& dSigmaValue, const int& nBlockSize)
{

    dSigma = dSigmaValue;

    // 마스크 만들기
    nHalf = 3.0 * dSigma;
    nFull = nHalf * 2 + 1;

    //create the conv. mask
    mKernelX.Create(nFull,nFull);
    mKernelY.Create(nFull,nFull);


    //compute the mask
    int i,j;
    int x,y;
    double dTmp,dScale=0.0;
    double dSigma2 = Square<double>(dSigma);


    // FDG
    y=-nHalf;
    for(i=0;i<nFull;i++)
    {
        x=-nHalf;
        dTmp = -y*exp(-(y*y)) / (2*dSigma2);
        for(j=0;j<nFull;j++)
        {
            mKernelY.at(i,j) = dTmp*exp(-(x*x)/dSigma2*PI);
            mKernelX.at(j,i) = mKernelY.at(i,j);
            if(mKernelY.at(i,j) > 0)
                dScale += mKernelY.at(i,j);
            x++;
        }
        y++;
    }
    mKernelX /= dScale;
    mKernelY /= dScale;

    /*
    // debug

    printf("\r\n [mKernelX] \r\n");
    for(i=0;i<nFull;i++)
    {
        for(j=0;j<nFull;j++)
        {
            printf("%.3lf ", mKernelX.at(i,j));
        }
        printf("\r\n");
    }
    */

    nHalfBlock   = nBlockSize/2;
    nFullBlock  = nHalfBlock*2 + 1;

    mBlockWeight.Create(nFullBlock,nFullBlock);

    double dSigmaBlock = nHalfBlock / 3.0;


    dSigma2 = 2.0*Square<double>(dSigmaBlock);
    dScale   = 0.0;

    y=-nHalfBlock;
    for(i=0;i<nFullBlock;i++)
    {
        x=-nHalfBlock;
        for(j=0;j<nFullBlock;j++)
        {
            mBlockWeight.at(i,j) = exp(-(x*x+y*y)/dSigma2);
            dScale += mBlockWeight.at(i,j);
            x++;
        }
        y++;
    }

    mBlockWeight /= dScale;

    /*
    // debug

    printf("\r\n [mBlockWeight] \r\n");
    for(i=0;i<nFullBlock;i++)
    {
        for(j=0;j<nFullBlock;j++)
        {
            printf("%.3lf ", mBlockWeight.at(i,j));
        }
        printf("\r\n");
    }
    */

}

void SHCorner::HarrisCorner(const double& dThresh, const SHImage& imgIn)
{
    if(pCornerData)
        delete pCornerData;

    nCol = imgIn.Col();
    nRow = imgIn.Row();

    pCornerData = new SHCornerImage[nCol * nRow];

    // GradX, GradY Convolution


    int i,j,k,l,x,y;


    // Gradient X Convolution
    for(i=nHalf;i<nRow - nHalf;i++)
    {
        for(j=nHalf; j<nCol-nHalf; j++)
        {
            double dValueTemp = 0.0;
            y = -nHalf;
            for(k=0;k<nFull;k++)
            {
                x = -nHalf;
                for(l=0;l<nFull;l++)
                {
                    dValueTemp += mKernelX.at(k,l) * imgIn.pixel(j+x,i+y);
                    x++;
                }
                y++;
            }
            CornerPixel(i,j).dGradX = dValueTemp;
        }
    }

    // Gradient Y Convolution
    for(i=nHalf;i<nRow - nHalf;i++)
    {
        for(j=nHalf; j<nCol-nHalf; j++)
        {
            double dValueTemp = 0.0;
            y = -nHalf;
            for(k=0;k<nFull;k++)
            {
                x = -nHalf;
                for(l=0;l<nFull;l++)
                {
                    dValueTemp += mKernelY.at(k,l) * imgIn.pixel(j+x,i+y);
                    x++;
                }
                y++;
            }
            CornerPixel(i,j).dGradY = dValueTemp;
        }
    }

    /*
    // debug

    printf("\r\n [dGradX] \r\n");
    for(i=nHalf;i<nRow - nHalf;i++)
    {
        for(j=nHalf;j<nCol - nHalf;j++)
        {
            printf("%5.2lf ", CornerPixel(i,j).dGradX);
        }
        printf("\r\n");
    }

    // debug

    printf("\r\n [dGradY] \r\n");
    for(i=nHalf;i<nRow - nHalf;i++)
    {
        for(j=nHalf;j<nCol - nHalf;j++)
        {
            printf("%5.2lf ", CornerPixel(i,j).dGradY);
        }
        printf("\r\n");
    }
    */


    for(i=nHalf;i<nRow - nHalf;i++)
    {
        for(j=nHalf; j<nCol-nHalf; j++)
        {
            CornerPixel(i,j).dGradXX = CornerPixel(i,j).dGradX * CornerPixel(i,j).dGradX;
            CornerPixel(i,j).dGradXY = CornerPixel(i,j).dGradX * CornerPixel(i,j).dGradY;
            CornerPixel(i,j).dGradYY = CornerPixel(i,j).dGradY * CornerPixel(i,j).dGradY;
        }
    }

    // M Matrix (Computer Vision #12, 9p)
    SHMatrix<double> mM(2,2);

    int nStartX = nHalf + nHalfBlock;
    int nStartY = nHalf + nHalfBlock;
    int nEndX = nCol - nHalf - nHalfBlock;
    int nEndY = nRow - nHalf - nHalfBlock;
    for(i=nStartY;i<nEndY;i++)
    {
        for(j=nStartX; j<nEndX; j++)
        {
            mM.Clear();
            y = -nHalfBlock;
            for(k=0;k<nFullBlock;k++)
            {
                x = -nHalfBlock;
                for(l=0;l<nFullBlock;l++)
                {
                    mM.at(0,0) += CornerPixel(i+y,j+x).dGradXX * mBlockWeight.at(k,l);
                    mM.at(1,0) += CornerPixel(i+y,j+x).dGradXY * mBlockWeight.at(k,l);
                    mM.at(1,1) += CornerPixel(i+y,j+x).dGradYY * mBlockWeight.at(k,l);
                    x++;
                }
                y++;
            }

            mM.at(0,1) = mM.at(1,0);

            double dDetM = (mM.at(0,0) * mM.at(1,1) - Square<double>(mM.at(0,1)));
            double dTrM = (mM.at(0,0) + mM.at(1,1));

            CornerPixel(i,j).dR = dDetM - 0.04 * (dTrM * dTrM);
        }
    }

    /*
    printf("\r\n [dR] \r\n");
    for(i=nStartY;i<nEndY;i++)
    {
        for(j=nStartX; j<nEndX; j++)
        {
            printf("%.2lf ", CornerPixel(i,j).dR);
        }
        printf("\r\n");
    }

    */
    vector::clear();

    SHCornerData oCornerData;

    nStartX++;
    nStartY++;
    nEndX--;
    nEndY--;

    for(i=nStartY;i<nEndY;i++)
    {
        for(j=nStartX;j<nEndX;j++)
        {
            if(CornerPixel(i,j).dR < dThresh) continue;
            bool isNear = false;
            for(k=-1;k<2;k++)
            {
                for(l=-1;l<2;l++)
                {
                    if(k==0 && l==0) continue;
                    if(CornerPixel(i,j).dR < CornerPixel(i+k,j+l).dR)
                        isNear = true;
                }
            }
            if(isNear == true) continue;
            oCornerData.u = j;
            oCornerData.v = i;
            oCornerData.R = CornerPixel(i,j).dR;
            vector::push_back(oCornerData);
        }
    }
}



int SHGaussianPiramid::Make(const SHImage& imgOrigin, double dSigma, int nOctave)
{
    this->dSigma = dSigma;
    this->nOctave = nOctave;

    // 마스크 만들기
    int nHalf = (int)(dSigma*2); // 95%
    int nFull = 2 * nHalf + 1;

    //최상위 레벨에서 커널 크기의 4배이상이 되도록
    if(nOctave == 0)
    {
        double dDim = _MIN(imgOrigin.Col(),imgOrigin.Row());
        nOctave     = (int)(log(dDim/(double)(nFull))/log(2.0));
    }

    SHImage imgBase(imgOrigin);
    vector::push_back(imgBase);

    for(int i=1; i<nOctave; i++)
    {
        SHImage imgTmp(vector::back());
        imgTmp.GaussianSmoothing(dSigma);
        SHImage imgPush;
        imgTmp.ToHalfSize(imgPush);

        vector::push_back(imgPush);
    }

    return nOctave;
}

SHOpticalData& SHOpticalData::operator = (int n)
{
    dDx = n;
    dDy = n;
    dDxDy = n;
    dDxDx = n;
    dDxDy = n;

    dDt = n;
    dDtDx = n;
    dDtDy = n;

    u = n;
    v = n;

    return *this;
}

void SHOpticalFlow::findFlow(SHImage& imgT0, SHImage& imgT1, int nWindows)
{
    nCol = imgT0.Col();
    nRow = imgT0.Row();

    nHalf = nWindows / 2;
    nFull = nHalf * 2+1;
    double dSigma = nFull / 2.0;

    int nLegnth = nCol * nRow;
    SHArray::Create(nLegnth);

    int i,j,k,l,x,y;


    // Weight Matrix (before diagonal)
    // Using Diagonal => mWeight[0~nWindows^2-1]
    SHMatrix<double> mWeight(nFull,nFull);

    double dScale = 0.0;
    double dSigma2 = Square(dSigma);

    double dConst = 1.0 / PI * 2 / dSigma2;

    y = -nHalf;
    for(i=0;i<nFull;i++)    // y
    {
        x = -nHalf;
        for(j=0;j<nFull;j++)    // x
        {
            mWeight.at(i,j) = dConst*
                    exp((Square<int>(x) + Square<int>(y))
                        * -1.0 / 2.0 / dSigma2);
            dScale += mWeight.at(i,j);
            x++;
        }
        y++;
    }
    mWeight /= dScale;

    // Diff
    for(i=1;i<nRow - 1;i++)
    {
        for(j=1; j<nCol-1; j++)
        {
            pixel(j,i).dDx = imgT1.pixel(j+1,i) - imgT1.pixel(j-1,i);
            pixel(j,i).dDy = imgT1.pixel(j,i+1) - imgT1.pixel(j,i-1);
            pixel(j,i).dDt = imgT1.pixel(j,i) - imgT0.pixel(j,i);
            pixel(j,i).dDxDx = pixel(j,i).dDx * pixel(j,i).dDx;
            pixel(j,i).dDxDy = pixel(j,i).dDx * pixel(j,i).dDy;
            pixel(j,i).dDyDy = pixel(j,i).dDy * pixel(j,i).dDy;
            pixel(j,i).dDtDx = pixel(j,i).dDt * pixel(j,i).dDx;
            pixel(j,i).dDtDy = pixel(j,i).dDt * pixel(j,i).dDy;
            //qDebug("%.3lf",pixel(j,i).dDt);
        }
    }

    /*
    for(i=0;i<nFull;i++)    // y
    {
        for(j=0;j<nFull;j++)    // x
        {
            qDebug("%d %d %.3lf",j,i,mWeight.at(i,j));
    }}
    */

    SHMatrix<double> mATWA(2,2);
    SHMatrix<double> mATWAInv(2,2);
    SHMatrix<double> mATWb(2,1); // row = 2 col = 1

    int nStartX = nHalf + 1; // MaskSize + Diff (1)
    int nStartY = nHalf + 1;
    int nEndX = nCol - nHalf - 1;
    int nEndY = nRow - nHalf - 1;

    for(i=nStartY;i<nEndY;i++)
    {
        for(j=nStartX; j<nEndX; j++)
        {
            // Calculate AWA, AWb
            mATWA.at(0,0) = 0.0;
            mATWA.at(0,1) = 0.0;
            mATWA.at(1,0) = 0.0;
            mATWA.at(1,1) = 0.0;

            mATWb.at(0,0) = 0.0;
            mATWb.at(1,0) = 0.0;
            y = -nHalf;
            for(k=0;k<nFull;k++)
            {
                x = -nHalf;
                for(l=0;l<nFull;l++)
                {
                    mATWA.at(0,0) += pixel(j+x,i+y).dDxDx * mWeight.at(l,k);
                    mATWA.at(0,1) += pixel(j+x,i+y).dDxDy * mWeight.at(l,k);
                    mATWA.at(1,1) += pixel(j+x,i+y).dDyDy * mWeight.at(l,k);

                    mATWb.at(0,0) += pixel(j+x,i+y).dDtDx * mWeight.at(l,k);
                    mATWb.at(1,0) += pixel(j+x,i+y).dDtDy * mWeight.at(l,k);
                    //qDebug("%.3lf %.3lf %.3lf",mATWA.at(0,0),mATWA.at(0,1),mATWA.at(1,1));
                    x++;
                }
                mATWA.at(1,0) = mATWA.at(0,1);
                y++;
            }
            //qDebug("%.3lf %.3lf %.3lf",mATWA.at(0,0),mATWA.at(0,1),mATWA.at(1,1));

            // Calculate inverse AWA
            double dDetATWA = mATWA.at(0,0) * mATWA.at(1,1) - Square(mATWA.at(0,1)); // ad-bc
            if(abs(dDetATWA) < 0.0001) // a inverse matrix is not exist
            {
                pixel(j,i).u = 0.0;
                pixel(j,i).v = 0.0;
                continue;
            }
            mATWAInv.at(0,0) = mATWA.at(1,1);
            mATWAInv.at(1,1) = mATWA.at(0,0);
            mATWAInv.at(0,1) = -mATWA.at(0,1);
            mATWAInv.at(1,0) = -mATWA.at(1,0);
            mATWAInv /= dDetATWA;

            // calculate d vector
            pixel(j,i).u =  mATWAInv.at(0,0) * mATWb.at(0,0) +
                            mATWAInv.at(0,1) * mATWb.at(1,0);
            pixel(j,i).v =  mATWAInv.at(1,0) * mATWb.at(0,0) +
                            mATWAInv.at(1,1) * mATWb.at(1,0);

        }
    }

}

void SHOpticalFlow::findFlowPyramid(SHImage& imgT0, SHImage& imgT1, int nWindows)
{

    nCol = imgT0.Col();
    nRow = imgT0.Row();

    nHalf = nWindows / 2;
    nFull = nHalf * 2+1;
    double dSigma = nFull / 2.0;

    int nLegnth = nCol * nRow;
    SHArray::Create(nLegnth);

    SHGaussianPiramid gpT0(imgT0,dSigma);
    SHGaussianPiramid gpT1(imgT1,dSigma);

    int nPyramidSize = gpT0.size();

    SHOpticalFlow ofProcessor[gpT0.size()];

    int i,x,y;
    ofProcessor[nPyramidSize-1].findFlow(gpT0[nPyramidSize-1], gpT1[nPyramidSize-1], nWindows); // Top level
    for(i=2;i<=nPyramidSize;i++)    // [nPyramidSize - i] -> [nPyramidSize - 1] ... [nPyramidSize - 2] ... [nPyramidSize - 3] ... [0]
    {
        ofProcessor[nPyramidSize-i].findFlow(gpT0[nPyramidSize-i], gpT1[nPyramidSize-i], nWindows); // Top level
        for(y=0; y<ofProcessor[nPyramidSize-i+1].Row(); y++)
        {
            for(x=0; x<ofProcessor[nPyramidSize-i+1].Col(); x++)
            {
                ofProcessor[nPyramidSize-i].pixel(x*2  ,y*2  ).u += ofProcessor[nPyramidSize-i+1].pixel(x,y).u * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2+1,y*2  ).u += ofProcessor[nPyramidSize-i+1].pixel(x,y).u * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2  ,y*2+1).u += ofProcessor[nPyramidSize-i+1].pixel(x,y).u * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2+1,y*2+1).u += ofProcessor[nPyramidSize-i+1].pixel(x,y).u * 2;

                ofProcessor[nPyramidSize-i].pixel(x*2  ,y*2  ).v += ofProcessor[nPyramidSize-i+1].pixel(x,y).v * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2+1,y*2  ).v += ofProcessor[nPyramidSize-i+1].pixel(x,y).v * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2  ,y*2+1).v += ofProcessor[nPyramidSize-i+1].pixel(x,y).v * 2;
                ofProcessor[nPyramidSize-i].pixel(x*2+1,y*2+1).v += ofProcessor[nPyramidSize-i+1].pixel(x,y).v * 2;
            }
        }
    }
    for(y=0;y<nRow;y++)
    {
        for(x=0;x<nCol;x++)
        {
            pixel(x,y).u = ofProcessor[0].pixel(x,y).u;
            pixel(x,y).v = ofProcessor[0].pixel(x,y).v;
        }
    }
}
