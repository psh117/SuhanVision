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

// 생성자

/**
 * @brief SHImage::SHImage
 */
SHImage::SHImage()
{
    pData = NULL;
    pHistogramData = NULL;
}

/**
 * @brief SHImage::SHImage
 * @param col 생성할 이미지의 칸 수 (너비)
 * @param row 생성할 이미지의 줄 수 (높이)
 * @param channel 생성할 이미지의 채널 수 (채널)
 * @param type 생성할 이미지의 타입 (RGB, HSV 등)
 */
SHImage::SHImage(int col, int row, int channel, IMG_TYPE type)
{
    pData = NULL;
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
    pData = NULL;
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


SHImage::SHImage(KImageGray& refImage)
{
    pData = NULL;
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
    if(pData != NULL)
        delete pData;
    if(pHistogramData != NULL)
        delete pHistogramData;
}

SHImage& SHImage::operator =(const SHImage& src)
{
    Create(src.Col(),src.Row(),src.Channel(),src.Type());
    nThreshold = src.Threshold();
    int adr = 0;

    for(int i=0;i<nRow;i++)
    {
        for(int j=0;j<nCol;j++)
        {
            for(int k=0;k<nChannel;k++)
            {
                pData[adr] = src[adr];
                adr++;
            }
        }
    }
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
    // delete previous data
    if(pData != NULL)
        delete pData;
    if(pHistogramData != NULL)
        delete pHistogramData;

    nCol = col;
    nRow = row;
    nChannel = channel;
    nType = type;
    nThreshold = -1;

    bHistogram = false;

    // allocate memory
    int length = nCol * nRow * nChannel;
    pData = new int[length];

    // reset all data
    for(int i=0; i<length; i++)
    {
        pData[i] = 0;
    }
}

/**
 * @brief SHImage::Reset
 */
void SHImage::Reset()
{
    if(pData != NULL)
        delete pData;

    pData = NULL;
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
                dstImage[i][j].r = tmp[nPtr * 3 + 0];
                dstImage[i][j].g = tmp[nPtr * 3 + 1];
                dstImage[i][j].b = tmp[nPtr * 3 + 2];
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
            dstImage[i][j] = pData[nPtr];
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

    for(i=0;i<nRow;i++)
    {
        for(j=0;j<nCol;j++)
        {
            int h,s,v;
            int r,g,b;

            int nPtr = i*nCol + j;

            switch (nType)
            {
            case HSV:
                h = pData[nPtr * 3 + 0];
                s = pData[nPtr * 3 + 1];
                v = pData[nPtr * 3 + 2];
                break;
            default:
                break;
            }

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
    int i,j;

    int nHalf = (int)(3.0*dSigma+0.3);
    int nMask[nHalf*2+1][nHalf*2+1];

    double dScale = 0.0;
    double dSigma2 = Square<double>(dSigma);

    double dConst = 1.0 / PI * 2 / dSigma2;



}

void SHImage::MedianFiltering(double dSigma)
{

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
