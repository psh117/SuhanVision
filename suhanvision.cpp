/**
  @file suhanvision.cpp
  @author ParkSuhan (psh117@gmail.com)
  @brief Suhan Vision Process Library

  */

#include "suhanvision.h"
#include "cmath"


// 생성자

/**
 * @brief SHImage::SHImage
 */
SHImage::SHImage()
{
    pData = NULL;
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
/**
 * @brief SHImage::~SHImage
 * 소멸자
 */
SHImage::~SHImage()
{
    if(pData != NULL)
        delete pData;
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

    nCol = col;
    nRow = row;
    nChannel = channel;
    nType = type;

    // allocate memory
    int length = nCol * nRow * nChannel;
    pData = new byte[length];

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
            int rgbMax = max<int>(3,r,g,b);
            int rgbMin = min<int>(3,r,g,b);
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
void SHImage::FillData(int channel, byte data)
{
    int length = nCol * nRow;
    for(int i=0; i<length; i++)
    {
        pData[i*nChannel + channel] = data;
    }
}
