/******************************************************************************
 *
 * Project:  GRIB Driver (using grib_api)
 * Purpose:  GDALDataset driver for GRIB translator for read support
 * Author:   Alberto Valverde, <alberto@meteogrid.com>
 *
 ******************************************************************************
 * Copyright (c) 2016, Alberto Valverde
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************
 *
 */

#include "cpl_multiproc.h"
#include "gdal_frmts.h"
#include "gdal_pam.h"
#include "ogr_spatialref.h"

#include <algorithm>

CPL_CVSID("$Id$");

static CPLMutex *hGRIB2Mutex = NULL;

/************************************************************************/
/* ==================================================================== */
/*                              GRIB2Dataset                             */
/* ==================================================================== */
/************************************************************************/

class GRIB2RasterBand;

class GRIB2Dataset : public GDALPamDataset
{
    friend class GRIB2RasterBand;

  public:
                GRIB2Dataset();
                ~GRIB2Dataset();

    static GDALDataset *Open( GDALOpenInfo * );
    static int          Identify( GDALOpenInfo * );

    CPLErr      GetGeoTransform( double * padfTransform );
    const char *GetProjectionRef();

  private:
    void SetGribMetaData(grib_MetaData* meta);
    VSILFILE    *fp;
    char  *pszProjection;
    // Calculate and store once as GetGeoTransform may be called multiple times.
    double adfGeoTransform[6];

    GIntBig  nCachedBytes;
    GIntBig  nCachedBytesThreshold;
    int      bCacheOnlyOneBand;
    GRIB2RasterBand* poLastUsedBand;
};

/************************************************************************/
/* ==================================================================== */
/*                            GRIB2RasterBand                             */
/* ==================================================================== */
/************************************************************************/

class GRIB2RasterBand : public GDALPamRasterBand
{
    friend class GRIB2Dataset;

public:
    GRIB2RasterBand( GRIB2Dataset*, int, inventoryType* );
    virtual ~GRIB2RasterBand();
    virtual CPLErr IReadBlock( int, int, void * );
    virtual const char *GetDescription() const;

    virtual double GetNoDataValue( int *pbSuccess = NULL );

    void    FindPDSTemplate();

    void    UncacheData();

private:
    CPLErr       LoadData();

    static void ReadGribData( DataSource &, sInt4, int, double**,
                              grib_MetaData** );
    sInt4 start;
    int subgNum;
    char *longFstLevel;

    double * m_Grib_Data;
    grib_MetaData* m_Grib_MetaData;

    int      nGribDataXSize;
    int      nGribDataYSize;
};

/************************************************************************/
/*                         ConvertUnitInText()                          */
/************************************************************************/

static CPLString ConvertUnitInText( bool bMetricUnits, const char* pszTxt )
{
    if( !bMetricUnits )
        return pszTxt;

    CPLString osRes(pszTxt);
    size_t iPos = osRes.find("[K]");
    if( iPos != std::string::npos )
        osRes = osRes.substr(0, iPos) + "[C]" + osRes.substr(iPos + 3);
    return osRes;
}

/************************************************************************/
/*                           GRIB2RasterBand()                            */
/************************************************************************/

GRIB2RasterBand::GRIB2RasterBand( GRIB2Dataset *poDSIn, int nBandIn,
                                inventoryType *psInv ) :
    start(psInv->start),
    subgNum(psInv->subgNum),
    longFstLevel(CPLStrdup(psInv->longFstLevel)),
    m_Grib_Data(NULL),
    m_Grib_MetaData(NULL),
    nGribDataXSize(poDSIn->nRasterXSize),
    nGribDataYSize(poDSIn->nRasterYSize)
{
    poDS = poDSIn;
    nBand = nBandIn;

    // Let user do -ot Float32 if needed for saving space, GRIB2 contains
    // Float64 (though not fully utilized most of the time).
    eDataType = GDT_Float64;

    nBlockXSize = poDSIn->nRasterXSize;
    nBlockYSize = 1;

    const char* pszGribNormalizeUnits =
        CPLGetConfigOption("GRIB_NORMALIZE_UNITS", "YES");
    bool bMetricUnits = CPLTestBool(pszGribNormalizeUnits);

    SetMetadataItem( "GRIB_UNIT",
                     ConvertUnitInText(bMetricUnits, psInv->unitName) );
    SetMetadataItem( "GRIB_COMMENT",
                     ConvertUnitInText(bMetricUnits, psInv->comment) );
    SetMetadataItem( "GRIB_ELEMENT", psInv->element );
    SetMetadataItem( "GRIB_SHORT_NAME", psInv->shortFstLevel );
    SetMetadataItem( "GRIB_REF_TIME",
                     CPLString().Printf("%12.0f sec UTC", psInv->refTime ) );
    SetMetadataItem( "GRIB_VALID_TIME",
                     CPLString().Printf("%12.0f sec UTC", psInv->validTime ) );
    SetMetadataItem( "GRIB_FORECAST_SECONDS",
                     CPLString().Printf("%.0f sec", psInv->foreSec ) );
}

/************************************************************************/
/*                          FindPDSTemplate()                           */
/*                                                                      */
/*      Scan the file for the PDS template info and represent it as     */
/*      metadata.                                                       */
/************************************************************************/

void GRIB2RasterBand::FindPDSTemplate()

{
    GRIB2Dataset *poGDS = reinterpret_cast<GRIB2Dataset *>( poDS );

/* -------------------------------------------------------------------- */
/*      Collect section 4 octet information ... we read the file        */
/*      ourselves since the GRIB2 API does not appear to preserve all    */
/*      this for us.                                                    */
/* -------------------------------------------------------------------- */
    GIntBig nOffset = VSIFTellL( poGDS->fp );

    VSIFSeekL( poGDS->fp, start+16, SEEK_SET );

    GByte abyHead[5] = { 0 };
    VSIFReadL( abyHead, 5, 1, poGDS->fp );

    GUInt32 nSectSize = 0;
    while( abyHead[4] != 4 )
    {
        memcpy( &nSectSize, abyHead, 4 );
        CPL_MSBPTR32( &nSectSize );

        if( VSIFSeekL( poGDS->fp, nSectSize-5, SEEK_CUR ) != 0
            || VSIFReadL( abyHead, 5, 1, poGDS->fp ) != 1 )
            break;
    }

    if( abyHead[4] == 4 )
    {
        memcpy( &nSectSize, abyHead, 4 );
        CPL_MSBPTR32( &nSectSize );

        GByte *pabyBody = static_cast<GByte *>( CPLMalloc(nSectSize - 5) );
        VSIFReadL( pabyBody, 1, nSectSize-5, poGDS->fp );

        GUInt16 nCoordCount = 0;
        memcpy( &nCoordCount, pabyBody + 5 - 5, 2 );
        CPL_MSBPTR16( &nCoordCount );

        GUInt16 nPDTN = 0;
        memcpy( &nPDTN, pabyBody + 7 - 5, 2 );
        CPL_MSBPTR16( &nPDTN );

        SetMetadataItem( "GRIB_PDS_PDTN",
                         CPLString().Printf( "%d", nPDTN ) );

        CPLString osOctet;
        for( int i = 9; i < static_cast<int>( nSectSize ); i++ )
        {
            char szByte[10] = { '\0' };

            if( i == 9 )
                snprintf( szByte, sizeof(szByte), "%d", pabyBody[i-5] );
            else
                snprintf( szByte, sizeof(szByte), " %d", pabyBody[i-5] );
            osOctet += szByte;
        }

        SetMetadataItem( "GRIB_PDS_TEMPLATE_NUMBERS", osOctet );

        CPLFree( pabyBody );
    }

    VSIFSeekL( poGDS->fp, nOffset, SEEK_SET );
}

/************************************************************************/
/*                         GetDescription()                             */
/************************************************************************/

const char * GRIB2RasterBand::GetDescription() const
{
    if( longFstLevel == NULL )
        return GDALPamRasterBand::GetDescription();

    return longFstLevel;
}

/************************************************************************/
/*                             LoadData()                               */
/************************************************************************/

CPLErr GRIB2RasterBand::LoadData()

{
    if( !m_Grib_Data )
    {
        GRIB2Dataset *poGDS = reinterpret_cast<GRIB2Dataset *>( poDS );

        if (poGDS->bCacheOnlyOneBand)
        {
            // In "one-band-at-a-time" strategy, if the last recently used
            // band is not that one, uncache it. We could use a smarter strategy
            // based on a LRU, but that's a bit overkill for now.
            poGDS->poLastUsedBand->UncacheData();
            poGDS->nCachedBytes = 0;
        }
        else
        {
            // Once we have cached more than nCachedBytesThreshold bytes, we
            // will switch to "one-band-at-a-time" strategy, instead of caching
            // all bands that have been accessed.
            if (poGDS->nCachedBytes > poGDS->nCachedBytesThreshold)
            {
                CPLDebug( "GRIB2",
                          "Maximum band cache size reached for this dataset. "
                          "Caching only one band at a time from now");
                for(int i=0;i<poGDS->nBands;i++)
                {
                    reinterpret_cast<GRIB2RasterBand*>(
                        poGDS->GetRasterBand(i+1))->UncacheData();
                }
                poGDS->nCachedBytes = 0;
                poGDS->bCacheOnlyOneBand = TRUE;
            }
        }

        FileDataSource grib_fp (poGDS->fp);

        // we don't seem to have any way to detect errors in this!
        ReadGribData(grib_fp, start, subgNum, &m_Grib_Data, &m_Grib_MetaData);
        if( !m_Grib_Data )
        {
            CPLError( CE_Failure, CPLE_AppDefined, "Out of memory." );
            return CE_Failure;
        }

/* -------------------------------------------------------------------- */
/*      Check that this band matches the dataset as a whole, size       */
/*      wise. (#3246)                                                   */
/* -------------------------------------------------------------------- */
        nGribDataXSize = m_Grib_MetaData->gds.Nx;
        nGribDataYSize = m_Grib_MetaData->gds.Ny;

        poGDS->nCachedBytes += nGribDataXSize * nGribDataYSize * sizeof(double);
        poGDS->poLastUsedBand = this;

        if( nGribDataXSize != nRasterXSize
            || nGribDataYSize != nRasterYSize )
        {
            CPLError( CE_Warning, CPLE_AppDefined,
                      "Band %d of GRIB2 dataset is %dx%d, while the first band "
                      "and dataset is %dx%d.  Georeferencing of band %d may "
                      "be incorrect, and data access may be incomplete.",
                      nBand,
                      nGribDataXSize, nGribDataYSize,
                      nRasterXSize, nRasterYSize,
                      nBand );
        }
    }

    return CE_None;
}

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr GRIB2RasterBand::IReadBlock( int /* nBlockXOff */,
                                   int nBlockYOff,
                                   void * pImage )

{
    CPLErr eErr = LoadData();
    if (eErr != CE_None)
        return eErr;

/* -------------------------------------------------------------------- */
/*      The image as read is always upside down to our normal           */
/*      orientation so we need to effectively flip it at this           */
/*      point.  We also need to deal with bands that are a different    */
/*      size than the dataset as a whole.                               */
/* -------------------------------------------------------------------- */
    if( nGribDataXSize == nRasterXSize
        && nGribDataYSize == nRasterYSize )
    {
        // Simple 1:1 case.
        memcpy(pImage,
               m_Grib_Data + nRasterXSize * (nRasterYSize - nBlockYOff - 1),
               nRasterXSize * sizeof(double));

        return CE_None;
    }

    memset( pImage, 0, sizeof(double) * nRasterXSize );

    if( nBlockYOff >= nGribDataYSize ) // off image?
        return CE_None;

    const int nCopyWords = std::min(nRasterXSize, nGribDataXSize);

    memcpy( pImage,
            m_Grib_Data + nGribDataXSize*(nGribDataYSize-nBlockYOff-1),
            nCopyWords * sizeof(double) );

    return CE_None;
}

/************************************************************************/
/*                           GetNoDataValue()                           */
/************************************************************************/

double GRIB2RasterBand::GetNoDataValue( int *pbSuccess )
{
    CPLErr eErr = LoadData();
    if (eErr != CE_None ||
        m_Grib_MetaData == NULL ||
        m_Grib_MetaData->gridAttrib.f_miss == 0)
    {
        if (pbSuccess)
            *pbSuccess = FALSE;
        return 0;
    }

    if (m_Grib_MetaData->gridAttrib.f_miss == 2)
    {
        /* what TODO ? */
        CPLDebug("GRIB2", "Secondary missing value also set for band %d : %f",
                 nBand, m_Grib_MetaData->gridAttrib.missSec);
    }

    if (pbSuccess)
        *pbSuccess = TRUE;
    return m_Grib_MetaData->gridAttrib.missPri;
}

/************************************************************************/
/*                            ReadGribData()                            */
/************************************************************************/

void GRIB2RasterBand::ReadGribData( DataSource & fp, sInt4 start, int subgNum,
                                   double** data, grib_MetaData** metaData)
{
}

/************************************************************************/
/*                            UncacheData()                             */
/************************************************************************/

void GRIB2RasterBand::UncacheData()
{
/************************************************************************/
/*                           ~GRIB2RasterBand()                          */
/************************************************************************/

GRIB2RasterBand::~GRIB2RasterBand()
{
}

/************************************************************************/
/* ==================================================================== */
/*                              GRIB2Dataset                             */
/* ==================================================================== */
/************************************************************************/

GRIB2Dataset::GRIB2Dataset() :
    fp(NULL),
    pszProjection(CPLStrdup("")),
    nCachedBytes(0),
    // Switch caching strategy once 100 MB threshold is reached.
    // Why 100 MB ? --> why not.
    nCachedBytesThreshold(
        static_cast<GIntBig>(atoi(CPLGetConfigOption("GRIB_CACHEMAX", "100")))
        * 1024 * 1024),
    bCacheOnlyOneBand(FALSE),
    poLastUsedBand(NULL)
{
  adfGeoTransform[0] = 0.0;
  adfGeoTransform[1] = 1.0;
  adfGeoTransform[2] = 0.0;
  adfGeoTransform[3] = 0.0;
  adfGeoTransform[4] = 0.0;
  adfGeoTransform[5] = 1.0;
}

/************************************************************************/
/*                            ~GRIB2Dataset()                             */
/************************************************************************/

GRIB2Dataset::~GRIB2Dataset()

{
    FlushCache();
    CPLFree( pszProjection );
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr GRIB2Dataset::GetGeoTransform( double * padfTransform )

{
    memcpy( padfTransform,  adfGeoTransform, sizeof(double) * 6 );
    return CE_None;
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char *GRIB2Dataset::GetProjectionRef()

{
    return pszProjection;
}

/************************************************************************/
/*                            Identify()                                */
/************************************************************************/

int GRIB2Dataset::Identify( GDALOpenInfo * poOpenInfo )
{
    return FALSE;
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *GRIB2Dataset::Open( GDALOpenInfo * poOpenInfo )

{
    if( !Identify(poOpenInfo) )
        return NULL;

/* -------------------------------------------------------------------- */
/*      Confirm the requested access is supported.                      */
/* -------------------------------------------------------------------- */
    if( poOpenInfo->eAccess == GA_Update )
    {
        CPLError( CE_Failure, CPLE_NotSupported,
                  "The GRIB2 driver does not support update access to existing"
                  " datasets.\n" );
        return NULL;
    }
/* -------------------------------------------------------------------- */
/*      Create a corresponding GDALDataset.                             */
/* -------------------------------------------------------------------- */
    GRIB2Dataset *poDS = new GRIB2Dataset();

/* -------------------------------------------------------------------- */
/*      Create band objects.                                            */
/* -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- */
/*      Initialize any PAM information.                                 */
/* -------------------------------------------------------------------- */
    poDS->SetDescription( poOpenInfo->pszFilename );

    // Release hGRIB2Mutex otherwise we'll deadlock with GDALDataset own
    // hGRIB2Mutex.
    CPLReleaseMutex(hGRIB2Mutex);
    poDS->TryLoadXML();

/* -------------------------------------------------------------------- */
/*      Check for external overviews.                                   */
/* -------------------------------------------------------------------- */
    poDS->oOvManager.Initialize( poDS, poOpenInfo->pszFilename,
                                 poOpenInfo->GetSiblingFiles() );
    CPLAcquireMutex(hGRIB2Mutex, 1000.0);

    return poDS;
}

/************************************************************************/
/*                            SetMetadata()                             */
/************************************************************************/

void GRIB2Dataset::SetGribMetaData(grib_MetaData* meta)
{
/*
    nRasterXSize = meta->gds.Nx;
    nRasterYSize = meta->gds.Ny;

// --------------------------------------------------------------------
//      Image projection.                                              
// --------------------------------------------------------------------
    OGRSpatialReference oSRS;

    switch(meta->gds.projType)
    {
      case GS3_LATLON:
      case GS3_GAUSSIAN_LATLON:
          // No projection, only latlon system (geographic)
          break;
      case GS3_MERCATOR:
        oSRS.SetMercator(meta->gds.meshLat, meta->gds.orientLon,
                         1.0, 0.0, 0.0);
        break;
      case GS3_POLAR:
        oSRS.SetPS(meta->gds.meshLat, meta->gds.orientLon,
                   meta->gds.scaleLat1,
                   0.0, 0.0);
        break;
      case GS3_LAMBERT:
        oSRS.SetLCC(meta->gds.scaleLat1, meta->gds.scaleLat2,
                    meta->gds.meshLat, meta->gds.orientLon,
                    0.0, 0.0); // set projection
        break;

      case GS3_ORTHOGRAPHIC:

        // oSRS.SetOrthographic( 0.0, meta->gds.orientLon,
        //                       meta->gds.lon2, meta->gds.lat2);

        // oSRS.SetGEOS( meta->gds.orientLon, meta->gds.stretchFactor,
        //               meta->gds.lon2, meta->gds.lat2);

        // TODO: Hardcoded for now. How to parse the meta->gds section?
        oSRS.SetGEOS(  0, 35785831, 0, 0 );
        break;
      case GS3_EQUATOR_EQUIDIST:
        break;
      case GS3_AZIMUTH_RANGE:
        break;
    }

// --------------------------------------------------------------------
//      Earth model                                                   
// --------------------------------------------------------------------
    double a = meta->gds.majEarth * 1000.0; // in meters
    double b = meta->gds.minEarth * 1000.0;
    if( a == 0 && b == 0 )
    {
        a = 6377563.396;
        b = 6356256.910;
    }

    if (meta->gds.f_sphere)
    {
        oSRS.SetGeogCS( "Coordinate System imported from GRIB file",
                        NULL,
                        "Sphere",
                        a, 0.0 );
    }
    else
    {
        const double fInv = a / (a - b);
        oSRS.SetGeogCS( "Coordinate System imported from GRIB file",
                        NULL,
                        "Spheroid imported from GRIB file",
                        a, fInv );
    }

    OGRSpatialReference oLL; // construct the "geographic" part of oSRS
    oLL.CopyGeogCSFrom( &oSRS );

    double rMinX = 0.0;
    double rMaxY = 0.0;
    double rPixelSizeX = 0.0;
    double rPixelSizeY = 0.0;
    if (meta->gds.projType == GS3_ORTHOGRAPHIC)
    {
        // This is what should work, but it doesn't .. Dx seems to have an
        // inverse relation with pixel size.
        // rMinX = -meta->gds.Dx * (meta->gds.Nx / 2);
        // rMaxY = meta->gds.Dy * (meta->gds.Ny / 2);
        // Hardcoded for now, assumption: GEOS projection, full disc (like MSG).
        const double geosExtentInMeters = 11137496.552;
        rMinX = -(geosExtentInMeters / 2);
        rMaxY = geosExtentInMeters / 2;
        rPixelSizeX = geosExtentInMeters / meta->gds.Nx;
        rPixelSizeY = geosExtentInMeters / meta->gds.Ny;
    }
    else if( oSRS.IsProjected() )
    {
        // Longitude in degrees, to be transformed to meters (or degrees in
        // case of latlon).
        rMinX = meta->gds.lon1;
        // Latitude in degrees, to be transformed to meters.
        rMaxY = meta->gds.lat1;
        OGRCoordinateTransformation *poTransformLLtoSRS =
            OGRCreateCoordinateTransformation( &(oLL), &(oSRS) );
        // Transform it to meters.
        if( (poTransformLLtoSRS != NULL) &&
            poTransformLLtoSRS->Transform( 1, &rMinX, &rMaxY ))
        {
            if (meta->gds.scan == GRIB2BIT_2) // Y is minY, GDAL wants maxY
            {
                // -1 because we GDAL needs the coordinates of the centre of
                // the pixel.
                rMaxY += (meta->gds.Ny - 1) * meta->gds.Dy;
            }
            rPixelSizeX = meta->gds.Dx;
            rPixelSizeY = meta->gds.Dy;
        }
        else
        {
            rMinX = 0.0;
            rMaxY = 0.0;

            rPixelSizeX = 1.0;
            rPixelSizeY = -1.0;

            oSRS.Clear();

            CPLError( CE_Warning, CPLE_AppDefined,
                      "Unable to perform coordinate transformations, so the "
                      "correct projected geotransform could not be deduced "
                      "from the lat/long control points.  "
                      "Defaulting to ungeoreferenced." );
        }
        delete poTransformLLtoSRS;
    }
    else
    {
        // Longitude in degrees, to be transformed to meters (or degrees in
        // case of latlon).
        rMinX = meta->gds.lon1;
        // Latitude in degrees, to be transformed to meters.
        rMaxY = meta->gds.lat1;

        double rMinY = meta->gds.lat2;
        if (meta->gds.lat2 > rMaxY)
        {
          rMaxY = meta->gds.lat2;
          rMinY = meta->gds.lat1;
        }

        if( meta->gds.Nx == 1 )
          rPixelSizeX = meta->gds.Dx;
        else if (meta->gds.lon1 > meta->gds.lon2)
          rPixelSizeX =
              (360.0 - (meta->gds.lon1 - meta->gds.lon2)) / (meta->gds.Nx - 1);
        else
          rPixelSizeX = (meta->gds.lon2 - meta->gds.lon1) / (meta->gds.Nx - 1);

        if( meta->gds.Ny == 1 )
            rPixelSizeY = meta->gds.Dy;
        else
            rPixelSizeY = (rMaxY - rMinY) / (meta->gds.Ny - 1);

        // Do some sanity checks for cases that can't be handled by the above
        // pixel size corrections. GRIB1 has a minimum precision of 0.001
        // for latitudes and longitudes, so we'll allow a bit higher than that.
        if (rPixelSizeX < 0 || fabs(rPixelSizeX - meta->gds.Dx) > 0.002)
          rPixelSizeX = meta->gds.Dx;

        if (rPixelSizeY < 0 || fabs(rPixelSizeY - meta->gds.Dy) > 0.002)
          rPixelSizeY = meta->gds.Dy;
    }

    // http://gdal.org/gdal_datamodel.html :
    //   we need the top left corner of the top left pixel.
    //   At the moment we have the center of the pixel.
    rMinX-=rPixelSizeX/2;
    rMaxY+=rPixelSizeY/2;

    adfGeoTransform[0] = rMinX;
    adfGeoTransform[3] = rMaxY;
    adfGeoTransform[1] = rPixelSizeX;
    adfGeoTransform[5] = -rPixelSizeY;

    CPLFree( pszProjection );
    pszProjection = NULL;
    oSRS.exportToWkt( &(pszProjection) );
***/
}

/************************************************************************/
/*                       GDALDeregister_GRIB()                          */
/************************************************************************/

static void GDALDeregister_GRIB2(GDALDriver* )
{
    if( hGRIB2Mutex != NULL )
    {
        CPLDestroyMutex(hGRIB2Mutex);
        hGRIB2Mutex = NULL;
    }
}

/************************************************************************/
/*                         GDALRegister_GRIB()                          */
/************************************************************************/

void GDALRegister_GRIB2()

{
    if( GDALGetDriverByName( "GRIB2" ) != NULL )
        return;

    GDALDriver *poDriver = new GDALDriver();

    poDriver->SetDescription( "GRIB2" );
    poDriver->SetMetadataItem( GDAL_DCAP_RASTER, "YES" );
    poDriver->SetMetadataItem( GDAL_DMD_LONGNAME, "GRIdded Binary (.grb)" );
    poDriver->SetMetadataItem( GDAL_DMD_HELPTOPIC, "frmt_grib2.html" );
    poDriver->SetMetadataItem( GDAL_DMD_EXTENSION, "grb" );
    //poDriver->SetMetadataItem( GDAL_DCAP_VIRTUALIO, "YES" );

    poDriver->pfnOpen = GRIB2Dataset::Open;
    poDriver->pfnIdentify = GRIB2Dataset::Identify;
    poDriver->pfnUnloadDriver = GDALDeregister_GRIB;

    GetGDALDriverManager()->RegisterDriver( poDriver );
}
