/******************************************************************************
 * $Id$
 *
 * Project:  KML Driver
 * Purpose:  Class for reading, parsing and handling a kmlfile.
 * Author:   Jens Oberender, j.obi@troja.net
 *
 ******************************************************************************
 * Copyright (c) 2007, Jens Oberender
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
 ****************************************************************************/
#include "kmlnode.h"
#include "kml.h"
#include "cpl_error.h"
#include "cpl_conv.h"
// std
#include <cstdio>
#include <cerrno>
#include <string>
#include <iostream>

KML::KML()
{
	nDepth_ = 0;
	validity = KML_VALIDITY_UNKNOWN;
	pKMLFile_ = NULL;
	sError_ = "";
	poTrunk_ = NULL;
	poCurrent_ = NULL;
	nNumLayers_ = -1;
}

KML::~KML()
{
    if( NULL != pKMLFile_ )
        VSIFCloseL(pKMLFile_);

    delete poTrunk_;
}

bool KML::open(const char * pszFilename)
{
    if( NULL != pKMLFile_ )
        VSIFCloseL( pKMLFile_ );

    pKMLFile_ = VSIFOpenL( pszFilename, "r" );
    if( NULL == pKMLFile_ )
    {
        return FALSE;
    }

    return TRUE;
}

void KML::parse()
{
    std::size_t nDone = 0;
    std::size_t nLen = 0;
    char aBuf[BUFSIZ] = { 0 };

    if( NULL == pKMLFile_ )
    {
        sError_ = "No file given";
        return;
    }

    if(poTrunk_ != NULL) {
        delete poTrunk_;
        poTrunk_ = NULL;
    }

    if(poCurrent_ != NULL)
    {
        delete poCurrent_;
        poCurrent_ = NULL;
    }

    XML_Parser oParser = XML_ParserCreate(NULL);
    XML_SetUserData(oParser, this);
    XML_SetElementHandler(oParser, startElement, endElement);
    XML_SetCharacterDataHandler(oParser, dataHandler);

    do
    {
        nLen = (int)VSIFReadL( aBuf, 1, sizeof(aBuf), pKMLFile_ );
        nDone = VSIFEofL(pKMLFile_);
        if (XML_Parse(oParser, aBuf, nLen, nDone) == XML_STATUS_ERROR)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                        "XML parsing of KML file failed : %s at line %d, column %d",
                        XML_ErrorString(XML_GetErrorCode(oParser)),
                        (int)XML_GetCurrentLineNumber(oParser),
                        (int)XML_GetCurrentColumnNumber(oParser));
            XML_ParserFree(oParser);
            VSIRewindL(pKMLFile_);
            return;
        }
    } while (!nDone && nLen > 0 );

    XML_ParserFree(oParser);
    VSIRewindL(pKMLFile_);
    poCurrent_ = NULL;
}

void KML::checkValidity()
{
    std::size_t nDone = 0;
    std::size_t nLen = 0;
    char aBuf[BUFSIZ] = { 0 };

    if(poTrunk_ != NULL)
    {
        delete poTrunk_;
        poTrunk_ = NULL;
    }

    if(poCurrent_ != NULL)
    {
        delete poCurrent_;
        poCurrent_ = NULL;
    }

    if(pKMLFile_ == NULL)
    {
        this->sError_ = "No file given";
        return;
    }

    XML_Parser oParser = XML_ParserCreate(NULL);
    XML_SetUserData(oParser, this);
    XML_SetElementHandler(oParser, startElementValidate, NULL);
    int nCount = 0;

    /* Parses the file until we find the first element */
    do
    {
        nLen = (int)VSIFReadL( aBuf, 1, sizeof(aBuf), pKMLFile_ );
        nDone = VSIFEofL(pKMLFile_);
        if (XML_Parse(oParser, aBuf, nLen, nDone) == XML_STATUS_ERROR)
        {
            if (nLen <= BUFSIZ-1)
                aBuf[nLen] = 0;
            else
                aBuf[BUFSIZ-1] = 0;
            if (strstr(aBuf, "<?xml") && strstr(aBuf, "<kml"))
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                        "XML parsing of KML file failed : %s at line %d, column %d",
                        XML_ErrorString(XML_GetErrorCode(oParser)),
                        (int)XML_GetCurrentLineNumber(oParser),
                        (int)XML_GetCurrentColumnNumber(oParser));
            }

            XML_ParserFree(oParser);
            VSIRewindL(pKMLFile_);
            return;
        }

        nCount ++;
        /* After reading 50 * BUFSIZE bytes, and not finding whether the file */
        /* is KML or not, we give up and fail silently */
    } while (!nDone && nLen > 0 && validity == KML_VALIDITY_UNKNOWN && nCount < 50);

    XML_ParserFree(oParser);
    VSIRewindL(pKMLFile_);
    poCurrent_ = NULL;
}

void XMLCALL KML::startElement(void* pUserData, const char* pszName, const char** ppszAttr)
{
	int i = 0;
    KMLNode* poMynew = NULL;
	Attribute* poAtt = NULL;
	
	if(((KML *)pUserData)->poTrunk_ == NULL 
        || (((KML *)pUserData)->poCurrent_->getName()).compare("description") != 0)
    {
    	poMynew = new KMLNode();
	    poMynew->setName(pszName);
    	poMynew->setLevel(((KML *)pUserData)->nDepth_);
	
    	for (i = 0; ppszAttr[i]; i += 2)
        {
            poAtt = new Attribute();
            poAtt->sName = ppszAttr[i];
            poAtt->sValue = ppszAttr[i + 1];
            poMynew->addAttribute(poAtt);
    	}

        if(((KML *)pUserData)->poTrunk_ == NULL)
            ((KML *)pUserData)->poTrunk_ = poMynew;
        if(((KML *)pUserData)->poCurrent_ != NULL)
            poMynew->setParent(((KML *)pUserData)->poCurrent_);
        ((KML *)pUserData)->poCurrent_ = poMynew;

    	((KML *)pUserData)->nDepth_++;
    }
    else
    {
        std::string sNewContent = "<";
        sNewContent += pszName;
    	for (i = 0; ppszAttr[i]; i += 2)
        {
            sNewContent += " ";
            sNewContent += ppszAttr[i];
            sNewContent += "=";
            sNewContent += ppszAttr[i + 1];
    	}
        sNewContent += ">";
        ((KML *)pUserData)->poCurrent_->addContent(sNewContent);
    }
}

void XMLCALL KML::startElementValidate(void* pUserData, const char* pszName, const char** ppszAttr)
{
    int i = 0;
    
    if (((KML *)pUserData)->validity != KML_VALIDITY_UNKNOWN)
        return;

    ((KML *)pUserData)->validity = KML_VALIDITY_INVALID;

    if(strcmp(pszName, "kml") == 0)
    {
        // Check all Attributes
        for (i = 0; ppszAttr[i]; i += 2)
        {
            // Find the namespace and determine the KML version
            if(strcmp(ppszAttr[i], "xmlns") == 0)
            {
                // Is it KML 2.2?
                if((strcmp(ppszAttr[i + 1], "http://earth.google.com/kml/2.2") == 0) || 
				   (strcmp(ppszAttr[i + 1], "http://www.opengis.net/kml/2.2") == 0))
                {
                    ((KML *)pUserData)->validity = KML_VALIDITY_VALID;
                    ((KML *)pUserData)->sVersion_ = "2.2";
                }
                else if(strcmp(ppszAttr[i + 1], "http://earth.google.com/kml/2.1") == 0)
                {
                    ((KML *)pUserData)->validity = KML_VALIDITY_VALID;
                    ((KML *)pUserData)->sVersion_ = "2.1";
                }
                else if(strcmp(ppszAttr[i + 1], "http://earth.google.com/kml/2.0") == 0)
                {
                    ((KML *)pUserData)->validity = KML_VALIDITY_VALID;
                    ((KML *)pUserData)->sVersion_ = "2.0";
                }
            }
        }
    }
}

void XMLCALL KML::endElement(void* pUserData, const char* pszName)
{
    KMLNode* poTmp = NULL;
    std::string sData;
    std::string sTmp;
    unsigned short nPos = 0;

	if(((KML *)pUserData)->poCurrent_ != NULL
        && 
	    ((KML *)pUserData)->poCurrent_->getName().compare(pszName) == 0)
    {
    	((KML *)pUserData)->nDepth_--;
        poTmp = ((KML *)pUserData)->poCurrent_;
        // Split the coordinates
        if(((KML *)pUserData)->poCurrent_->getName().compare("coordinates") == 0)
        {
            sData = ((KML *)pUserData)->poCurrent_->getContent(0);
            while(sData.length() > 0)
            {
                // Cut off whitespaces
           	    while((sData[nPos] == ' ' || sData[nPos] == '\n'
                    || sData[nPos] == '\r' || 
        	        sData[nPos] == '\t' ) && sData.length() > 0)
                        sData = sData.substr(1, sData.length()-1);
                
                // Get content
                while(sData[nPos] != ' ' && sData[nPos] != '\n' && sData[nPos] != '\r' && 
                    sData[nPos] != '\t' && nPos < sData.length()) 
                        nPos++;
                sTmp = sData.substr(0, nPos);

            	if(sTmp.length() > 0)
        	        ((KML *)pUserData)->poCurrent_->addContent(sTmp);
                // Cut the content from the rest
        	    if(nPos < sData.length())
                    sData = sData.substr(nPos, sData.length() - nPos);
                else
                    break;
                nPos = 0;
            }
            if(((KML *)pUserData)->poCurrent_->numContent() > 1)
                ((KML *)pUserData)->poCurrent_->deleteContent(0);
        }
        
        if(((KML *)pUserData)->poCurrent_->getParent() != NULL)
            ((KML *)pUserData)->poCurrent_ = ((KML *)pUserData)->poCurrent_->getParent();
        else
            ((KML *)pUserData)->poCurrent_ = NULL;

    	if(!((KML *)pUserData)->isHandled(pszName))
        {
    	    CPLDebug("KML", "Not handled: %s", pszName);
    	    delete poTmp;
	    }
        else
        {
	        if(((KML *)pUserData)->poCurrent_ != NULL)
               	((KML *)pUserData)->poCurrent_->addChildren(poTmp);
	    }
    }
    else if(((KML *)pUserData)->poCurrent_ != NULL)
    {
        std::string sNewContent = "</";
        sNewContent += pszName;
        sNewContent += ">";
        ((KML *)pUserData)->poCurrent_->addContent(sNewContent);
    }
}

void XMLCALL KML::dataHandler(void* pUserData, const char* pszData, int nLen)
{
    if(nLen < 1 || ((KML *)pUserData)->poCurrent_ == NULL)
        return;
    std::string sData(pszData, nLen);
    std::string sTmp;
    
	if(((KML *)pUserData)->poCurrent_->getName().compare("coordinates") == 0)
	{
	    if(((KML *)pUserData)->poCurrent_->numContent() == 0)
	        ((KML *)pUserData)->poCurrent_->addContent(sData);
	    else
	        ((KML *)pUserData)->poCurrent_->appendContent(sData);
    }
    else
    {
    	while(sData[0] == ' ' || sData[0] == '\n' || sData[0] == '\r' || sData[0] == '\t')
        {
    		sData = sData.substr(1, sData.length()-1);
    	}
       	if(sData.length() > 0)
	        ((KML *)pUserData)->poCurrent_->addContent(sData);
	}
}

bool KML::isValid()
{
    checkValidity();

    if( validity == KML_VALIDITY_VALID )
        CPLDebug("KML", "Valid: %d Version: %s", 
                 validity == KML_VALIDITY_VALID, sVersion_.c_str());

    return validity == KML_VALIDITY_VALID;
}

std::string KML::getError() const
{
	return sError_;
}

void KML::classifyNodes()
{
    poTrunk_->classify(this);
}

void KML::eliminateEmpty()
{
    poTrunk_->eliminateEmpty(this);
}

void KML::print(unsigned short nNum)
{
    if( poTrunk_ != NULL )
        poTrunk_->print(nNum);
}

bool KML::isHandled(std::string const& elem) const
{
    if( isLeaf(elem) || isFeature(elem) || isFeatureContainer(elem)
        || isContainer(elem) || isRest(elem) )
    {
        return true;
    }
    return false;
}

bool KML::isLeaf(std::string const& elem) const
{
    return false;
};

bool KML::isFeature(std::string const& elem) const
{
    return false;
};

bool KML::isFeatureContainer(std::string const& elem) const
{
    return false;
};

bool KML::isContainer(std::string const& elem) const
{
    return false;
};

bool KML::isRest(std::string const& elem) const
{
    return false;
};

void KML::findLayers(KMLNode* poNode)
{
    // idle
};

int KML::getNumLayers() const
{
    return nNumLayers_;
}

bool KML::selectLayer(unsigned short nNum) {
    if(this->nNumLayers_ < 1 || (short)nNum >= this->nNumLayers_)
        return FALSE;
    poCurrent_ = poTrunk_->getLayer(nNum);
    if(poCurrent_ == NULL)
        return FALSE;
    else
        return TRUE;
}

std::string KML::getCurrentName() const
{
    std::string tmp;
    if( poCurrent_ != NULL )
    {
        tmp = poCurrent_->getNameElement();
    }
    return tmp;
}

Nodetype KML::getCurrentType() const
{
    if(poCurrent_ != NULL)
        return poCurrent_->getType();
    else
        return Unknown;
}

int KML::getNumFeatures() const
{
    if(poCurrent_ != NULL)
        return static_cast<int>(poCurrent_->getNumFeatures());
    else
        return -1;
}

Feature* KML::getFeature(std::size_t nNum) const
{
    if(poCurrent_ != NULL)
        return poCurrent_->getFeature(nNum);
    else
        return NULL;
}

bool KML::getExtents(double& pdfXMin, double& pdfXMax, double& pdfYMin, double& pdfYMax) const
{
    if( poCurrent_ != NULL )
    {
        Extent const* poXT = poCurrent_->getExtents();
        pdfXMin = poXT->dfX1;
        pdfXMax = poXT->dfX2;
        pdfYMin = poXT->dfY1;
        pdfYMax = poXT->dfY2;

        return true;
    }
    return false;
}

