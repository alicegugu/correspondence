//=============================================================================
// Classes for Simple OpenGL Texture Management
// (c) 2005-2006 by Martin Christen
//=============================================================================

#include "texture.h"
#include <iostream>
using namespace std;


using namespace cwc;

#if (defined(BYTE_ORDER) && BYTE_ORDER==BIG_ENDIAN) || \
	(defined(__BYTE_ORDER) && __BYTE_ORDER==__BIG_ENDIAN) || \
	defined(__BIG_ENDIAN__)
#define CWC_BIG_ENDIAN
#define sENDIAN "Big Endian" 
#else
#define CWC_LITTLE_ENDIAN
#define sENDIAN "Little Endian"
#endif

//-----------------------------------------------------------------------------

namespace cwc
{
   inline bool isPower2(int x)
   {
	   return ( (x > 0) && ((x & (x - 1)) == 0) );
   }
}

//-----------------------------------------------------------------------------

TextureBase::TextureBase()
{
   _bCreated = false;
   _bHasMultitexture = glewGetExtension("GL_ARB_multitexture");
   _bHasRectTexture =  glewGetExtension("GL_ARB_texture_non_power_of_two");
   _bHasFBO = glewGetExtension("GL_EXT_framebuffer_object");
   _bOldRectTextureExt = glewGetExtension("GL_EXT_texture_rectangle");
   _bOldRectTextureArb = glewGetExtension("GL_ARB_texture_rectangle");
   _bOldRectTextureNv = glewGetExtension("GL_NV_texture_rectangle");
}

//-----------------------------------------------------------------------------

TextureBase::~TextureBase()
{
   if (_bCreated)
      glDeleteTextures(1, &_texname);
}

//-----------------------------------------------------------------------------

void TextureBase::bind(int textureUnit)
{
   //glBindTexture(GL_TEXTURE_RECTANGLE_ARB, _texname);
   
   if (_bHasMultitexture)
      glActiveTextureARB(GL_TEXTURE0_ARB+textureUnit);
   
   glBindTexture(GL_TEXTURE_2D, _texname);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

TextureBase* TextureFactory::CreateTextureFromFile(char* sFilename)
{
   TextureBase* pTexture = 0;

   FIBITMAP* pBitmap = ImageLoader(sFilename, 0); 
   if (pBitmap == 0)
   {
      return pTexture;
   }   
   FREE_IMAGE_TYPE ImageType = FreeImage_GetImageType(pBitmap);
      
   /*FIT_UNKNOWN = 0,	// unknown type
	FIT_BITMAP  = 1,	// standard image			: 1-, 4-, 8-, 16-, 24-, 32-bit
	FIT_UINT16	= 2,	// array of unsigned short	: unsigned 16-bit
	FIT_INT16	= 3,	// array of short			: signed 16-bit
	FIT_UINT32	= 4,	// array of unsigned long	: unsigned 32-bit
	FIT_INT32	= 5,	// array of long			: signed 32-bit
	FIT_FLOAT	= 6,	// array of float			: 32-bit IEEE floating point
	FIT_DOUBLE	= 7,	// array of double			: 64-bit IEEE floating point
	FIT_COMPLEX	= 8,	// array of FICOMPLEX		: 2 x 64-bit IEEE floating point
	FIT_RGB16	= 9,	// 48-bit RGB image			: 3 x 16-bit
	FIT_RGBA16	= 10,	// 64-bit RGBA image		: 4 x 16-bit
	FIT_RGBF	= 11,	// 96-bit RGB float image	: 3 x 32-bit IEEE floating point
	FIT_RGBAF	= 12	// 128-bit RGBA float image	: 4 x 32-bit IEEE floating point
   */
   
   unsigned int nBPP = FreeImage_GetBPP(pBitmap);
   // Bits per pixel: 
   
   unsigned int nWidth = FreeImage_GetWidth(pBitmap);
   unsigned int nHeight = FreeImage_GetHeight(pBitmap);   
   unsigned int nPitch = FreeImage_GetPitch(pBitmap);
   
   // Do we need a rectangular texture ?
   bool bRectTexture = true;
   if (isPower2(nHeight) && isPower2(nWidth))
      bRectTexture = false;
   
   if (bRectTexture)
   {
      if (!glewGetExtension("GL_ARB_texture_non_power_of_two"))
      {
         std::cout << "Warning: GL_ARB_texture_non_power_of_two is not available!\n"; 
         FreeImage_Unload(pBitmap);
         return 0;
      }
      
      pTexture = new Texture();
   }
   else
   {
      pTexture = new Texture();
   }
   
   void* bits = (void*)FreeImage_GetBits(pBitmap);
   unsigned char* texturedata = 0;
   
   
   GLenum textureFormat;
   GLenum datatype;
   GLint  components;
      
   if((ImageType == FIT_BITMAP) && (FreeImage_GetBPP(pBitmap) == 24))
   { 
      textureFormat = GL_RGBA;
      datatype = GL_UNSIGNED_BYTE;
      components = GL_RGBA;
      
      texturedata = (unsigned char*)malloc(nWidth*nHeight*4);
      
      int red,green,blue,alpha;
         
#ifdef CWC_BIG_ENDIAN
      red = 0; green = 1; blue = 2; alpha = 3;
#else
      red = 2; green = 1; blue = 0; alpha = 3;
#endif
      
      int offset = 0;
      int offset_img =0;
      
      for (unsigned int y=0;y<nHeight;y++)
      {
         for (unsigned int x=0;x<nWidth;x++)
         {         
            texturedata[offset+red]    = ((unsigned char*)bits)[offset_img+0];
            texturedata[offset+green]  = ((unsigned char*)bits)[offset_img+1];
            texturedata[offset+blue]   = ((unsigned char*)bits)[offset_img+2];
            texturedata[offset+alpha]  = 0;
            offset += 4; 
            offset_img += 3;
         }
         offset_img = y*nPitch;
      }
      
   }
   else if ((ImageType == FIT_BITMAP) && (FreeImage_GetBPP(pBitmap) == 32))
   { 
      textureFormat = GL_RGBA;
      datatype = GL_UNSIGNED_BYTE;
      components = GL_RGBA;
      
      texturedata = (unsigned char*)malloc(nWidth*nHeight*4);
      int red,green,blue,alpha;
            
// untested under Linux... if you get wrong colors you have to change it here.

#ifdef OS_MACOSX            
#ifdef CWC_BIG_ENDIAN
      red = 0; green = 1; blue = 2; alpha = 3;
#else
      red = 2; green = 1; blue = 0; alpha = 3;
#endif
#else
      red = 0; green = 1; blue = 2; alpha = 3;
#endif

      int offset = 0;
      int offset_img;
      
      for (unsigned int y=0;y<nHeight;y++)
      {
         for (unsigned int x=0;x<nWidth;x++)
         {         
            texturedata[offset+red]    = ((unsigned char*)bits)[offset_img+0];
            texturedata[offset+green]  = ((unsigned char*)bits)[offset_img+1];
            texturedata[offset+blue]   = ((unsigned char*)bits)[offset_img+2];
            texturedata[offset+alpha]  = ((unsigned char*)bits)[offset_img+3];
            offset += 4;
            offset_img +=4; 
         }
         offset_img = y*nPitch;
      }  
   }
   else
   {
      std::cout << "Warning: Image Format is currently not supported!\n";
      FreeImage_Unload(pBitmap);
      return 0; 
   }
   
   pTexture->Create(nWidth, nHeight, components, textureFormat, datatype, texturedata);  
   
   FreeImage_Unload(pBitmap);
   
   return pTexture;
}

TextureBase* TextureFactory::CreateCheckBoxTexture(int cbtextureType, int cbwidth, int cbheight)
{
	//allocate the memory of texture
	GLubyte* checkImage = (GLubyte*)malloc(sizeof(GLubyte)*cbheight*cbwidth*4);
	if (!checkImage)
	{
		cerr<<"No memory available for creating checkbox texture"<<endl;
		return NULL;
	}
	
	int color;
	switch(cbtextureType)
	{
	case CheckBoxTextureType::WHITE_BLACK_CHECK_BOX:

		//Set value of texture memory as black and white check box texture
		for (int i = 0; i<cbheight; i++)
		{
			for (int j = 0; j<cbwidth; j++)
			{
				color = (((i&0x8)==0)^((j&0x8)==0))*255;
				checkImage[i*4*cbwidth+j*4+0] = (GLubyte)color;
				checkImage[i*4*cbwidth+j*4+1] = (GLubyte)color;
				checkImage[i*4*cbwidth+j*4+2] = (GLubyte)color;
				checkImage[i*4*cbwidth+j*4+3] = (GLubyte)255;
			}
		}
		break;
		
	case CheckBoxTextureType::COLOR:
		cout<<"No color texture generated."<<endl;
		free(checkImage);
		return NULL;
		break;

	}
	Texture* pTexture = new Texture();
	pTexture->Create(cbwidth, cbheight, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);  
	free(checkImage);
	return (TextureBase*)pTexture;
}

void TextureBase::Create(int nWidth, int nHeight, GLint components, GLshort format, GLshort datatype, void* data)
{
   // if texture was already created, destroy it.
    if (_bCreated)
    glDeleteTextures(1, &_texname);

    glGenTextures(1, &_texname);
    glBindTexture(GL_TEXTURE_2D, _texname);    
    glTexImage2D(GL_TEXTURE_2D, 0, components, nWidth, nHeight,
                    0, format, GL_UNSIGNED_BYTE, data);
      
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	// Linear Filtering
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	// Linear Filtering
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT ); 
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT ); 

   _bCreated = true;
}

//-----------------------------------------------------------------------------

 FIBITMAP* TextureFactory::ImageLoader(const char* lpszPathName, int flag) 
{
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;

	// check the file signature and deduce its format
	// (the second argument is currently not used by FreeImage)
	fif = FreeImage_GetFileType(lpszPathName, 0);
	if(fif == FIF_UNKNOWN) 
	{
		// no signature ?
		// try to guess the file format from the file extension
		fif = FreeImage_GetFIFFromFilename(lpszPathName);
	}
	// check that the plugin has reading capabilities ...
	if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) 
	{
		// ok, let's load the file
		return FreeImage_Load(fif, lpszPathName, flag);
	}
	
	return 0;
}

//-----------------------------------------------------------------------------

bool TextureFactory::ImageWriter(FIBITMAP* dib, const char* lpszPathName, int flag) 
{
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
	BOOL bSuccess = FALSE;

	if(dib) 
	{
		// try to guess the file format from the file extension
		fif = FreeImage_GetFIFFromFilename(lpszPathName);
		if(fif != FIF_UNKNOWN ) 
		{
			// check that the plugin has sufficient writing and export capabilities ...
			WORD bpp = FreeImage_GetBPP(dib);
			if(FreeImage_FIFSupportsWriting(fif) && FreeImage_FIFSupportsExportBPP(fif, bpp)) 
			{
				// ok, we can save the file
				bSuccess = FreeImage_Save(fif, dib, lpszPathName, flag);
				// unless an abnormal bug, we are done !
			}
		}
	}
	return (bSuccess == TRUE) ? true : false;
}

//-----------------------------------------------------------------------------