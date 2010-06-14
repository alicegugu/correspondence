//=============================================================================
// Classes for Simple OpenGL Texture Management
// (c) 2005-2006 by Martin Christen
// edited by Gu Huiqin
//=============================================================================

#ifndef CWC_RENDERTARGET_H
#define CWC_RENDERTARGET_H

#include <GL/glew.h>
#include "FreeImage.h"

//-----------------------------------------------------------------------------
namespace cwc
{

	class TextureBase;
	enum CheckBoxTextureType
	{
		WHITE_BLACK_CHECK_BOX = 0,
		COLOR = 1
	};

	class TextureFactory
	{
	public:
		//! Generate rectangular or power of 2 texture from file:
		static TextureBase* CreateTextureFromFile(char* sFilename);

		//! Generate a texture from memory
		//static TextureBase* CreateTextureFromMemory( ... )

		//! Generate a dynamic texture from a functor
		//static TextureBase* CreateDynamicTexture(DynamicClass);
		static TextureBase* CreateCheckBoxTexture(int cbtextureType, int cbwidth, int cbheight);

	private:   
		static FIBITMAP* ImageLoader(const char* lpszPathName, int flag);
		static bool ImageWriter(FIBITMAP* dib, const char* lpszPathName, int flag);
	};

	//-----------------------------------------------------------------------------
	// Texture Interface

	class TextureBase
	{
	public:
		TextureBase();
		virtual ~TextureBase();

		GLuint   getTexture(void) {return _texname;}
		void     bind(int textureUnit);

		virtual bool   RequiresNormalizedCoords(){return false;}   // returns true if texture coordinates must be normalized

		// "sprite-like" Draw to Screen
		virtual void     drawToScreen(int x, int y, int w, int h) = 0;

		// TODO: RENAME
		//virtual float     xcoord(float nc) = 0;    // convert texture coordinate
		//virtual float     ycoord(float nc) = 0;

		virtual void      Create(int nWidth, int nHeight, GLint components, GLshort format, GLshort datatype, void* data);

		int         GetWidth(){return _width;}
		int         GetHeight(){return _height;}

	protected:
		bool      _bCreated;
		GLuint    _texname;
		int       _width;
		int       _height;
		int       _mode;
		GLboolean _bHasMultitexture;  // true if "GL_ARB_multitexture" extension is available.
		GLboolean _bHasFBO;           // true if "GL_EXT_framebuffer_object" extension is available
		GLboolean _bHasRectTexture;   // true if "GL_ARB_texture_non_power_of_two" extension is available

		GLboolean _bOldRectTextureExt;   // true if "GL_EXT_texture_rectangle" extension is available
		GLboolean _bOldRectTextureArb;   // true if "GL_ARB_texture_rectangle" extension is available
		GLboolean _bOldRectTextureNv;   // true if "GL_NV_texture_rectangle" extension is available
	};

	//-----------------------------------------------------------------------------
	// Powerof2 Texture and rectangular texture using "GL_ARB_texture_non_power_of_two"
	class Texture : public TextureBase
	{
	public:
		virtual void     drawToScreen(int x, int y, int w, int h)
		{
			glBegin(GL_QUADS);
			glTexCoord2f(1.0f,1.0f); glVertex3f(-1.0f, 1.0f, 1.0f);
			glTexCoord2f(1.0f,0.0f); glVertex3f( 1.0f, 1.0f, 1.0f);
			glTexCoord2f(0.0f,0.0f); glVertex3f( 1.0f, 1.0f,-1.0f);
			glTexCoord2f(0.0f,1.0f); glVertex3f(-1.0f, 1.0f,-1.0f);
			glEnd();
		}


	};

	//-----------------------------------------------------------------------------
	// Non Power of 2 Texture using "GL_EXT_texture_rectangle" extension.
	// For the sake of compatibility, don't use if not necessary!

	class RectTextureExt : public TextureBase
	{
	};

	//-----------------------------------------------------------------------------
	// Non Power of 2 Texture using "GL_ARB_texture_rectangle" extension
	// For the sake of compatibility, don't use if not necessary!

	class RectTextureArb : public TextureBase
	{
	};


	//-----------------------------------------------------------------------------
	// Non Power of 2 Texture using "GL_NV_texture_rectangle" extension
	// For the sake of compatibility, don't use if not necessary!

	class RectTextureNv : public TextureBase
	{
	};

	//-----------------------------------------------------------------------------
	// requires: GL_EXT_framebuffer_object

	class RenderTarget : public TextureBase
	{
	public:
		//RenderTarget(int width, int height, int mode=0);
		RenderTarget();
		virtual ~RenderTarget();

		virtual void      Create(int nWidth, int nHeight, GLint components, GLshort format, GLshort datatype, void* data);

		/* void     drawToScreen(int x, int y, int w, int h)
		{
		//if (_tex) _tex->drawToScreen(x,y,w,h);
		}*/

		// returns true if GL_EXT_framebuffer_object is available
		static bool checkFBO(void);

	private:
		void SetRenderTarget();
		void ClearRenderTarget();    // draw to screen again!

		void CheckFramebufferStatus();

		GLuint fb;
		GLuint color_tex;
		GLuint depth_rb;

		int _width, _height, _mode;
		//RectTexture*   _tex;
	};

}
#endif