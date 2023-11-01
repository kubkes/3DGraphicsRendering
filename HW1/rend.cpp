#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"

// iostream
#include <iostream>
#include <bitset>


GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */
	 // set display resolution
	 // check boundary constants
	 // typecast int to unsigned short
	if (xRes >= MAXXRES)
	{
		xres = (unsigned short)MAXXRES;
		std::cout << "x resolution set to max" << std::endl;
	}
	else
		xres = (unsigned short)xRes;

	if (yRes >= MAXYRES)
	{
		yres = (unsigned short)MAXYRES;
		std::cout << "y resolution set to max" << std::endl;
	}
	else
		yres = (unsigned short)yRes;

	// frame buffer's size 
	int sizeFB = 3 * xres * yres * sizeof(char);
	// allocate memory for frame buffer
	framebuffer = new char[sizeFB];

	// pixel buffer's size
	int sizePB = xres * yres;
	// allocate memory for pixel buffer
	pixelbuffer = new GzPixel[sizePB];
}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */

		// free frame buffer memory
	delete[] framebuffer;

	// free pixel buffer memory
	delete[] pixelbuffer;
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */

	// RGB values are for a light blue color (8-bit version: 174,230,234)
	GzIntensity r = 174 * 16;
	GzIntensity g = 230 * 16;
	GzIntensity b = 234 * 16;

	GzPixel tempPixel = {
		r,		// red
		g,		// green
		b,		// blue
		1,		// alpha
		0		// z
	};

	for (int i = 0; i < xres; i++)
	{
		for (int j = 0; j < yres; j++)
		{
			pixelbuffer[ARRAY(i, j)] = tempPixel;
			framebuffer[ARRAY(i, j) * 3] = (char)b;
			framebuffer[ARRAY(i, j) * 3 + 1] = (char)g;
			framebuffer[ARRAY(i, j) * 3 + 2] = (char)r;
		}
	}

	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */

	GzPixel tempPixel = {
		r,		// red
		g,		// green
		b,		// blue
		a,		// alpha
		z		// z
	};

	// check for pixels outside of boundaries
	if (i < xres && i >= 0 && j < yres && j >= 0)
		pixelbuffer[ARRAY(i, j)] = tempPixel;

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

		// retrieve the GzPixel at pixel (i,j)
	GzPixel tempPixel = pixelbuffer[ARRAY(i, j)];

	// pass the values inside the tempPixel to the r, g, b, a, z pointers
	*r = tempPixel.red;
	*g = tempPixel.green;
	*b = tempPixel.blue;
	*a = tempPixel.alpha;
	*z = tempPixel.z;

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	// header info for ppm output
	fprintf(outfile, "P6 %d %d 255\r", xres, yres);

	// temporary red, green, blue
	GzIntensity r;
	GzIntensity g;
	GzIntensity b;

	// maximum 12-bit value
	int MAX12B = 4095;

	// write each pixel's RGB to the ppm file
	// write row by row so the outer loop is over y values
	for (int i = 0; i < yres; i++)
	{
		for (int j = 0; j < xres; j++)
		{
			// read values into the temporary rgb
			r = pixelbuffer[ARRAY(j, i)].red;
			g = pixelbuffer[ARRAY(j, i)].green;
			b = pixelbuffer[ARRAY(j, i)].blue;

			// clamp values to 0-4095
			if (r > MAX12B)
				r = MAX12B;
			if (r < 0)
				r = 0;
			if (g > MAX12B)
				g = MAX12B;
			if (g < 0)
				g = 0;
			if (b > MAX12B)
				b = MAX12B;
			if (b < 0)
				b = 0;

			// turn 12-bit values to 8-bit chars by shifting
			char* rgb8bit = new char[3];
			rgb8bit[0] = (r >> 4);
			rgb8bit[1] = (g >> 4);
			rgb8bit[2] = (b >> 4);

			// write clamped values to the file
			fwrite(rgb8bit, sizeof(char), 3, outfile);
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/

	// temporary red, green, blue
	GzIntensity r;
	GzIntensity g;
	GzIntensity b;

	// maximum 12-bit value
	int MAX12B = 4095;

	// write each pixel's RGB to the frame buffer file
	for (int i = 0; i < xres; i++)
	{
		for (int j = 0; j < yres; j++)
		{
			// read values into the temporary rgb
			r = pixelbuffer[ARRAY(i, j)].red;
			g = pixelbuffer[ARRAY(i, j)].green;
			b = pixelbuffer[ARRAY(i, j)].blue;

			// clamp values to 0-4095
			if (r > MAX12B)
				r = MAX12B;
			if (r < 0)
				r = 0;
			if (g > MAX12B)
				g = MAX12B;
			if (g < 0)
				g = 0;
			if (b > MAX12B)
				b = MAX12B;
			if (b < 0)
				b = 0;

			// turn 12-bit values to 8-bit chars by shifting
			char r8bit = (char)(r >> 4);
			char g8bit = (char)(g >> 4);
			char b8bit = (char)(b >> 4);

			// write clamped values to the frame buffer (BGR)
			framebuffer[ARRAY(i, j) * 3] = b8bit;
			framebuffer[ARRAY(i, j) * 3 + 1] = g8bit;
			framebuffer[ARRAY(i, j) * 3 + 2] = r8bit;
		}
	}
	return GZ_SUCCESS;
}