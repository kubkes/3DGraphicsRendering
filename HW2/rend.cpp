#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

// iostream
#include <iostream>
#include <bitset>

void sortVerticesCW(GzCoord*, GzCoord*);
void calculatePlaneNormal(GzCoord*, float*);

/***********************************************/
/* HW1 methods: copy here the methods from HW1 */

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
		MAXINT	// z
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


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
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


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList) 
{
/* HW 2.1
-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
-- In later homeworks set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++)
	{
		if (nameList[i] == GZ_RGB_COLOR)
		{
			// typecast to GzColor* pointer then access the values inside
			flatcolor[0] = (* (GzColor*) valueList[i])[0]; 
			flatcolor[1] = (* (GzColor*) valueList[i])[1];
			flatcolor[2] = (* (GzColor*) valueList[i])[2];
		}

	}
	
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken *nameList, GzPointer *valueList) 
/* numParts - how many names and values */
{
/* HW 2.2
-- Pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
-- Invoke the rastrizer/scanline framework
-- Return error code
*/
	for (int i = 0; i < numParts; i++)
	{
		if (nameList[i] == GZ_POSITION)
		{
			// Create a vertices array for the triangle
			GzCoord vert[3];
			for (int k = 0; k < 3; k++)
			{
				vert[k][X] = ((GzCoord*)valueList[i])[k][0]; // X
				vert[k][Y] = ((GzCoord*)valueList[i])[k][1]; // Y
				vert[k][Z] = ((GzCoord*)valueList[i])[k][2]; // Z
			}

			// Create an array of sorted vertices
			GzCoord sort_vert[3];

			// Sort the vertices clockwise V0, V1, V2 and write them into sort_vert
			sortVerticesCW(vert, sort_vert);

			// Compute edge parameters A, B, C such that Ax+By+C=0
			float edge1_param[3]; // Edge 1: Tail:V0 Head:V1
			float edge2_param[3]; // Edge 2: Tail:V1 Head:V2
			float edge3_param[3]; // Edge 3: Tail:V2 Head:V0

			edge1_param[0] = sort_vert[1][Y] - sort_vert[0][Y]; // A = dY
			edge1_param[1] = sort_vert[0][X] - sort_vert[1][X]; // B = -dX
			edge1_param[2] = - edge1_param[1] * sort_vert[0][Y] - edge1_param[0] * sort_vert[0][X]; // C = dX*Y-dY*X

			edge2_param[0] = sort_vert[2][Y] - sort_vert[1][Y]; // A = dY
			edge2_param[1] = sort_vert[1][X] - sort_vert[2][X]; // B = -dX
			edge2_param[2] = - edge2_param[1] * sort_vert[1][Y] - edge2_param[0] * sort_vert[1][X]; // C = dX*Y-dY*X

			edge3_param[0] = sort_vert[0][Y] - sort_vert[2][Y]; // A = dY
			edge3_param[1] = sort_vert[2][X] - sort_vert[0][X]; // B = -dX
			edge3_param[2] = - edge3_param[1] * sort_vert[2][Y] - edge3_param[0] * sort_vert[2][X]; // C = dX*Y-dY*X

			// Find the pixels corresponding to the bounding box including the triangle
			int box_minX = (int) floor(min(min(vert[0][X], vert[1][X]), vert[2][X]));
			int box_minY = (int) floor(min(min(vert[0][Y], vert[1][Y]), vert[2][Y]));
			int box_maxX = (int) ceil(max(max(vert[0][X], vert[1][X]), vert[2][X]));
			int box_maxY = (int) ceil(max(max(vert[0][Y], vert[1][Y]), vert[2][Y]));

			// For each pixel in the bounding box
			for (int pixel_x = max(box_minX,0); (pixel_x <= box_maxX && i < xres); pixel_x++)
			{
				for (int pixel_y = max(box_minY,0); (pixel_y <= box_maxY && pixel_y < yres); pixel_y++)
				{
					// Evaluate the edge line equation at pixel i,j
					float edge1_eval = edge1_param[0] * (float)pixel_x + edge1_param[1] * (float)pixel_y + edge1_param[2];
					float edge2_eval = edge2_param[0] * (float)pixel_x + edge2_param[1] * (float)pixel_y + edge2_param[2];
					float edge3_eval = edge3_param[0] * (float)pixel_x + edge3_param[1] * (float)pixel_y + edge3_param[2];

					// Determine whether signes of all evaluations are consistent
					bool signs_consistent = (edge1_eval > 0.0f && edge2_eval > 0.0f && edge3_eval > 0.0f) || 
						(edge1_eval < 0.0f && edge2_eval < 0.0f && edge3_eval < 0.0f);

					// If the pixel is just at the edge3 and other signs are consistent, put it inside the triangle 
					bool on_edge3 = edge3_eval == 0.0f && edge1_eval*edge2_eval > 0.0f;

					// If the pixel is just at the edge2 and other signs are consistent 
					// and if it the triangle is to the right of the pixel, put it inside the right triangle
					// this includes the case where the edge2 is horizontal, in that case pixel is put inside
					// the top triangle
					bool on_edge2 = edge2_eval == 0.0f && edge1_eval*edge3_eval > 0.0f && sort_vert[1][Y] >= sort_vert[2][Y];

					// if the pixel inside the triangle
					if (signs_consistent || on_edge3 || on_edge2)
					{
						float plane_normal[3];

						// Calculate the plane normal parameters (A,B,C)
						calculatePlaneNormal(vert, plane_normal);

						// Solve for D
						float plane_D = -(plane_normal[0] * vert[0][X] + plane_normal[1] * vert[0][Y] + plane_normal[2] * vert[0][Z]);
							
						// Find z value via interpolation
						float z_interp = -(plane_D + plane_normal[0] * (float)pixel_x + plane_normal[1] * (float)pixel_y) / plane_normal[2];

						// Temporary variables
						GzIntensity r;
						GzIntensity g;
						GzIntensity b;
						GzDepth z_buffer;
						GzIntensity a;

						// Get current z value on the buffer
						GzGet(pixel_x, pixel_y, &r, &g, &b, &a, &z_buffer);

						// If z value < z buffer value update the color
						if (z_interp <= (float)z_buffer)
						{
							// Update the color
							r = ctoi(flatcolor[0]);
							g = ctoi(flatcolor[1]);
							b = ctoi(flatcolor[2]);

							// Pass the values to the buffer
							GzPut(pixel_x, pixel_y, r, g, b, a, (int)z_interp);
						}
					}

				}
			}
				
		}
	}

	return GZ_SUCCESS;
}


void sortVerticesCW(GzCoord* vert, GzCoord* sort_vert)
{
	// -- Vertex sorting --
	int vert_TBsort_index[3]; // top to bottom sort list according to Y values
	int vert_sort_index[3]; // sort list according to CW rotation starting from the top (1) then the right (2) and then the left (2)

	// sort Y positions first
	if (vert[0][Y] > vert[1][Y])
	{
		if (vert[1][Y] > vert[2][Y])
		{
			vert_TBsort_index[0] = 2;
			vert_TBsort_index[1] = 1;
			vert_TBsort_index[2] = 0;
		}
		else if (vert[2][Y] > vert[1][Y])
		{
			vert_TBsort_index[0] = 1;

			if (vert[0][Y] > vert[2][Y])
			{
				vert_TBsort_index[1] = 2;
				vert_TBsort_index[2] = 0;
			}
			else if (vert[2][Y] > vert[0][Y])
			{
				vert_TBsort_index[1] = 0;
				vert_TBsort_index[2] = 2;
			}
			else //(vert[2][Y] == vert[0][Y])
			{
				// if equal put the one on the right (higher X value) first in the array
				if (vert[2][X] > vert[0][X])
				{
					vert_TBsort_index[1] = 2;
					vert_TBsort_index[2] = 0;
				}
				else
				{
					vert_TBsort_index[1] = 0;
					vert_TBsort_index[2] = 2;
				}
			}
		}
		else //(vert[2][Y] == vert[1][Y])
		{
			// if equal put the one on the right (higher X value) first in the array
			if (vert[2][X] > vert[1][X])
			{
				vert_TBsort_index[0] = 2;
				vert_TBsort_index[1] = 1;
				vert_TBsort_index[2] = 0;
			}
			else
			{
				vert_TBsort_index[0] = 1;
				vert_TBsort_index[1] = 2;
				vert_TBsort_index[2] = 0;
			}
		}
	}
	else if (vert[1][Y] > vert[0][Y])
	{
		if (vert[0][Y] > vert[2][Y])
		{
			vert_TBsort_index[0] = 2;
			vert_TBsort_index[1] = 0;
			vert_TBsort_index[2] = 1;
		}
		else if (vert[2][Y] > vert[0][Y])
		{
			vert_TBsort_index[0] = 0;

			if (vert[1][Y] > vert[2][Y])
			{
				vert_TBsort_index[1] = 2;
				vert_TBsort_index[2] = 1;
			}
			else if (vert[2][Y] > vert[1][Y])
			{
				vert_TBsort_index[1] = 1;
				vert_TBsort_index[2] = 2;
			}
			else //(vert[2][Y] == vert[1][Y])
			{
				// if equal put the one on the right (higher X value) first in the array
				if (vert[2][X] > vert[1][X])
				{
					vert_TBsort_index[1] = 2;
					vert_TBsort_index[2] = 1;
				}
				else
				{
					vert_TBsort_index[1] = 1;
					vert_TBsort_index[2] = 2;
				}
			}
		}
		else //(vert[2][Y] == vert[0][Y])
		{
			// if equal put the one on the right (higher X value) first in the array
			if (vert[2][X] > vert[0][X])
			{
				vert_TBsort_index[0] = 2;
				vert_TBsort_index[1] = 0;
				vert_TBsort_index[2] = 1;
			}
			else
			{
				vert_TBsort_index[0] = 0;
				vert_TBsort_index[1] = 2;
				vert_TBsort_index[2] = 1;
			}
		}
	}
	else //(vert[0][Y] == vert[1][Y])
	{
		if (vert[2][Y] > vert[0][Y])
		{
			vert_TBsort_index[2] = 2;

			if (vert[0][X] > vert[1][X])
			{
				vert_TBsort_index[0] = 0;
				vert_TBsort_index[1] = 1;
			}
			else
			{
				vert_TBsort_index[0] = 1;
				vert_TBsort_index[1] = 0;
			}
		}
		else
		{
			vert_TBsort_index[0] = 2;

			if (vert[0][X] > vert[1][X])
			{
				vert_TBsort_index[1] = 0;
				vert_TBsort_index[2] = 1;
			}
			else
			{
				vert_TBsort_index[1] = 1;
				vert_TBsort_index[2] = 0;
			}
		}
	}

	// Look at slopes to sort L/R

	// Long and short edge slopes
	float long_edge_slope = (vert[vert_TBsort_index[0]][Y] - vert[vert_TBsort_index[2]][Y]) / (vert[vert_TBsort_index[0]][X] - vert[vert_TBsort_index[2]][X]);
	float short_edge_slope = (vert[vert_TBsort_index[0]][Y] - vert[vert_TBsort_index[1]][Y]) / (vert[vert_TBsort_index[0]][X] - vert[vert_TBsort_index[1]][X]);

	// The bigger the slope the righter it is
	if (long_edge_slope > short_edge_slope)
	{
		vert_sort_index[0] = vert_TBsort_index[0];
		vert_sort_index[1] = vert_TBsort_index[2];
		vert_sort_index[2] = vert_TBsort_index[1];
	}
	else //(short_edge_slope >= long_edge_slope)
	{
		vert_sort_index[0] = vert_TBsort_index[0];
		vert_sort_index[1] = vert_TBsort_index[1];
		vert_sort_index[2] = vert_TBsort_index[2];
	}

	for (int i = 0; i < 3; i++)
	{
		sort_vert[i][X] = vert[vert_sort_index[i]][X];
		sort_vert[i][Y] = vert[vert_sort_index[i]][Y];
		sort_vert[i][Z] = vert[vert_sort_index[i]][Z];
	}

}

void calculatePlaneNormal(GzCoord* vert, float* plane_normal) 
{
	float u[3];
	float v[3];
	
	for (int i = 0; i < 3; i++)
	{
		u[i] = vert[1][i] - vert[0][i];
		v[i] = vert[2][i] - vert[0][i];
	}

	plane_normal[0] = u[1] * v[2] - u[2] * v[1];
	plane_normal[1] = u[2] * v[0] - u[0] * v[2];
	plane_normal[2] = u[0] * v[1] - u[1] * v[0];
}
