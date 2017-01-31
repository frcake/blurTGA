#ifndef TGA_H
#define TGA_H

/******************************************************************************/

typedef struct {
	int w, h;
	unsigned char pixels[1];
} TGA;

/******************************************************************************/

#pragma pack(push, 1)
typedef struct {
	unsigned char id_lenght; /* size of image id */
	unsigned char colormap_type; /* 1 is has a colormap */
	unsigned char image_type; /* compression type */

	short cm_first_entry; /* colormap origin */
	short cm_length; /* colormap length */
	unsigned char cm_size; /* colormap size */

	short x_origin; /* bottom left x coord origin */
	short y_origin; /* bottom left y coord origin */

	short width; /* picture width (in pixels) */
	short height; /* picture height (in pixels) */

	unsigned char pixel_depth; /* bits per pixel: 8, 16, 24 or 32 */
	unsigned char image_descriptor; /* 24 bits = 0x00; 32 bits = 0x80 */
} TGA_Header;
#pragma pack(pop)

/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

TGA *load_tga(const char *filename, TGA_Header *header);
void save_tga(char *filename, TGA_Header header, TGA *image);
void free_tga(TGA *image);

#ifdef __cplusplus
};
#endif

/******************************************************************************/

#endif
