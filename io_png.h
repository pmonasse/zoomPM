/*
 * SPDX-License-Identifier: GPL-3.0-or-later or BSD-2-Clause
 * @file io_png.h
 * @brief PNG input/output
 *
 * Copyright (c) 2010-2011, Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 * Copyright (c) 2021, Pascal Monasse <pascal.monasse@enpc.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef _IO_PNG_H
#define _IO_PNG_H

#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION "0.20210709"

#include <stddef.h>


/* io_png.c */
char *io_png_info(void);
unsigned char *io_png_read_u8(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
unsigned char *io_png_read_u8_rgb(const char *fname, size_t *nxp, size_t *nyp);
unsigned char *io_png_read_u8_gray(const char *fname, size_t *nxp, size_t *nyp);
float *io_png_read_f32(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
float *io_png_read_f32_rgb(const char *fname, size_t *nxp, size_t *nyp);
float *io_png_read_f32_gray(const char *fname, size_t *nxp, size_t *nyp);
int io_png_write_u8(const char *fname, const unsigned char *data, size_t nx, size_t ny, size_t nc);
int io_png_write_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc);

float rgb_to_gray(float r, float g, float b);

#ifdef __cplusplus
}
#endif

#endif /* !_IO_PNG_H */
