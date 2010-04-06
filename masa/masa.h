/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PSDNS Development Team
 *
 * Please see https://wiki.ices.utexas.edu/PSDNS for more information.
 *
 * This file is part of the esio library.
 *
 * esio is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * esio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with esio.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * esio.h: Main header for all esio functionality
 *
 * $Id: esio.h 4643 2009-09-04 19:40:42Z nick $
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef _ESIO_H
#define _ESIO_H

/* External header dependencies here, e.g. #include <stdlib.h> */

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

void esio_test_(int *passed);
void esio_init_(char FILE1[], int *ny,int *nx,int *nz,int *nc,int *passed,int fname_len);
void esio_read_();
void esio_read_to_var_();
void esio_write_();
void esio_write_field_(char FILE1[], int *ny,int *xist,int *zjst, int *xisz, int *zjsz, int fname_len);
void esio_read_field_();
void esio_timer_();
void esio_diff_();

__END_DECLS

#endif /* _ESIO_H */
