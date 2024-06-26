/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i < fac->num; i++)
        fac->exp[i] *= exp;
}
