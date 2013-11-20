/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("mullow... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with truncated product of a and b */
    for (i = 0; i < 2000; i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, d;
        long n;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_init)(d, ctx);

        TEMPLATE(T, poly_randtest)(a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest)(b, state, n_randint(state, 100), ctx);
        n = n_randint(state, 100);

        TEMPLATE(T, poly_mullow)(c, a, b, n, ctx);
        TEMPLATE(T, poly_mul)(d, a, b, ctx);
        TEMPLATE(T, poly_truncate)(d, n, ctx);

        result = (TEMPLATE(T, poly_equal)(c, d, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty)(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty)(b, "X", ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty)(c, "X", ctx), flint_printf("\n");
            flint_printf("d = "), TEMPLATE(T, poly_print_pretty)(d, "X", ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);
        TEMPLATE(T, poly_clear)(d, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif