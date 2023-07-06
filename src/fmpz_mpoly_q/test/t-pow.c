
/*
  Copyright (C) 2023 David Einstein

  This file is part of Calcium.

  Calcium is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "fmpz_mpoly_q.h"

void fmpz_mpoly_q_pow_naive(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t f,
                          slong n, const fmpz_mpoly_ctx_t ctx)
{
  if (n == 0) {
    fmpz_mpoly_q_one(res, ctx);
  }
  else if (fmpz_mpoly_q_is_zero(f, ctx)){
    fmpz_mpoly_q_zero(res, ctx);
  }
  else if (n == 1) {
    fmpz_mpoly_q_set(res, f, ctx);
  }
  else if (n == -1) {
    fmpz_mpoly_q_inv(res, f, ctx);
  }
  else
    {
      slong n1 = (n < 0 ? -n : n);
      slong i;
      fmpz_mpoly_q_t pow;

      fmpz_mpoly_q_init(pow, ctx);
      fmpz_mpoly_q_set(pow, f, ctx);

      for (i = 1; i < n1 - 1; i++)
        fmpz_mpoly_q_mul(pow, pow, f, ctx);

      fmpz_mpoly_q_mul(res, pow, f, ctx);
      if (n < 0)
        fmpz_mpoly_q_inv(res, res, ctx);
    }
}

int main(void)
{
    slong i, j, tmul = 20;
    FLINT_TEST_INIT(state);

    flint_printf("fmpz_mpoly_q pow_ui....");
    fflush(stdout);

    /* Check pow_ui against pow_naive */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t f, g, h;
        slong len, len1;
        slong pow;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_q_init(f, ctx);
        fmpz_mpoly_q_init(g, ctx);
        fmpz_mpoly_q_init(h, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);

        pow = n_randint(state, 1 + 50/(len1 + 2));

        exp_bound = n_randint(state, 600) + 2;
        exp_bound1 = n_randint(state, 600) + 10;
        exp_bound1 = n_randint(state, exp_bound1) + 2; /* increase chances of lower values */

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_q_randtest(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_q_randtest(g, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_q_randtest(h, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_q_pow_ui(g, f, pow, ctx);
            if (!fmpz_mpoly_q_is_canonical(g, ctx))
              flint_throw(FLINT_ERROR, "Non canonical result");
            fmpz_mpoly_q_pow_naive(h, f, pow, ctx);
            if (!fmpz_mpoly_q_is_canonical(h, ctx))
              flint_throw(FLINT_ERROR, "Non canonical result");

            if (!fmpz_mpoly_q_equal(g, h, ctx))
            {
                flint_printf("FAIL: Check pow_ui against pow_naive\n");
                flint_printf("i = %d, j = %d\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_q_clear(f, ctx);
        fmpz_mpoly_q_clear(g, ctx);
        fmpz_mpoly_q_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    /* Check pow_si against pow_naive */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t f, g, h;
        slong len, len1;
        slong pow;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_q_init(f, ctx);
        fmpz_mpoly_q_init(g, ctx);
        fmpz_mpoly_q_init(h, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);

        pow = n_randint(state, 1 + 50/(len1 + 2));
        pow -= 25 / (len1 + 2); 

        exp_bound = n_randint(state, 600) + 2;
        exp_bound1 = n_randint(state, 600) + 10;
        exp_bound1 = n_randint(state, exp_bound1) + 2; /* increase chances of lower values */

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_q_randtest(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_q_randtest(g, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_q_randtest(h, state, len, coeff_bits, exp_bound, ctx);

            if (pow < 0 && fmpz_mpoly_q_is_zero(f, ctx)){
              if (fmpz_mpoly_q_pow_si(g, f, pow, ctx) != 0) {
                flint_printf("FAIL: Check pow_si negative exponent of 0 does not fail\n");
              }
            }
            else {
              fmpz_mpoly_q_pow_si(g, f, pow, ctx);
              if (!fmpz_mpoly_q_is_canonical(g, ctx)) {
                flint_throw(FLINT_ERROR, "Non canonical result (pow si)");
              }
              fmpz_mpoly_q_pow_naive(h, f, pow, ctx);
              if (!fmpz_mpoly_q_is_canonical(h, ctx)){
                flint_throw(FLINT_ERROR, "Non canonical result (naive)");
              }
              if (!fmpz_mpoly_q_equal(g, h, ctx))
                {
                  flint_printf("FAIL: Check pow_si against pow_naive\n");
                  flint_printf("i = %d, j = %d\n", i, j);
                  flint_printf("pow = %d\n", pow);
                  fflush(stdout);
                  flint_abort();
                }
            }
        }

        fmpz_mpoly_q_clear(f, ctx);
        fmpz_mpoly_q_clear(g, ctx);
        fmpz_mpoly_q_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    /* Check pow_fmpz against pow_naive */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t f, g, h;
        slong len, len1;
        slong pow;
        fmpz_t pow1;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_q_init(f, ctx);
        fmpz_mpoly_q_init(g, ctx);
        fmpz_mpoly_q_init(h, ctx);
        fmpz_init(pow1);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);

        pow = n_randint(state, 1 + 50/(len1 + 2));
        pow -= 25 / (len1 + 2);
        fmpz_set_si(pow1, pow);

        exp_bound = n_randint(state, 600) + 2;
        exp_bound1 = n_randint(state, 600) + 10;
        exp_bound1 = n_randint(state, exp_bound1) + 2; /* increase chances of lower values */

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_q_randtest(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_q_randtest(g, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_q_randtest(h, state, len, coeff_bits, exp_bound, ctx);

            if (pow < 0 && fmpz_mpoly_q_is_zero(f, ctx)){
              if (fmpz_mpoly_q_pow_fmpz(g, f, pow1, ctx) != 0) {
                flint_printf("FAIL: Check pow_fmpz negative exponent of 0 does not fail\n");
              }
            }
            else {
              fmpz_mpoly_q_pow_fmpz(g, f, pow1, ctx);
              if (!fmpz_mpoly_q_is_canonical(g, ctx)) {
                flint_throw(FLINT_ERROR, "Non canonical result (pow si)");
              }
              fmpz_mpoly_q_pow_naive(h, f, pow, ctx);
              if (!fmpz_mpoly_q_is_canonical(h, ctx)){
                flint_throw(FLINT_ERROR, "Non canonical result (naive)");
              }
              if (!fmpz_mpoly_q_equal(g, h, ctx))
                {
                  flint_printf("FAIL: Check pow_fmpz against pow_naive\n");
                  flint_printf("i = %d, j = %d\n", i, j);
                  flint_printf("pow = %d\n", pow);
                  fflush(stdout);
                  flint_abort();
                }
            }
        }

        fmpz_mpoly_q_clear(f, ctx);
        fmpz_mpoly_q_clear(g, ctx);
        fmpz_mpoly_q_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
        fmpz_clear(pow1);
    }
}
