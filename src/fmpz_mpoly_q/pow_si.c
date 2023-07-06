/*
  Copyright (C) 2023 David Einstein

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int fmpz_mpoly_q_pow_si(fmpz_mpoly_q_t A, const fmpz_mpoly_q_t B,
                        slong k, const fmpz_mpoly_ctx_t ctx)
{
  if (k < 0) {
    int status = 1;
    if (fmpz_mpoly_q_is_zero(B, ctx)) {
      return 0;
    }
    status = fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(A), fmpz_mpoly_q_numref(B),
                               -k, ctx);
    if (status) {
      status = fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(A), fmpz_mpoly_q_denref(B),
                                 -k, ctx);
      /* The numerator and denominator should be relatively prime, so all we need
         to do is make the denominator positive */
      if (fmpz_sgn(fmpz_mpoly_q_denref(A)->coeffs) < 0){
        fmpz_mpoly_neg(fmpz_mpoly_q_denref(A), fmpz_mpoly_q_denref(A), ctx);
        fmpz_mpoly_neg(fmpz_mpoly_q_numref(A), fmpz_mpoly_q_numref(A), ctx);
      }
    }
    return status;
  } else {
    return fmpz_mpoly_q_pow_ui(A, B, k, ctx);
  }
}
