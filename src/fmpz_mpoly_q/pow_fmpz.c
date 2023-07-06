/*
  Copyright (C) 2023 David Einstein

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int fmpz_mpoly_q_pow_fmpz(fmpz_mpoly_q_t A, const fmpz_mpoly_q_t B,
                          const fmpz_t k, const fmpz_mpoly_ctx_t ctx)
{
  if (fmpz_fits_si(k)) {
    return fmpz_mpoly_q_pow_si(A, B, fmpz_get_si(k), ctx);
  }

  /* the exponent is ridiculous, so we insist that B is either zero or a fraction of monomials with unit coefficients */
  if (fmpz_mpoly_q_is_zero(B,ctx)){
    fmpz_mpoly_q_zero(A, ctx);
    return 1;
  }
  if (B->num.length != 1 || B->den.length != 1) {
    return 0;
  }
  if (!fmpz_is_pm1(B->num.coeffs + 0) || !fmpz_is_pm1(B->den.coeffs + 0)){
    return 0;
  }

  int status=1;
  if (fmpz_sgn(k) < 0){
    fmpz_t temp_k;
    fmpz_init(temp_k);
    status = fmpz_mpoly_pow_fmpz(&(A->den), &(B->num), temp_k, ctx);
    if (status){
      status = fmpz_mpoly_pow_fmpz(&(A->num), &(B->den), temp_k, ctx);
    }
    fmpz_clear(temp_k);
  } else {
    status = fmpz_mpoly_pow_fmpz(&(A->num), &(B->num), k, ctx);
    if (status){
      status = fmpz_mpoly_pow_fmpz(&(A->den), &(B->den), k, ctx);
    }
  }
  return status;
}
