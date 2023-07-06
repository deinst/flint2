/*
  Copyright (C) 2023 David Einstein

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int fmpz_mpoly_q_pow_ui(fmpz_mpoly_q_t A, const fmpz_mpoly_q_t B,
                      ulong k, const fmpz_mpoly_ctx_t ctx)
{
  int status = 1;
  status = fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(A), fmpz_mpoly_q_numref(B),
                             k, ctx);
  if (status) {
    status = fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(A), fmpz_mpoly_q_denref(B),
                               k, ctx);
  }
  return status;
}
