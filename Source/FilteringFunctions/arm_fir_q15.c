/* ----------------------------------------------------------------------
 * Project:      CMSIS DSP Library
 * Title:        arm_fir_q15.c
 * Description:  Q15 FIR filter processing function
 *
 * $Date:        23 April 2021
 * $Revision:    V1.9.0
 *
 * Target Processor: Cortex-M and Cortex-A cores
 * -------------------------------------------------------------------- */
/*
 * Copyright (C) 2010-2021 ARM Limited or its affiliates. All rights reserved.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "dsp/filtering_functions.h"

/**
  @ingroup groupFilters
 */

/**
  @addtogroup FIR
  @{
 */

/**
  @brief         Processing function for the Q15 FIR filter.
  @param[in]     S          points to an instance of the Q15 FIR filter structure
  @param[in]     pSrc       points to the block of input data
  @param[out]    pDst       points to the block of output data
  @param[in]     blockSize  number of samples to process

  @par           Scaling and Overflow Behavior
                   The function is implemented using a 64-bit internal accumulator.
                   Both coefficients and state variables are represented in 1.15 format and multiplications yield a 2.30 result.
                   The 2.30 intermediate results are accumulated in a 64-bit accumulator in 34.30 format.
                   There is no risk of internal overflow with this approach and the full precision of intermediate multiplications is preserved.
                   After all additions have been performed, the accumulator is truncated to 34.15 format by discarding low 15 bits.
                   Lastly, the accumulator is saturated to yield a result in 1.15 format.

  @remark
                   Refer to \ref arm_fir_fast_q15() for a faster but less precise implementation of this function.
 */
#if defined(ARM_MATH_MVEI) && !defined(ARM_MATH_AUTOVECTORIZE)

#define MVE_ASRL_SAT16(acc, shift) ((sqrshrl_sat48(acc, -(32 - shift)) >> 32) & 0xffffffff)

#define FIR_Q15_CORE(pOutput, nbAcc, nbVecTaps, pSample, vecCoeffs) \
    for (int j = 0; j < nbAcc; j++)                                 \
    {                                                               \
        const q15_t *pSmp = &pSample[j];                            \
        q63_t acc[4];                                               \
                                                                    \
        acc[j] = 0;                                                 \
        for (int i = 0; i < nbVecTaps; i++)                         \
        {                                                           \
            vecIn0 = vld1q(pSmp + 8 * i);                           \
            acc[j] = vmlaldavaq(acc[j], vecIn0, vecCoeffs[i]);      \
        }                                                           \
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc[j], 15);             \
    }

#define FIR_Q15_MAIN_CORE()                                                                  \
    {                                                                                        \
        q15_t *pState = S->pState;         /* State pointer */                               \
        const q15_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */                         \
        q15_t *pStateCur;                  /* Points to the current sample of the state */   \
        const q15_t *pSamples;             /* Temporary pointer to the sample buffer */      \
        q15_t *pOutput;                    /* Temporary pointer to the output buffer */      \
        const q15_t *pTempSrc;             /* Temporary pointer to the source data */        \
        q15_t *pTempDest;                  /* Temporary pointer to the destination buffer */ \
        uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */ \
        int32_t blkCnt;                                                                      \
        q15x8_t vecIn0;                                                                      \
                                                                                             \
        /*                                                                                   \
         * load coefs                                                                        \
         */                                                                                  \
        q15x8_t vecCoeffs[NBVECTAPS];                                                        \
                                                                                             \
        for (int i = 0; i < NBVECTAPS; i++)                                                  \
            vecCoeffs[i] = vldrhq_s16(pCoeffs + 8 * i);                                      \
                                                                                             \
        /*                                                                                   \
         * pState points to state array which contains previous frame (numTaps - 1) samples  \
         * pStateCur points to the location where the new input data should be written       \
         */                                                                                  \
        pStateCur = &(pState[(numTaps - 1u)]);                                               \
        pTempSrc = pSrc;                                                                     \
        pSamples = pState;                                                                   \
        pOutput = pDst;                                                                      \
                                                                                             \
        blkCnt = blockSize >> 2;                                                             \
        while (blkCnt > 0)                                                                   \
        {                                                                                    \
            /*                                                                               \
             * Save 4 input samples in the history buffer                                    \
             */                                                                              \
            vstrhq_s32(pStateCur, vldrhq_s32(pTempSrc));                                     \
            pStateCur += 4;                                                                  \
            pTempSrc += 4;                                                                   \
                                                                                             \
            FIR_Q15_CORE(pOutput, 4, NBVECTAPS, pSamples, vecCoeffs);                        \
            pSamples += 4;                                                                   \
                                                                                             \
            blkCnt--;                                                                        \
        }                                                                                    \
                                                                                             \
        /* tail */                                                                           \
        int32_t residual = blockSize & 3;                                                    \
                                                                                             \
        for (int i = 0; i < residual; i++)                                                   \
            *pStateCur++ = *pTempSrc++;                                                      \
                                                                                             \
        FIR_Q15_CORE(pOutput, residual, NBVECTAPS, pSamples, vecCoeffs);                     \
                                                                                             \
        /*                                                                                   \
         * Copy the samples back into the history buffer start                               \
         */                                                                                  \
        pTempSrc = &pState[blockSize];                                                       \
        pTempDest = pState;                                                                  \
                                                                                             \
        /* current compiler limitation */                                                    \
        blkCnt = (numTaps - 1) >> 3;                                                         \
        while (blkCnt > 0)                                                                   \
        {                                                                                    \
            vstrhq_s16(pTempDest, vldrhq_s16(pTempSrc));                                     \
            pTempSrc += 8;                                                                   \
            pTempDest += 8;                                                                  \
            blkCnt--;                                                                        \
        }                                                                                    \
        blkCnt = (numTaps - 1) & 7;                                                          \
        if (blkCnt > 0)                                                                      \
        {                                                                                    \
            mve_pred16_t p = vctp16q(blkCnt);                                                \
            vstrhq_p_s16(pTempDest, vldrhq_z_s16(pTempSrc, p), p);                           \
        }                                                                                    \
    }

static void arm_fir_q15_25_32_mve(const arm_fir_instance_q15 *S,
                                  const q15_t *__restrict pSrc,
                                  q15_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 32
#define NBVECTAPS (NBTAPS / 8)
    FIR_Q15_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q15_17_24_mve(const arm_fir_instance_q15 *S,
                                  const q15_t *__restrict pSrc,
                                  q15_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 24
#define NBVECTAPS (NBTAPS / 8)
    FIR_Q15_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q15_9_16_mve(const arm_fir_instance_q15 *S,
                                 const q15_t *__restrict pSrc,
                                 q15_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 16
#define NBVECTAPS (NBTAPS / 8)
    FIR_Q15_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q15_1_8_mve(const arm_fir_instance_q15 *S,
                                const q15_t *__restrict pSrc,
                                q15_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 8
#define NBVECTAPS (NBTAPS / 8)
    FIR_Q15_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

ARM_DSP_ATTRIBUTE void arm_fir_q15(
    const arm_fir_instance_q15 *S,
    const q15_t *pSrc,
    q15_t *pDst,
    uint32_t blockSize)
{
    q15_t *pState = S->pState;         /* State pointer */
    const q15_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q15_t *pStateCur;                  /* Points to the current sample of the state */
    const q15_t *pSamples;             /* Temporary pointer to the sample buffer */
    q15_t *pOutput;                    /* Temporary pointer to the output buffer */
    const q15_t *pTempSrc;             /* Temporary pointer to the source data */
    q15_t *pTempDest;                  /* Temporary pointer to the destination buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t blkCnt;
    q15x8_t vecIn0;
    uint32_t tapsBlkCnt = (numTaps + 7) / 8;
    q63_t acc0, acc1, acc2, acc3;

    int32_t nbTaps = (numTaps + 7) >> 3;

    switch (nbTaps)
    {

    case 1:
        arm_fir_q15_1_8_mve(S, pSrc, pDst, blockSize);
        return;
    case 2:
        arm_fir_q15_9_16_mve(S, pSrc, pDst, blockSize);
        return;
    case 3:
        arm_fir_q15_17_24_mve(S, pSrc, pDst, blockSize);
        return;
    case 4:
        arm_fir_q15_25_32_mve(S, pSrc, pDst, blockSize);
        return;
    }
    /*
     * pState points to state array which contains previous frame (numTaps - 1) samples
     * pStateCur points to the location where the new input data should be written
     */
    pStateCur = &(pState[(numTaps - 1u)]);
    pTempSrc = pSrc;
    pSamples = pState;
    pOutput = pDst;
    blkCnt = blockSize >> 2;

    while (blkCnt > 0U)
    {
        const q15_t *pCoeffsTmp = pCoeffs;
        const q15_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;
        acc2 = 0LL;
        acc3 = 0LL;

        /*
         * Save 8 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 8;
        pTempSrc += 8;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            /*
             * load 8 coefs
             */
            q15x8_t vecCoeffs = *(q15x8_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vmlaldavaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[1]);
            acc1 = vmlaldavaq(acc1, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[2]);
            acc2 = vmlaldavaq(acc2, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[3]);
            acc3 = vmlaldavaq(acc3, vecIn0, vecCoeffs);

            pSamplesTmp += 8;
            pCoeffsTmp += 8;
            /*
             * Decrement the taps block loop counter
             */
            i--;
        }

        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc0, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc1, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc2, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc3, 15);

        pSamples += 4;
        /*
         * Decrement the sample block loop counter
         */
        blkCnt--;
    }

    uint32_t residual = blockSize & 3;
    switch (residual)
    {
    case 3:
    {
        const q15_t *pCoeffsTmp = pCoeffs;
        const q15_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;
        acc2 = 0LL;

        /*
         * Save 8 input samples in the history buffer
         */
        *(q15x8_t *)pStateCur = *(q15x8_t *)pTempSrc;
        pStateCur += 8;
        pTempSrc += 8;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            /*
             * load 8 coefs
             */
            q15x8_t vecCoeffs = *(q15x8_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vmlaldavaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[2]);
            acc1 = vmlaldavaq(acc1, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[4]);
            acc2 = vmlaldavaq(acc2, vecIn0, vecCoeffs);

            pSamplesTmp += 8;
            pCoeffsTmp += 8;
            /*
             * Decrement the taps block loop counter
             */
            i--;
        }

        acc0 = asrl(acc0, 15);
        acc1 = asrl(acc1, 15);
        acc2 = asrl(acc2, 15);

        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc0, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc1, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc2, 15);
    }
    break;

    case 2:
    {
        const q15_t *pCoeffsTmp = pCoeffs;
        const q15_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;
        /*
         * Save 8 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 8;
        pTempSrc += 8;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            /*
             * load 8 coefs
             */
            q15x8_t vecCoeffs = *(q15x8_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vmlaldavaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[2]);
            acc1 = vmlaldavaq(acc1, vecIn0, vecCoeffs);

            pSamplesTmp += 8;
            pCoeffsTmp += 8;
            /*
             * Decrement the taps block loop counter
             */
            i--;
        }

        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc0, 15);
        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc1, 15);
    }
    break;

    case 1:
    {
        const q15_t *pCoeffsTmp = pCoeffs;
        const q15_t *pSamplesTmp = pSamples;

        acc0 = 0LL;

        /*
         * Save 8 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 8;
        pTempSrc += 8;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            /*
             * load 8 coefs
             */
            q15x8_t vecCoeffs = *(q15x8_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vmlaldavaq(acc0, vecIn0, vecCoeffs);

            pSamplesTmp += 8;
            pCoeffsTmp += 8;
            /*
             * Decrement the taps block loop counter
             */
            i--;
        }

        *pOutput++ = (q15_t)MVE_ASRL_SAT16(acc0, 15);
    }
    break;
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;

    blkCnt = numTaps >> 3;
    while (blkCnt > 0U)
    {
        vst1q(pTempDest, vld1q(pTempSrc));
        pTempSrc += 8;
        pTempDest += 8;
        blkCnt--;
    }
    blkCnt = numTaps & 7;
    if (blkCnt > 0U)
    {
        mve_pred16_t p0 = vctp16q(blkCnt);
        vstrhq_p_s16(pTempDest, vld1q(pTempSrc), p0);
    }
}

#else

#if defined(ARM_MATH_NEON)

#if defined(SVE)

#define UPDATE_SUMS(COEF)                              \
    tempLO = vmull_n_s16(vget_low_s16(xa), (COEF));    \
    tempHI = vmull_n_s16(vget_high_s16(xa), (COEF));   \
    sumLL0 = vaddw_s32(sumLL0, vget_low_s32(tempLO));  \
    sumLH0 = vaddw_s32(sumLH0, vget_high_s32(tempLO)); \
    sumHL0 = vaddw_s32(sumHL0, vget_low_s32(tempHI));  \
    sumHH0 = vaddw_s32(sumHH0, vget_high_s32(tempHI)); \
    tempLO = vmull_n_s16(vget_low_s16(xb), (COEF));    \
    tempHI = vmull_n_s16(vget_high_s16(xb), (COEF));   \
    sumLL1 = vaddw_s32(sumLL1, vget_low_s32(tempLO));  \
    sumLH1 = vaddw_s32(sumLH1, vget_high_s32(tempLO)); \
    sumHL1 = vaddw_s32(sumHL1, vget_low_s32(tempHI));  \
    sumHH1 = vaddw_s32(sumHH1, vget_high_s32(tempHI));

#define SHIFT_VECTORS(d)       \
    xa = vextq_u16(x0, x1, d); \
    xb = vextq_u16(x1, x2, d);

#define UPDATE_SUMS_SVE(VOUT0, VOUT1, COEF, xa, xb) \
    VOUT0 = svdot(VOUT0, xa, COEF);                 \
    VOUT1 = svdot(VOUT1, xb, COEF);

#define SHIFT_VECTORS_SVE(d, x0, x1, x2)           \
    xa = int16x8_to_svint16(vextq_u16(x0, x1, d)); \
    xb = int16x8_to_svint16(vextq_u16(x1, x2, d));

#define UPDATE_VECTOR(xa, xb, x0, x1, x2)                    \
    UPDATE_SUMS_SVE(sve_sumLL0, sve_sumLL1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(1, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLH0, sve_sumLH1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(2, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHL0, sve_sumHL1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(3, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHH0, sve_sumHH1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(4, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLL0, sve_sumLL1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(5, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLH0, sve_sumLH1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(6, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHL0, sve_sumHL1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(7, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHH0, sve_sumHH1, vCoef1, xa, xb);

#define START_OF_FUNCTION_COMMON()                                                            \
    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */ \
    /* pStateCurnt points to the location where the new input data should be written */       \
    pStateCurnt = &(S->pState[(numTaps - 1U)]);                                               \
    int nTail = numTaps & 0x7;                                                                \
    int nTaps = numTaps >> 3;                                                                 \
    /* Initialize blkCnt with blockSize */                                                    \
    blkCnt = blockSize >> 4;                                                                  \
    while (blkCnt > 0U)                                                                       \
    {                                                                                         \
        /* Copy one sample at a time into state buffer */                                     \
        samples = vld1q_s16(pSrc);                                                            \
        vst1q_s16(pStateCurnt, samples);                                                      \
                                                                                              \
        pStateCurnt += 8;                                                                     \
        pSrc += 8;                                                                            \
                                                                                              \
        samples = vld1q_s16(pSrc);                                                            \
        vst1q_s16(pStateCurnt, samples);                                                      \
                                                                                              \
        pStateCurnt += 8;                                                                     \
        pSrc += 8;                                                                            \
                                                                                              \
        /* Set the accumulator to zero */                                                     \
                                                                                              \
        sve_sumLL0 = svdup_s64(0);                                                            \
        sve_sumLH0 = svdup_s64(0);                                                            \
        sve_sumHL0 = svdup_s64(0);                                                            \
        sve_sumHH0 = svdup_s64(0);                                                            \
                                                                                              \
        sve_sumLL1 = svdup_s64(0);                                                            \
        sve_sumLH1 = svdup_s64(0);                                                            \
        sve_sumHL1 = svdup_s64(0);                                                            \
        sve_sumHH1 = svdup_s64(0);                                                            \
                                                                                              \
        /* Initialize state pointer */                                                        \
        px = pState;

#define END_OF_FUNCTION_COMMON()                                                                                                         \
    sumLL0 = svint64_to_int64x2(sve_sumLL0);                                                                                             \
    sumLH0 = svint64_to_int64x2(sve_sumLH0);                                                                                             \
    sumHL0 = svint64_to_int64x2(sve_sumHL0);                                                                                             \
    sumHH0 = svint64_to_int64x2(sve_sumHH0);                                                                                             \
    sumLL1 = svint64_to_int64x2(sve_sumLL1);                                                                                             \
    sumLH1 = svint64_to_int64x2(sve_sumLH1);                                                                                             \
    sumHL1 = svint64_to_int64x2(sve_sumHL1);                                                                                             \
    sumHH1 = svint64_to_int64x2(sve_sumHH1);                                                                                             \
                                                                                                                                         \
    tempLLL = vqshrn_n_s64(sumLL0, 15);                                                                                                  \
    tempLHL = vqshrn_n_s64(sumLH0, 15);                                                                                                  \
    tempHLL = vqshrn_n_s64(sumHL0, 15);                                                                                                  \
    tempHHL = vqshrn_n_s64(sumHH0, 15);                                                                                                  \
    tempLLH = vqshrn_n_s64(sumLL1, 15);                                                                                                  \
    tempLHH = vqshrn_n_s64(sumLH1, 15);                                                                                                  \
    tempHLH = vqshrn_n_s64(sumHL1, 15);                                                                                                  \
    tempHHH = vqshrn_n_s64(sumHH1, 15);                                                                                                  \
                                                                                                                                         \
    tempL0 = vcombine_s32(tempLLL, tempLLH);                                                                                             \
    tempH0 = vcombine_s32(tempLHL, tempLHH);                                                                                             \
    tempL1 = vcombine_s32(tempHLL, tempHLH);                                                                                             \
    tempH1 = vcombine_s32(tempHHL, tempHHH);                                                                                             \
                                                                                                                                         \
    temp0L.val[0] = vqmovn_s32(tempL0);                                                                                                  \
    temp0L.val[1] = vqmovn_s32(tempH0);                                                                                                  \
    temp0L.val[2] = vqmovn_s32(tempL1);                                                                                                  \
    temp0L.val[3] = vqmovn_s32(tempH1);                                                                                                  \
                                                                                                                                         \
    vst4_s16(pDst, temp0L);                                                                                                              \
    pDst += 16;                                                                                                                          \
                                                                                                                                         \
    /* The result is stored in the destination buffer. */                                                                                \
                                                                                                                                         \
    /* Advance state pointer by 1 for the next sample */                                                                                 \
    pState = pState + 16;                                                                                                                \
    blkCnt--;                                                                                                                            \
    }                                                                                                                                    \
                                                                                                                                         \
    blkCnt = blockSize & 0xF;                                                                                                            \
    while (blkCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        /* Copy one sample at a time into state buffer */                                                                                \
        *pStateCurnt++ = *pSrc++;                                                                                                        \
                                                                                                                                         \
        /* Set the accumulator to zero */                                                                                                \
        acc = 0;                                                                                                                         \
                                                                                                                                         \
        /* Initialize state pointer */                                                                                                   \
        px = pState;                                                                                                                     \
                                                                                                                                         \
        /* Initialize Coefficient pointer */                                                                                             \
        pb = pCoeffs;                                                                                                                    \
                                                                                                                                         \
        i = numTaps;                                                                                                                     \
                                                                                                                                         \
        /* Perform the multiply-accumulates */                                                                                           \
        do                                                                                                                               \
        {                                                                                                                                \
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
            acc += (q31_t) * px++ * *pb++;                                                                                               \
            i--;                                                                                                                         \
        } while (i > 0U);                                                                                                                \
                                                                                                                                         \
        /* The result is stored in the destination buffer. */                                                                            \
        *pDst++ = (q15_t)__SSAT((acc >> 15U), 16);                                                                                       \
                                                                                                                                         \
        /* Advance state pointer by 1 for the next sample */                                                                             \
        pState = pState + 1;                                                                                                             \
                                                                                                                                         \
        blkCnt--;                                                                                                                        \
    }                                                                                                                                    \
                                                                                                                                         \
    /* Processing is complete.                                                                                                           \
    ** Now copy the last numTaps - 1 samples to the starting of the state buffer.                                                        \
    ** This prepares the state buffer for the next function call. */                                                                     \
                                                                                                                                         \
    /* Points to the start of the state buffer */                                                                                        \
    pStateCurnt = S->pState;                                                                                                             \
                                                                                                                                         \
    /* Copy numTaps number of values */                                                                                                  \
    tapCnt = numTaps - 1U;                                                                                                               \
                                                                                                                                         \
    /* Copy data */                                                                                                                      \
    while (tapCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        *pStateCurnt++ = *pState++;                                                                                                      \
                                                                                                                                         \
        /* Decrement the loop counter */                                                                                                 \
        tapCnt--;                                                                                                                        \
    }

#define UPDATE_HALF_SUM_TAIL_1(x0, x1, x2)                   \
    UPDATE_SUMS_SVE(sve_sumLL0, sve_sumLL1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(1, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLH0, sve_sumLH1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(2, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHL0, sve_sumHL1, vCoef0, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(3, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHH0, sve_sumHH1, vCoef0, xa, xb);

#define UPDATE_HALF_SUM_TAIL_2(x0, x1, x2)                   \
    SHIFT_VECTORS_SVE(4, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLL0, sve_sumLL1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(5, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumLH0, sve_sumLH1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(6, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHL0, sve_sumHL1, vCoef1, xa, xb); \
                                                             \
    SHIFT_VECTORS_SVE(7, x0, x1, x2);                        \
    UPDATE_SUMS_SVE(sve_sumHH0, sve_sumHH1, vCoef1, xa, xb);

#define TAIL_FIR_Q15_S(x0, x1, x2)                                   \
    switch (nTail)                                                   \
    {                                                                \
    case 1:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        svbool_t sv_selector = svdupq_n_b16(1, 0, 0, 0, 1, 0, 0, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef0 = svsel_s16(sv_selector, vCoef0, sve_zero);           \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        break;                                                       \
    }                                                                \
    case 2:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        svbool_t sv_selector = svdupq_n_b16(1, 1, 0, 0, 1, 1, 0, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef0 = svsel_s16(sv_selector, vCoef0, sve_zero);           \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
        break;                                                       \
    }                                                                \
    case 3:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        svbool_t sv_selector = svdupq_n_b16(1, 1, 1, 0, 1, 1, 1, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef0 = svsel_s16(sv_selector, vCoef0, sve_zero);           \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
        break;                                                       \
    }                                                                \
    case 4:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
        break;                                                       \
    }                                                                \
    case 5:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
                                                                     \
        svbool_t sv_selector = svdupq_n_b16(1, 0, 0, 0, 1, 0, 0, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef1 = svsel_s16(sv_selector, vCoef1, sve_zero);           \
                                                                     \
        UPDATE_HALF_SUM_TAIL_2(x0, x1, x2);                          \
        break;                                                       \
    }                                                                \
    case 6:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
                                                                     \
        svbool_t sv_selector = svdupq_n_b16(1, 1, 0, 0, 1, 1, 0, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef1 = svsel_s16(sv_selector, vCoef1, sve_zero);           \
                                                                     \
        UPDATE_HALF_SUM_TAIL_2(x0, x1, x2);                          \
        break;                                                       \
    }                                                                \
    case 7:                                                          \
    {                                                                \
        x2 = vld1q_s16(px + 16);                                     \
        xa = int16x8_to_svint16(x0);                                 \
        xb = int16x8_to_svint16(x1);                                 \
        UPDATE_HALF_SUM_TAIL_1(x0, x1, x2);                          \
        ;                                                            \
                                                                     \
        svbool_t sv_selector = svdupq_n_b16(1, 1, 1, 0, 1, 1, 1, 0); \
        svint32_t sve_zero = svdup_n_s32(0);                         \
        vCoef1 = svsel_s16(sv_selector, vCoef1, sve_zero);           \
                                                                     \
        UPDATE_HALF_SUM_TAIL_2(x0, x1, x2);                          \
        break;                                                       \
    }                                                                \
    default:                                                         \
        break;                                                       \
    }

#define INNER_LOOP_GENERIC()                                                                                                         \
    /* Initialize Coefficient pointer */                                                                                             \
    pb = pCoeffs;                                                                                                                    \
    i = numTaps >> 3;                                                                                                                \
                                                                                                                                     \
    /* Perform the multiply-accumulates */                                                                                           \
    x0 = vld1q_s16(px);                                                                                                              \
    x1 = vld1q_s16(px + 8);                                                                                                          \
                                                                                                                                     \
    b[0] = vld1q_s16(pb);                                                                                                            \
    svint16_t vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));                                     \
    svint16_t vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0])));                                   \
                                                                                                                                     \
    /* svptest_any(svptrue_b16(), pg) */                                                                                             \
    while (i > 0)                                                                                                                    \
    {                                                                                                                                \
        /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
                                                                                                                                     \
        x2 = vld1q_s16(px + 16);                                                                                                     \
        xa = int16x8_to_svint16(x0);                                                                                                 \
        xb = int16x8_to_svint16(x1);                                                                                                 \
                                                                                                                                     \
        UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                                                           \
                                                                                                                                     \
        pb += 8;                                                                                                                     \
        b[0] = vld1q_s16(pb);                                                                                                        \
        vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));                                           \
        vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0])));                                         \
        x0 = x1;                                                                                                                     \
        x1 = x2;                                                                                                                     \
        px += 8;                                                                                                                     \
        i--;                                                                                                                         \
    }                                                                                                                                \
    TAIL_FIR_Q15_S(x0, x1, x2);

#define INNER_LOOP_COEF_8_15()                                                                     \
    /* Initialize Coefficient pointer */                                                           \
    /* Perform the multiply-accumulates */                                                         \
    x0 = vld1q_s16(px);                                                                            \
    x1 = vld1q_s16(px + 8);                                                                        \
                                                                                                   \
    svint16_t vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));   \
    svint16_t vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0]))); \
                                                                                                   \
    x2 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x0);                                                                   \
    xb = int16x8_to_svint16(x1);                                                                   \
                                                                                                   \
    UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                             \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[1]), vget_low_s16(b[1])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[1]), vget_high_s16(b[1])));           \
    px += 8;                                                                                       \
    TAIL_FIR_Q15_S(x1, x2, x0);

#define INNER_LOOP_COEF_16_23()                                                                    \
    /* Initialize Coefficient pointer */                                                           \
    /* Perform the multiply-accumulates */                                                         \
    x0 = vld1q_s16(px);                                                                            \
    x1 = vld1q_s16(px + 8);                                                                        \
                                                                                                   \
    svint16_t vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));   \
    svint16_t vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0]))); \
                                                                                                   \
    x2 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x0);                                                                   \
    xb = int16x8_to_svint16(x1);                                                                   \
                                                                                                   \
    UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                             \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[1]), vget_low_s16(b[1])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[1]), vget_high_s16(b[1])));           \
    px += 8;                                                                                       \
    x0 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x1);                                                                   \
    xb = int16x8_to_svint16(x2);                                                                   \
    UPDATE_VECTOR(xa, xb, x1, x2, x0);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[2]), vget_low_s16(b[2])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[2]), vget_high_s16(b[2])));           \
    px += 8;                                                                                       \
    TAIL_FIR_Q15_S(x2, x0, x1);

#define INNER_LOOP_COEF_24_31()                                                                    \
    /* Initialize Coefficient pointer */                                                           \
    /* Perform the multiply-accumulates */                                                         \
    x0 = vld1q_s16(px);                                                                            \
    x1 = vld1q_s16(px + 8);                                                                        \
                                                                                                   \
    svint16_t vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));   \
    svint16_t vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0]))); \
                                                                                                   \
    x2 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x0);                                                                   \
    xb = int16x8_to_svint16(x1);                                                                   \
                                                                                                   \
    UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                             \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[1]), vget_low_s16(b[1])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[1]), vget_high_s16(b[1])));           \
    px += 8;                                                                                       \
    x0 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x1);                                                                   \
    xb = int16x8_to_svint16(x2);                                                                   \
    UPDATE_VECTOR(xa, xb, x1, x2, x0);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[2]), vget_low_s16(b[2])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[2]), vget_high_s16(b[2])));           \
    px += 8;                                                                                       \
    x1 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x2);                                                                   \
    xb = int16x8_to_svint16(x0);                                                                   \
    UPDATE_VECTOR(xa, xb, x2, x0, x1);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[3]), vget_low_s16(b[3])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[3]), vget_high_s16(b[3])));           \
    px += 8;                                                                                       \
    TAIL_FIR_Q15_S(x0, x1, x2);

#define INNER_LOOP_COEF_32_39()                                                                    \
    /* Initialize Coefficient pointer */                                                           \
    /* Perform the multiply-accumulates */                                                         \
    x0 = vld1q_s16(px);                                                                            \
    x1 = vld1q_s16(px + 8);                                                                        \
                                                                                                   \
    svint16_t vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[0]), vget_low_s16(b[0])));   \
    svint16_t vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[0]), vget_high_s16(b[0]))); \
                                                                                                   \
    x2 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x0);                                                                   \
    xb = int16x8_to_svint16(x1);                                                                   \
                                                                                                   \
    UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                             \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[1]), vget_low_s16(b[1])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[1]), vget_high_s16(b[1])));           \
    px += 8;                                                                                       \
    x0 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x1);                                                                   \
    xb = int16x8_to_svint16(x2);                                                                   \
    UPDATE_VECTOR(xa, xb, x1, x2, x0);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[2]), vget_low_s16(b[2])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[2]), vget_high_s16(b[2])));           \
    px += 8;                                                                                       \
    x1 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x2);                                                                   \
    xb = int16x8_to_svint16(x0);                                                                   \
    UPDATE_VECTOR(xa, xb, x2, x0, x1);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[3]), vget_low_s16(b[3])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[3]), vget_high_s16(b[3])));           \
    px += 8;                                                                                       \
    x2 = vld1q_s16(px + 16);                                                                       \
    xa = int16x8_to_svint16(x0);                                                                   \
    xb = int16x8_to_svint16(x1);                                                                   \
    UPDATE_VECTOR(xa, xb, x0, x1, x2);                                                             \
                                                                                                   \
    vCoef0 = int16x8_to_svint16(vcombine_s16(vget_low_s16(b[4]), vget_low_s16(b[4])));             \
    vCoef1 = int16x8_to_svint16(vcombine_s16(vget_high_s16(b[4]), vget_high_s16(b[4])));           \
    px += 8;                                                                                       \
    TAIL_FIR_Q15_S(x1, x2, x0);

void arm_fir_q15(const arm_fir_instance_q15 *S,
                 const q15_t *pSrc,
                 q15_t *pDst,
                 uint32_t blockSize)
{
    q15_t *pState = S->pState;         /* State pointer */
    const q15_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q15_t *pStateCurnt;                /* Points to the current sample of the state */
    q15_t *px;                         /* Temporary pointers for state buffer */
    const q15_t *pb;                   /* Temporary pointers for coefficient buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t i, tapCnt, blkCnt;        /* Loop counters */

    int32x2_t tempLLL, tempLHL, tempHLL, tempHHL;
    int32x2_t tempLLH, tempLHH, tempHLH, tempHHH;
    int32x4_t tempL0, tempH1, tempL1, tempH0;

    int16x4x4_t temp0L;
    int16x8_t temp;
    int32x4_t tempLO, tempHI;
    int64x2_t sumLL0, sumLH0, sumHL0, sumHH0;
    int64x2_t sumLL1, sumLH1, sumHL1, sumHH1;

    svint64_t sve_sumLL0, sve_sumLH0, sve_sumHL0, sve_sumHH0, sve_sumLL1, sve_sumLH1, sve_sumHL1, sve_sumHH1;
    svint16_t xa, xb;

    int16x8_t samples, x0, x1, x2;
    q63_t acc;
    int32_t nbTaps = (numTaps) >> 3;
    int nTail = S->numTaps & 0x7;
    switch (nbTaps)
    {
    case 1:
    {
        pb = pCoeffs;
        int16x8_t b[2];
        b[0] = vld1q_s16(pb);
        if (nTail > 0)
            b[1] = vld1q_s16(pb + 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_8_15();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 2:
    {
        pb = pCoeffs;
        int16x8_t b[3];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[2] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_16_23();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 3:
    {
        pb = pCoeffs;
        int16x8_t b[4];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[3] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_24_31();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 4:
    {
        pb = pCoeffs;
        int16x8_t b[5];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[4] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_32_39();
        END_OF_FUNCTION_COMMON();
        return;
    }
    default:
    {
        int32x4_t b[1];
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_GENERIC();

        END_OF_FUNCTION_COMMON();
        return;
    }
    }
}

#else

#define UPDATE_SUMS(COEF)                              \
    tempLO = vmull_n_s16(vget_low_s16(xa), (COEF));    \
    tempHI = vmull_n_s16(vget_high_s16(xa), (COEF));   \
    sumLL0 = vaddw_s32(sumLL0, vget_low_s32(tempLO));  \
    sumLH0 = vaddw_s32(sumLH0, vget_high_s32(tempLO)); \
    sumHL0 = vaddw_s32(sumHL0, vget_low_s32(tempHI));  \
    sumHH0 = vaddw_s32(sumHH0, vget_high_s32(tempHI)); \
    tempLO = vmull_n_s16(vget_low_s16(xb), (COEF));    \
    tempHI = vmull_n_s16(vget_high_s16(xb), (COEF));   \
    sumLL1 = vaddw_s32(sumLL1, vget_low_s32(tempLO));  \
    sumLH1 = vaddw_s32(sumLH1, vget_high_s32(tempLO)); \
    sumHL1 = vaddw_s32(sumHL1, vget_low_s32(tempHI));  \
    sumHH1 = vaddw_s32(sumHH1, vget_high_s32(tempHI));

#define SHIFT_VECTORS(x0, x1, x2, d) \
    xa = vextq_u16(x0, x1, d);       \
    xb = vextq_u16(x1, x2, d);

#define UPDATE_VECTOR(xa, xb, x0, x1, x2, b) \
    UPDATE_SUMS(b[0]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 1);            \
    UPDATE_SUMS(b[1]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 2);            \
    UPDATE_SUMS(b[2]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 3);            \
    UPDATE_SUMS(b[3]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 4);            \
    UPDATE_SUMS(b[4]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 5);            \
    UPDATE_SUMS(b[5]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 6);            \
    UPDATE_SUMS(b[6]);                       \
                                             \
    SHIFT_VECTORS(x0, x1, x2, 7);            \
    UPDATE_SUMS(b[7]);

#define START_OF_FUNCTION_COMMON()                                                            \
    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */ \
    /* pStateCurnt points to the location where the new input data should be written */       \
    pStateCurnt = &(S->pState[(numTaps - 1U)]);                                               \
                                                                                              \
    /* Initialize blkCnt with blockSize */                                                    \
    blkCnt = blockSize >> 4;                                                                  \
    while (blkCnt > 0U)                                                                       \
    {                                                                                         \
        /* Copy one sample at a time into state buffer */                                     \
        samples = vld1q_s16(pSrc);                                                            \
        vst1q_s16(pStateCurnt, samples);                                                      \
                                                                                              \
        pStateCurnt += 8;                                                                     \
        pSrc += 8;                                                                            \
                                                                                              \
        samples = vld1q_s16(pSrc);                                                            \
        vst1q_s16(pStateCurnt, samples);                                                      \
                                                                                              \
        pStateCurnt += 8;                                                                     \
        pSrc += 8;                                                                            \
                                                                                              \
        /* Set the accumulator to zero */                                                     \
                                                                                              \
        sumLL0 = vdupq_n_s64(0);                                                              \
        sumLH0 = vdupq_n_s64(0);                                                              \
        sumHL0 = vdupq_n_s64(0);                                                              \
        sumHH0 = vdupq_n_s64(0);                                                              \
                                                                                              \
        sumLL1 = vdupq_n_s64(0);                                                              \
        sumLH1 = vdupq_n_s64(0);                                                              \
        sumHL1 = vdupq_n_s64(0);                                                              \
        sumHH1 = vdupq_n_s64(0);                                                              \
                                                                                              \
        /* Initialize state pointer */                                                        \
        px = pState;                                                                          \
                                                                                              \
        /* Initialize Coefficient pointer */                                                  \
        pb = pCoeffs;                                                                         \
                                                                                              \
        i = numTaps >> 3;

#define END_OF_FUNCTION_COMMON()                                                                                                         \
    tempLL = vqshrn_n_s64(sumLL0, 15);                                                                                                   \
    tempLH = vqshrn_n_s64(sumLH0, 15);                                                                                                   \
    tempHL = vqshrn_n_s64(sumHL0, 15);                                                                                                   \
    tempHH = vqshrn_n_s64(sumHH0, 15);                                                                                                   \
                                                                                                                                         \
    tempL = vcombine_s32(tempLL, tempLH);                                                                                                \
    tempH = vcombine_s32(tempHL, tempHH);                                                                                                \
                                                                                                                                         \
    temp0 = vqmovn_s32(tempL);                                                                                                           \
    temp1 = vqmovn_s32(tempH);                                                                                                           \
    temp = vcombine_s16(temp0, temp1);                                                                                                   \
                                                                                                                                         \
    /* The result is stored in the destination buffer. */                                                                                \
    vst1q_s16(pDst, temp);                                                                                                               \
                                                                                                                                         \
    pDst += 8;                                                                                                                           \
                                                                                                                                         \
    tempLL = vqshrn_n_s64(sumLL1, 15);                                                                                                   \
    tempLH = vqshrn_n_s64(sumLH1, 15);                                                                                                   \
    tempHL = vqshrn_n_s64(sumHL1, 15);                                                                                                   \
    tempHH = vqshrn_n_s64(sumHH1, 15);                                                                                                   \
                                                                                                                                         \
    tempL = vcombine_s32(tempLL, tempLH);                                                                                                \
    tempH = vcombine_s32(tempHL, tempHH);                                                                                                \
                                                                                                                                         \
    temp0 = vqmovn_s32(tempL);                                                                                                           \
    temp1 = vqmovn_s32(tempH);                                                                                                           \
    temp = vcombine_s16(temp0, temp1);                                                                                                   \
                                                                                                                                         \
    /* The result is stored in the destination buffer. */                                                                                \
    vst1q_s16(pDst, temp);                                                                                                               \
                                                                                                                                         \
    pDst += 8;                                                                                                                           \
                                                                                                                                         \
    /* Advance state pointer by 1 for the next sample */                                                                                 \
    pState = pState + 16;                                                                                                                \
                                                                                                                                         \
    blkCnt--;                                                                                                                            \
    }                                                                                                                                    \
                                                                                                                                         \
    blkCnt = blockSize & 0xF;                                                                                                            \
    while (blkCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        /* Copy one sample at a time into state buffer */                                                                                \
        *pStateCurnt++ = *pSrc++;                                                                                                        \
                                                                                                                                         \
        /* Set the accumulator to zero */                                                                                                \
        acc = 0;                                                                                                                         \
                                                                                                                                         \
        /* Initialize state pointer */                                                                                                   \
        px = pState;                                                                                                                     \
                                                                                                                                         \
        /* Initialize Coefficient pointer */                                                                                             \
        pb = pCoeffs;                                                                                                                    \
                                                                                                                                         \
        i = numTaps;                                                                                                                     \
                                                                                                                                         \
        /* Perform the multiply-accumulates */                                                                                           \
        do                                                                                                                               \
        {                                                                                                                                \
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
            acc += (q31_t) * px++ * *pb++;                                                                                               \
            i--;                                                                                                                         \
        } while (i > 0U);                                                                                                                \
                                                                                                                                         \
        /* The result is stored in the destination buffer. */                                                                            \
        *pDst++ = (q15_t)__SSAT((acc >> 15U), 16);                                                                                       \
                                                                                                                                         \
        /* Advance state pointer by 1 for the next sample */                                                                             \
        pState = pState + 1;                                                                                                             \
                                                                                                                                         \
        blkCnt--;                                                                                                                        \
    }                                                                                                                                    \
                                                                                                                                         \
    /* Processing is complete.                                                                                                           \
    ** Now copy the last numTaps - 1 samples to the starting of the state buffer.                                                        \
    ** This prepares the state buffer for the next function call. */                                                                     \
                                                                                                                                         \
    /* Points to the start of the state buffer */                                                                                        \
    pStateCurnt = S->pState;                                                                                                             \
                                                                                                                                         \
    /* Copy numTaps number of values */                                                                                                  \
    tapCnt = numTaps - 1U;                                                                                                               \
                                                                                                                                         \
    /* Copy data */                                                                                                                      \
    while (tapCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        *pStateCurnt++ = *pState++;                                                                                                      \
                                                                                                                                         \
        /* Decrement the loop counter */                                                                                                 \
        tapCnt--;                                                                                                                        \
    }

#define TAIL_FIR_Q15_S(x0, x1, x2)             \
    switch (nTail)                             \
    {                                          \
    case 1:                                    \
    {                                          \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
                                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
        break;                                 \
    }                                          \
    case 2:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
        break;                                 \
    }                                          \
    case 3:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 2);          \
        UPDATE_SUMS(coefs1_32_1[2]);           \
        break;                                 \
    }                                          \
    case 4:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 2);          \
        UPDATE_SUMS(coefs1_32_1[2]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 3);          \
        UPDATE_SUMS(coefs1_32_1[3]);           \
        break;                                 \
    }                                          \
    case 5:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 2);          \
        UPDATE_SUMS(coefs1_32_1[2]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 3);          \
        UPDATE_SUMS(coefs1_32_1[3]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 4);          \
        UPDATE_SUMS(coefs1_32_1[4]);           \
        break;                                 \
    }                                          \
    case 6:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 2);          \
        UPDATE_SUMS(coefs1_32_1[2]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 3);          \
        UPDATE_SUMS(coefs1_32_1[3]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 4);          \
        UPDATE_SUMS(coefs1_32_1[4]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 5);          \
        UPDATE_SUMS(coefs1_32_1[5]);           \
        break;                                 \
    }                                          \
    case 7:                                    \
    {                                          \
        int16x8_t coefs1_32_1 = vld1q_s16(pb); \
        x2 = vld1q_s16(px + 16);               \
        xa = x0;                               \
        xb = x1;                               \
        UPDATE_SUMS(coefs1_32_1[0]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 1);          \
        UPDATE_SUMS(coefs1_32_1[1]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 2);          \
        UPDATE_SUMS(coefs1_32_1[2]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 3);          \
        UPDATE_SUMS(coefs1_32_1[3]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 4);          \
        UPDATE_SUMS(coefs1_32_1[4]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 5);          \
        UPDATE_SUMS(coefs1_32_1[5]);           \
                                               \
        SHIFT_VECTORS(x0, x1, x2, 6);          \
        UPDATE_SUMS(coefs1_32_1[6]);           \
                                               \
        break;                                 \
    }                                          \
    default:                                   \
        break;                                 \
    }

#define INNER_LOOP_GENERIC()                                                                                                         \
    /* Perform the multiply-accumulates */                                                                                           \
    x0 = vld1q_s16(px);                                                                                                              \
    x1 = vld1q_s16(px + 8);                                                                                                          \
    b = vld1q_s16(pb);                                                                                                               \
    while (i > 0)                                                                                                                    \
    {                                                                                                                                \
        /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
        x2 = vld1q_s16(px + 16);                                                                                                     \
        xa = x0;                                                                                                                     \
        xb = x1;                                                                                                                     \
        UPDATE_VECTOR(xa, xb, x0, x1, x2, b)                                                                                         \
        UPDATE_SUMS(b[7]);                                                                                                           \
                                                                                                                                     \
        pb += 8;                                                                                                                     \
        b = vld1q_s16(pb);                                                                                                           \
        x0 = x1;                                                                                                                     \
        x1 = x2;                                                                                                                     \
        px += 8;                                                                                                                     \
        i--;                                                                                                                         \
    }                                                                                                                                \
    TAIL_FIR_Q15_S(x0, x1, x2);

#define INNER_LOOP_COEF_8_15()                                                                                                   \
    /* Initialize Coefficient pointer */                                                                                         \
    /* Perform the multiply-accumulates */                                                                                       \
    x0 = vld1q_s16(px);                                                                                                          \
    x1 = vld1q_s16(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x2 = vld1q_s16(px + 16);                                                                                                     \
    xa = x0;                                                                                                                     \
    xb = x1;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x0, x1, x2, b[0]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    TAIL_FIR_Q15_S(x1, x2, x0);

#define INNER_LOOP_COEF_16_23()                                                                                                  \
    /* Initialize Coefficient pointer */                                                                                         \
    /* Perform the multiply-accumulates */                                                                                       \
    x0 = vld1q_s16(px);                                                                                                          \
    x1 = vld1q_s16(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x2 = vld1q_s16(px + 16);                                                                                                     \
    xa = x0;                                                                                                                     \
    xb = x1;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x0, x1, x2, b[0]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x0 = vld1q_s16(px + 16);                                                                                                     \
    xa = x1;                                                                                                                     \
    xb = x2;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x1, x2, x0, b[1]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    TAIL_FIR_Q15_S(x2, x0, x1);

#define INNER_LOOP_COEF_24_31()                                                                                                  \
    /* Initialize Coefficient pointer */                                                                                         \
    /* Perform the multiply-accumulates */                                                                                       \
    x0 = vld1q_s16(px);                                                                                                          \
    x1 = vld1q_s16(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x2 = vld1q_s16(px + 16);                                                                                                     \
    xa = x0;                                                                                                                     \
    xb = x1;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x0, x1, x2, b[0]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x0 = vld1q_s16(px + 16);                                                                                                     \
    xa = x1;                                                                                                                     \
    xb = x2;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x1, x2, x0, b[1]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x1 = vld1q_s16(px + 16);                                                                                                     \
    xa = x2;                                                                                                                     \
    xb = x0;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x2, x0, x1, b[2]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    TAIL_FIR_Q15_S(x0, x1, x2);

#define INNER_LOOP_COEF_32_39()                                                                                                  \
    /* Initialize Coefficient pointer */                                                                                         \
    /* Perform the multiply-accumulates */                                                                                       \
    x0 = vld1q_s16(px);                                                                                                          \
    x1 = vld1q_s16(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x2 = vld1q_s16(px + 16);                                                                                                     \
    xa = x0;                                                                                                                     \
    xb = x1;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x0, x1, x2, b[0]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x0 = vld1q_s16(px + 16);                                                                                                     \
    xa = x1;                                                                                                                     \
    xb = x2;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x1, x2, x0, b[1]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x1 = vld1q_s16(px + 16);                                                                                                     \
    xa = x2;                                                                                                                     \
    xb = x0;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x2, x0, x1, b[2]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    x2 = vld1q_s16(px + 16);                                                                                                     \
    xa = x0;                                                                                                                     \
    xb = x1;                                                                                                                     \
    UPDATE_VECTOR(xa, xb, x0, x1, x2, b[3]);                                                                                     \
                                                                                                                                 \
    pb += 8;                                                                                                                     \
    px += 8;                                                                                                                     \
    TAIL_FIR_Q15_S(x1, x2, x0);

void arm_fir_q15(const arm_fir_instance_q15 *S, const q15_t *pSrc, q15_t *pDst, uint32_t blockSize)
{
    q15_t *pState = S->pState;         /* State pointer */
    const q15_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q15_t *pStateCurnt;                /* Points to the current sample of the state */
    q15_t *px;                         /* Temporary pointers for state buffer */
    const q15_t *pb;                   /* Temporary pointers for coefficient buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t i, tapCnt, blkCnt;        /* Loop counters */

    int32x2_t tempLL, tempLH, tempHL, tempHH;
    int32x4_t tempL, tempH;
    int16x4_t temp0, temp1;
    int16x8_t temp;
    int32x4_t tempLO, tempHI;
    int64x2_t sumLL0, sumLH0, sumHL0, sumHH0;
    int64x2_t sumLL1, sumLH1, sumHL1, sumHH1;

    int16x8_t samples, x0, x1, x2, xa, xb;
    q63_t acc;

    int32_t nbTaps = (numTaps) >> 3;
    int nTail = S->numTaps & 0x7;
    switch (nbTaps)
    {
    case 1:
    {
        pb = pCoeffs;
        int16x8_t b[2];
        b[0] = vld1q_s16(pb);
        if (nTail > 0)
            b[1] = vld1q_s16(pb + 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_8_15();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 2:
    {
        pb = pCoeffs;
        int16x8_t b[3];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[2] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_16_23();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 3:
    {
        pb = pCoeffs;
        int16x8_t b[4];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[3] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_24_31();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 4:
    {
        pb = pCoeffs;
        int16x8_t b[5];

        for (int i = 0; i < nbTaps; i++)
            b[i] = vld1q_s16(pb + i * 8);
        if (nTail > 0)
            b[4] = vld1q_s16(pb + nbTaps * 8);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_32_39();
        END_OF_FUNCTION_COMMON();
        return;
    }
    default:
    {
        int16x8_t b;
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_GENERIC();
        END_OF_FUNCTION_COMMON();
        return;
    }
    }
}

#endif

#else
ARM_DSP_ATTRIBUTE void arm_fir_q15(
    const arm_fir_instance_q15 *S,
    const q15_t *pSrc,
    q15_t *pDst,
    uint32_t blockSize)
{
    q15_t *pState = S->pState;         /* State pointer */
    const q15_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q15_t *pStateCurnt;                /* Points to the current sample of the state */
    q15_t *px;                         /* Temporary pointer for state buffer */
    const q15_t *pb;                   /* Temporary pointer for coefficient buffer */
    q63_t acc0;                        /* Accumulators */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t tapCnt, blkCnt;           /* Loop counters */

#if defined(ARM_MATH_LOOPUNROLL)
    q63_t acc1, acc2, acc3; /* Accumulators */
    q31_t x0, x1, x2, c0;   /* Temporary variables to hold state and coefficient values */
#endif

    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */
    /* pStateCurnt points to the location where the new input data should be written */
    pStateCurnt = &(S->pState[(numTaps - 1U)]);

#if defined(ARM_MATH_LOOPUNROLL)

    /* Loop unrolling: Compute 4 output values simultaneously.
     * The variables acc0 ... acc3 hold output values that are being computed:
     *
     *    acc0 =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0]
     *    acc1 =  b[numTaps-1] * x[n-numTaps]   + b[numTaps-2] * x[n-numTaps-1] + b[numTaps-3] * x[n-numTaps-2] +...+ b[0] * x[1]
     *    acc2 =  b[numTaps-1] * x[n-numTaps+1] + b[numTaps-2] * x[n-numTaps]   + b[numTaps-3] * x[n-numTaps-1] +...+ b[0] * x[2]
     *    acc3 =  b[numTaps-1] * x[n-numTaps+2] + b[numTaps-2] * x[n-numTaps+1] + b[numTaps-3] * x[n-numTaps]   +...+ b[0] * x[3]
     */
    blkCnt = blockSize >> 2U;

    while (blkCnt > 0U)
    {
        /* Copy 4 new input samples into the state buffer. */
        *pStateCurnt++ = *pSrc++;
        *pStateCurnt++ = *pSrc++;
        *pStateCurnt++ = *pSrc++;
        *pStateCurnt++ = *pSrc++;

        /* Set all accumulators to zero */
        acc0 = 0;
        acc1 = 0;
        acc2 = 0;
        acc3 = 0;

        /* Typecast q15_t pointer to q31_t pointer for state reading in q31_t */
        px = pState;

        /* Typecast q15_t pointer to q31_t pointer for coefficient reading in q31_t */
        pb = pCoeffs;

        /* Read the first two samples from the state buffer:  x[n-N], x[n-N-1] */
        x0 = read_q15x2_ia(&px);

        /* Read the third and forth samples from the state buffer: x[n-N-2], x[n-N-3] */
        x2 = read_q15x2_ia(&px);

        /* Loop over the number of taps.  Unroll by a factor of 4.
           Repeat until we've computed numTaps-(numTaps%4) coefficients. */
        tapCnt = numTaps >> 2U;

        while (tapCnt > 0U)
        {
            /* Read the first two coefficients using SIMD:  b[N] and b[N-1] coefficients */
            c0 = read_q15x2_ia(&pb);

            /* acc0 +=  b[N] * x[n-N] + b[N-1] * x[n-N-1] */
            acc0 = __SMLALD(x0, c0, acc0);

            /* acc2 +=  b[N] * x[n-N-2] + b[N-1] * x[n-N-3] */
            acc2 = __SMLALD(x2, c0, acc2);

            /* pack  x[n-N-1] and x[n-N-2] */
#ifndef ARM_MATH_BIG_ENDIAN
            x1 = __PKHBT(x2, x0, 0);
#else
            x1 = __PKHBT(x0, x2, 0);
#endif

            /* Read state x[n-N-4], x[n-N-5] */
            x0 = read_q15x2_ia(&px);

            /* acc1 +=  b[N] * x[n-N-1] + b[N-1] * x[n-N-2] */
            acc1 = __SMLALDX(x1, c0, acc1);

            /* pack  x[n-N-3] and x[n-N-4] */
#ifndef ARM_MATH_BIG_ENDIAN
            x1 = __PKHBT(x0, x2, 0);
#else
            x1 = __PKHBT(x2, x0, 0);
#endif

            /* acc3 +=  b[N] * x[n-N-3] + b[N-1] * x[n-N-4] */
            acc3 = __SMLALDX(x1, c0, acc3);

            /* Read coefficients b[N-2], b[N-3] */
            c0 = read_q15x2_ia(&pb);

            /* acc0 +=  b[N-2] * x[n-N-2] + b[N-3] * x[n-N-3] */
            acc0 = __SMLALD(x2, c0, acc0);

            /* Read state x[n-N-6], x[n-N-7] with offset */
            x2 = read_q15x2_ia(&px);

            /* acc2 +=  b[N-2] * x[n-N-4] + b[N-3] * x[n-N-5] */
            acc2 = __SMLALD(x0, c0, acc2);

            /* acc1 +=  b[N-2] * x[n-N-3] + b[N-3] * x[n-N-4] */
            acc1 = __SMLALDX(x1, c0, acc1);

            /* pack  x[n-N-5] and x[n-N-6] */
#ifndef ARM_MATH_BIG_ENDIAN
            x1 = __PKHBT(x2, x0, 0);
#else
            x1 = __PKHBT(x0, x2, 0);
#endif

            /* acc3 +=  b[N-2] * x[n-N-5] + b[N-3] * x[n-N-6] */
            acc3 = __SMLALDX(x1, c0, acc3);

            /* Decrement tap count */
            tapCnt--;
        }

        /* If the filter length is not a multiple of 4, compute the remaining filter taps.
           This is always be 2 taps since the filter length is even. */
        if ((numTaps & 0x3U) != 0U)
        {
            /* Read last two coefficients */
            c0 = read_q15x2_ia(&pb);

            /* Perform the multiply-accumulates */
            acc0 = __SMLALD(x0, c0, acc0);
            acc2 = __SMLALD(x2, c0, acc2);

            /* pack state variables */
#ifndef ARM_MATH_BIG_ENDIAN
            x1 = __PKHBT(x2, x0, 0);
#else
            x1 = __PKHBT(x0, x2, 0);
#endif

            /* Read last state variables */
            x0 = read_q15x2(px);

            /* Perform the multiply-accumulates */
            acc1 = __SMLALDX(x1, c0, acc1);

            /* pack state variables */
#ifndef ARM_MATH_BIG_ENDIAN
            x1 = __PKHBT(x0, x2, 0);
#else
            x1 = __PKHBT(x2, x0, 0);
#endif

            /* Perform the multiply-accumulates */
            acc3 = __SMLALDX(x1, c0, acc3);
        }

        /* The results in the 4 accumulators are in 2.30 format. Convert to 1.15 with saturation.
           Then store the 4 outputs in the destination buffer. */
#ifndef ARM_MATH_BIG_ENDIAN
        write_q15x2_ia(&pDst, __PKHBT(__SSAT((acc0 >> 15), 16), __SSAT((acc1 >> 15), 16), 16));
        write_q15x2_ia(&pDst, __PKHBT(__SSAT((acc2 >> 15), 16), __SSAT((acc3 >> 15), 16), 16));
#else
        write_q15x2_ia(&pDst, __PKHBT(__SSAT((acc1 >> 15), 16), __SSAT((acc0 >> 15), 16), 16));
        write_q15x2_ia(&pDst, __PKHBT(__SSAT((acc3 >> 15), 16), __SSAT((acc2 >> 15), 16), 16));
#endif /* #ifndef ARM_MATH_BIG_ENDIAN */

        /* Advance the state pointer by 4 to process the next group of 4 samples */
        pState = pState + 4U;

        /* Decrement loop counter */
        blkCnt--;
    }

    /* Loop unrolling: Compute remaining output samples */
    blkCnt = blockSize % 0x4U;

#else

    /* Initialize blkCnt with number of taps */
    blkCnt = blockSize;

#endif /* #if defined (ARM_MATH_LOOPUNROLL) */

    while (blkCnt > 0U)
    {
        /* Copy two samples into state buffer */
        *pStateCurnt++ = *pSrc++;

        /* Set the accumulator to zero */
        acc0 = 0;

        /* Use SIMD to hold states and coefficients */
        px = pState;
        pb = pCoeffs;

        tapCnt = numTaps >> 1U;

        while (tapCnt > 0U)
        {
            acc0 += (q31_t)*px++ * *pb++;
            acc0 += (q31_t)*px++ * *pb++;

            tapCnt--;
        }

        /* The result is in 2.30 format. Convert to 1.15 with saturation.
           Then store the output in the destination buffer. */
        *pDst++ = (q15_t)(__SSAT((acc0 >> 15), 16));

        /* Advance state pointer by 1 for the next sample */
        pState = pState + 1U;

        /* Decrement loop counter */
        blkCnt--;
    }

    /* Processing is complete.
       Now copy the last numTaps - 1 samples to the start of the state buffer.
       This prepares the state buffer for the next function call. */

    /* Points to the start of the state buffer */
    pStateCurnt = S->pState;

#if defined(ARM_MATH_LOOPUNROLL)

    /* Loop unrolling: Compute 4 taps at a time */
    tapCnt = (numTaps - 1U) >> 2U;

    /* Copy data */
    while (tapCnt > 0U)
    {
        *pStateCurnt++ = *pState++;
        *pStateCurnt++ = *pState++;
        *pStateCurnt++ = *pState++;
        *pStateCurnt++ = *pState++;

        /* Decrement loop counter */
        tapCnt--;
    }

    /* Calculate remaining number of copies */
    tapCnt = (numTaps - 1U) % 0x4U;

#else

    /* Initialize tapCnt with number of taps */
    tapCnt = (numTaps - 1U);

#endif /* #if defined (ARM_MATH_LOOPUNROLL) */

    /* Copy remaining data */
    while (tapCnt > 0U)
    {
        *pStateCurnt++ = *pState++;

        /* Decrement loop counter */
        tapCnt--;
    }
}
#endif /* defined(ARM_MATH_NEON) */
#endif /* defined(ARM_MATH_MVEI) */

/**
  @} end of FIR group
 */
